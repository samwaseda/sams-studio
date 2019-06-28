#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "sstream"
#include "iostream"
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "corr.h"

#define lll 2.855312531
#define tau 1.0e-13

using namespace std;

double A_base[3][3]={
		{1.0/1.7320508, -1.0/1.4142135, 1.0/2.449489},
		{1.0/1.7320508, 0,          -2.0/2.449489},
		{1.0/1.7320508, 1.0/1.4142135,  1.0/2.449489}
	};
double A(int i, int j){
	if(i<3 && j<3)	return A_base[i][j];
	else {cout<<"ERROR: A["<<i<<"]["<<j<<"] doesn't exist!"<<endl; abort();}
}
double Ortho(double * x, int i){
	return (x[0]*A(i,0)+x[1]*A(i,1)+x[2]*A(i,2))/lll;
}
double lntau=log(tau);
double const critical_r=0.0;
double fq(double x){return x*x;}
double C_concentration, disloc_density,T;
double sc_prod(double* x, double* y){ return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
bool kriterium(double*, double*, double*, double*);
bool finalizer=true;
void normalize(double *x){
	double laenge=sqrt(sc_prod(x,x));
	for(int i=0; i<3; i++) x[i] /= laenge;
}
bool hors_LinCoSS(double* x){
	if(fq(x[0]/40)+fq((x[1])/20)<1) return true;
	else return false;
}
char const* FPN(int, double);
struct chem{
	int ID;
	double dx[3];
};
double calc_energy(double*);
void ID_update(int, int*, int*);
bool chem_en(double*, double *, int, int);

class C {
	private:
		int * Fe_ID, ** c_ID, * c_update_count, c_KMC_ID[2], NC, NFe, Nmax;
		double * dE, * S, dt, dt_accum, ** c_dE, ** c_expdE,  * F_X, * C_X, * c_dist_core, KMC_time, * E_current, * E_modified, boxlo[3], boxhi[3], history[5], r_c_before, r_c, k_BT, Diff_coeff, vf_MCM, kappa_tot;//, S_err, S_ave;
		fstream positionen;
		FILE* output;
		vector<int> Fe_skin, KMC_ID, HL_ID, matrix_ID, c_update, *corr_counter, inserted_ID;
		clock_t begin_time;
		char * file;
		stringstream logfile, dumpfile;
	public:
		double * c_x(int c_id){
			if (c_id>=NC) {cout<<"C ID "<<c_id<<" greater than NC "<<endl; abort();}
			return &(C_X[3*c_id]);}
		double * f_x(int f_id){
			if (f_id>=NFe) {cout<<"Fe ID "<<f_id<<" greater than NFe "<<endl; abort();}
			return &(F_X[3*f_id]);}
		bool core_region(double* x){
			if(fq(x[0])+fq(x[1])<fq(r_c)) return true;
			else return false;
		}
		int natoms(){return NFe+NC;}
		double lx(int i){if(i<3) return boxhi[i]-boxlo[i]; else{cout<<"ERROR: lx"<<endl; abort();}}
		int Fe_index(int ID, int i){if(ID<NC && i<5) return Fe_ID[15*c_ID[ID][1]+5*c_ID[ID][2]+i]; else {cout<<"ERROR: Fe_ID"<<endl; abort();}}
		double cdist(int ID){return sqrt(fq(c_x(ID)[0])+fq(c_x(ID)[1]));}
		double fdist(int ID){return sqrt(fq(f_x(ID)[0])+fq(f_x(ID)[1]));}
		int core_count(){int count=HL_ID.size()+inserted_ID.size(); for(int i=0; i<KMC_ID.size(); i++) count += core_region(c_x(KMC_ID.at(i))); if(count>=Nmax){finalizer=false;} return count;}
		int interactable_size(){return inserted_ID.size()+KMC_ID.size();}
		int interactable(int ID){if(ID<KMC_ID.size()){return KMC_ID.at(ID);}else{return inserted_ID.at(ID-KMC_ID.size());}}
		void vec_update(){
			HL_ID.clear();
			KMC_ID.clear();
			matrix_ID.clear();
			for(int i=0; i<NC; i++)
			{
				if(c_ID[i][3]==0) HL_ID.push_back(i);
				if(c_ID[i][3]==1) KMC_ID.push_back(i);
				if(c_ID[i][3]==2) matrix_ID.push_back(i);
			}
			while(KMC_ID.size()==0){dt_accum=1.0/(4.0*Diff_coeff); outside_mapping();}
		}
		double dist(double * xx, double * yy){
			double laenge=fq(xx[2]-yy[2]-lx(2)*int(1.5*(xx[2]-yy[2])/lx(2)));
			for(int i=0; i<2; i++) laenge += fq(xx[i]-yy[i]);
			//for(int i=0; i<3; i++) laenge += fq(xx[i]-yy[i]-lx(i)*int(1.5*(xx[i]-yy[i])/lx(i)));
			return sqrt(laenge);
		}
		void dx_PBC (double* xx, double* yy, double* zz){
			zz[2]=xx[2]-yy[2]-lx(2)*int(1.5*(xx[2]-yy[2])/lx(2));
			for(int i=0; i<2; i++) zz[i]=xx[i]-yy[i];
			//for(int i=0; i<3; i++){
				//zz[i]=xx[i]-yy[i]-lx(i)*int(1.5*(xx[i]-yy[i])/lx(i));
				//if(abs(zz[i])>lx(i)){cout<<"ERROR ZZ"<<endl; abort();}
			//}
		}
		void x_mid_PBC (double* xx, double* yy, double* zz){
			dx_PBC(xx, yy, zz);
			zz[2] = yy[2]+0.5*zz[2]-lx(2)*floor((yy[2]+0.5*zz[2]-boxlo[2])/lx(2));
			for(int i=0; i<2; i++) zz[i] = yy[i]+0.5*zz[i];
		}
		int rel_type(int ID_self, int ID_other){ return (c_ID[ID_self][2]-c_ID[ID_other][2]+3)%3;}
		C(char * file_in, double r_c_b_in, double r_c_in, int Nmax_in) : dt_accum(0), KMC_time(0), file(file_in), r_c_before(r_c_b_in), r_c(r_c_in), Nmax(Nmax_in), k_BT(T*8.617e-5)
		{
			Diff_coeff=5.0e13*exp(-0.8153/k_BT);
			vf_MCM=sqrt(4.0*Diff_coeff);
			cout<<"Starting initialization"<<endl;
			int ID_in, ID_type_in, ID_type_count=0, Fe_ID_in, Fe_orient_ID_in, x_in[3], KMCcount=0, corecount=0;
			double s_in;
			dumpfile<<file<<".xyz";
			positionen.open(dumpfile.str().c_str(), ios::out);
			logfile<<file<<".dat";
			output = fopen(logfile.str().c_str(), "w");
			string line;
			fstream eingabe;
			eingabe.open("/home_nfs/mateis/owaseda/PHD/AKMC/init.dat", ios::in);
			if(!eingabe)
				{cout<<"ERROR: init.dat not found "<<endl; abort();}
			else
				cout<<"init.dat opened"<<endl;
			for(int i=0; i<3; i++)
				getline(eingabe,line);
			eingabe>>NFe;
			cout<<NFe<<" Fe atoms"<<endl;
			for(int i=0; i<2; i++)
				getline(eingabe,line);
			for(int i=0; i<3; i++)
				eingabe>>boxlo[i]>>boxhi[i];
			for(int i=0; i<2; i++)
				getline(eingabe,line);
			boxlo[0]=-sqrt(disloc_density)*0.5;
			boxlo[1]=-sqrt(disloc_density)*0.5;
			boxhi[0]=sqrt(disloc_density)*0.5;
			boxhi[1]=sqrt(disloc_density)*0.5;
			Fe_ID = new int [15*NFe];
			S = new double [6*NFe];
			F_X = new double[NFe*3];
			cout<<"Reading in Fe information"<<endl;
			for(int i=0; i<NFe; i++)
			{
				eingabe>>Fe_ID_in;
				if(Fe_ID_in<0 || Fe_ID_in>=NFe)
					{cout<<"ERROR with Fe ID "<<Fe_ID_in<<endl; abort();}
				for(int j=0; j<3; j++)
					eingabe>>F_X[3*Fe_ID_in+j];
				for(int j=0; j<6; j++)
					eingabe>>S[6*Fe_ID_in+j];
				for(int j=0; j<15; j++)
					eingabe>>Fe_ID[15*Fe_ID_in+j];
			}
			for(int i=0; i<NFe; i++)
				if(fdist(i)>145.0 && fdist(i)<150.0)
					Fe_skin.push_back(i);
			eingabe.close();
			NC=(int)(C_concentration*2.0*(lx(0)/lll)*(lx(1)/lll)*(lx(2)/lll));
			cout<<"Number of carbon atoms: "<<NC<<endl;
			C_X = new double [3*NC];
			c_ID = new int* [NC];			// ID of C atoms, ID of Fe attributed & vector
			for(int i=0; i<NC; i++)			// 0: ID 1: Fe_ID 2: Fe_ID_orient
				c_ID[i] = new int [4];		// 3: core/KMC/matrix/inserted
			corr_counter = new vector<int> [NC];
			eingabe.open("C_init.dat", ios::in);
			if(eingabe)
			{
				cout<<"starting to read C_init.dat"<<endl;
				while(eingabe>>ID_in)
				{
					c_ID[ID_in][0]=ID_in;
					eingabe>>c_x(ID_in)[0]>>c_x(ID_in)[1]>>c_x(ID_in)[2]>>c_ID[ID_in][1]>>c_ID[ID_in][2];
					if(cdist(ID_in)<r_c_before)
					{
						c_ID[ID_in][3]=3;
						corecount++;
						inserted_ID.push_back(ID_in);
					}
					else if(cdist(ID_in)<150)
					{
						c_ID[ID_in][3]=1;
						KMCcount++;
					}
					else
						c_ID[ID_in][3]=2;
				}
			}
			else
			{
				cout<<"Creating new C configuration"<<endl;
				for(int i=0; i<NC; i++)
				{
					c_ID[i][0]=i;
					for(int j=0; j<3; j++)
						c_x(i)[j] = rand()/(double)RAND_MAX*lx(j)+boxlo[j];
					if(cdist(i)<150)
					{
						for(bool criterion=true; criterion; )
						{
							criterion=false;
							c_ID[i][1]=rand()%NFe;
							c_ID[i][2]=rand()%3;
							for(int j=0; j<3; j++)
								c_x(i)[j]=f_x(c_ID[i][1])[j]+0.5*lll*A(c_ID[i][2],j);
							if(hors_LinCoSS(c_x(i)))
								c_ID[i][3]=0;
							else if(cdist(i)<150)
								c_ID[i][3]=1;
							else
								c_ID[i][3]=2;
							for(int j=0; j<i; j++)
								if(dist(c_x(i),c_x(j))<5.0)
									criterion=true;
							if(!criterion)
							{
								if(hors_LinCoSS(c_x(i)))
									corecount++;
								if(cdist(i)<150)
									KMCcount++;
							}
						}
					}
					else
						c_ID[i][3]=2;
				}
			}
			cout<<KMCcount<<" C atoms inside KMC zone. "<<corecount<<" atoms inside core"<<endl;
			c_dE = new double *[NC];		// Energies associated with jumps. 4 values for each C
			c_expdE = new double *[NC];		// Energies associated with jumps. 4 values for each C
			c_update_count = new int [NC];
			for(int i=0; i<NC; i++)
				c_update_count[i]=0;
			for(int i=0; i<NC; i++)
			{
				c_dE[i] = new double [4];
				c_expdE[i] = new double [4];
				if(c_ID[i][3]==1)
					c_update.push_back(i);
			}
			for(int i=0; i<NC; i++)
				for(int j=0; j<4; j++)
					c_expdE[i][j]=0;
			eingabe.close();
			vec_update();
			begin_time = clock();
		}
		~C(){
			positionen.close();
			dumpfile.str("");
			dumpfile<<file<<"_end.xyz";
			positionen.open(dumpfile.str().c_str(), ios::out);
			KMC_output(0);
			cout<<"#Finishing KMC"<<endl;
		}

		void kappa_tot_actualization(){
			kappa_tot=0;

			for(int i=0; i<KMC_ID.size(); i++)
				for(int j=0; j<4; j++)
				{
					c_expdE[KMC_ID.at(i)][j]=exp(-c_dE[KMC_ID.at(i)][j]/k_BT-lntau);
					kappa_tot+=c_expdE[KMC_ID.at(i)][j];
				}
		}

		bool KMC(int MCStep){ // FUNCTION
			for(int i=0; i<c_update.size(); i++)
			{
				if(c_ID[c_update.at(i)][3]==1)
					LinCoSS(c_update.at(i));
				else
					{cout<<"ERROR c_dE not defined for atom "<<c_update.at(i)<<" of Fe "<<c_ID[c_update.at(i)][1]<<endl; abort();}
				for(int j=0; j<4; j++)
				{
					kappa_tot -= c_expdE[c_update.at(i)][j];
					c_expdE[c_update.at(i)][j] = exp(-c_dE[c_update.at(i)][j]/k_BT-lntau);
					kappa_tot += c_expdE[c_update.at(i)][j];
				}
			}
			c_update.clear();
			dt = KMC_choice();

			if(c_dE[c_KMC_ID[0]][c_KMC_ID[1]]>10)
				{cout<<"ERROR choosing dE for atom "<<c_KMC_ID[0]<<endl; abort();}
			KMC_time += dt;
			dt_accum += dt;


			c_update_count[c_KMC_ID[0]]++;
			ID_update(c_KMC_ID[1], c_ID[c_KMC_ID[0]], &(Fe_ID[15*c_ID[c_KMC_ID[0]][1]+5*c_ID[c_KMC_ID[0]][2]]));
			history[MCStep%5]=c_KMC_ID[0];

			int output_frequency=1000;
			if(MCStep%output_frequency==0)
			{
				fprintf(output, "%d s %.2f",MCStep,KMC_time);
				fprintf(output, " s/MC %.3e",double(clock()-begin_time)/CLOCKS_PER_SEC/MCStep);
				fprintf(output, " MCs/Rs %.3e",KMC_time/(double(clock()-begin_time)/CLOCKS_PER_SEC));
				fprintf(output, " c_ID%s%d",FPN(6,c_ID[c_KMC_ID[0]][0]),c_ID[c_KMC_ID[0]][0]);
				//fprintf(output, " c_ID%s%d Dist%s%.2f",string(FPN(7,c_ID[c_KMC_ID[0]][0]),' ').c_str(),c_ID[c_KMC_ID[0]][0],string(FPN(4,cdist(c_KMC_ID[0])),' ').c_str(),cdist(c_KMC_ID[0]));
				//fprintf(output, " deg%s%.2f",string(FPN(4,atan2(c_x(c_KMC_ID[0])[0],c_x(c_KMC_ID[0])[1])*180/M_PI),' ').c_str(),atan2(c_x(c_KMC_ID[0])[0],c_x(c_KMC_ID[0])[1])*180/M_PI);
				fprintf(output, " UD %s%d dE %.4f",FPN(3,c_update_count[c_KMC_ID[0]]),c_update_count[c_KMC_ID[0]],c_dE[c_KMC_ID[0]][c_KMC_ID[1]]);
				//fprintf(output, " S_ave%s%.4f dS %.4f",string(FPN(5,S_ave),' ').c_str(),S_ave,S_err);
				fprintf(output, " NKMC %d",(int)KMC_ID.size());
				fprintf(output, " NC %d",core_count());
				fprintf(output, " HL %d",(int)HL_ID.size());
				fprintf(output, " NGH %d",(int)corr_counter[c_KMC_ID[0]].size());
				fprintf(output, " x_o%s%.3f%s%.3f%s%.3f",FPN(4,c_x(c_KMC_ID[0])[0]),c_x(c_KMC_ID[0])[0],FPN(4,c_x(c_KMC_ID[0])[1]),c_x(c_KMC_ID[0])[1],FPN(4,c_x(c_KMC_ID[0])[2]),c_x(c_KMC_ID[0])[2]);
			}
			x_mid_PBC(f_x(Fe_index(c_KMC_ID[0],4)), f_x(c_ID[c_KMC_ID[0]][1]), c_x(c_KMC_ID[0]));
			if(MCStep%output_frequency==0)
			{
				fprintf(output, " x_n%s%.3f%s%.3f%s%.3f",FPN(4,c_x(c_KMC_ID[0])[0]),c_x(c_KMC_ID[0])[0],FPN(4,c_x(c_KMC_ID[0])[1]),c_x(c_KMC_ID[0])[1],FPN(4,c_x(c_KMC_ID[0])[2]),c_x(c_KMC_ID[0])[2]);
				fprintf(output, " st %d", c_ID[c_KMC_ID[0]][2]);
				fprintf(output, " dE %.4f %.4f %.4f %.4f",c_dE[c_KMC_ID[0]][0],c_dE[c_KMC_ID[0]][1],c_dE[c_KMC_ID[0]][2],c_dE[c_KMC_ID[0]][3]);
				fprintf(output, "\n");
				fflush(output);

			}
			if(c_update_count[c_KMC_ID[0]]==5000 || count(history, history+5, c_KMC_ID[0])==5)
			{
				cout<<"Irregularity with C ID "<<c_KMC_ID[0]<<" at "<<MCStep<<" putting back to somewhere in camembert"<<endl;
				for(int ii=0; ii<interactable_size(); ii++)
				{
					int i=interactable(ii);
					if(c_KMC_ID[0]!=i && dist(c_x(i), c_x(c_KMC_ID[0]))<10.0 && c_ID[i][3]%2==1)
					{
						double dx_temp[3];
						dx_PBC(c_x(c_KMC_ID[0]), c_x(i), dx_temp);
						//cout<<c_ID[c_KMC_ID[0]][2]<<" "<<i<<" "<<c_ID[i][2]<<" "<<dx_temp[0]<<" "<<dx_temp[1]<<" "<<dx_temp[2]<<endl;
					}
				}
				for(bool criterion=true; criterion; )
				{
					criterion=false;
					c_ID[c_KMC_ID[0]][1]=Fe_skin.at(rand()%Fe_skin.size());
					c_ID[c_KMC_ID[0]][2]=rand()%3;
					for(int j=0; j<3; j++)
						c_x(c_KMC_ID[0])[j] = f_x(c_ID[c_KMC_ID[0]][1])[j]+0.5*lll*A(c_ID[c_KMC_ID[0]][2],j);
					if(cdist(c_KMC_ID[0])>150)
						criterion=true;
					for(int j=0; j<KMC_ID.size() && !criterion; j++)
						if(dist(c_x(c_KMC_ID[0]),c_x(KMC_ID.at(j)))<10.0 && KMC_ID.at(j)!=c_KMC_ID[0])
							criterion=true;
				}
				c_ID[c_KMC_ID[0]][3]=1;
				c_update_count[c_KMC_ID[0]]=0;
				c_update.push_back(c_KMC_ID[0]);
				vec_update();
			}

			if(cdist(c_KMC_ID[0])>150)
			{
				//cout<<"Atom "<<c_KMC_ID[0]<<" at "<<c_x(c_KMC_ID[0])[0]<<" "<<c_x(c_KMC_ID[0])[1]<<" "<<c_x(c_KMC_ID[0])[2]<<" got out of cylinder"<<endl;;
				c_ID[c_KMC_ID[0]][3]=2;
				c_update_count[c_KMC_ID[0]]=0;
				vec_update();
				kappa_tot_actualization();
			}
			else if(hors_LinCoSS(c_x(c_KMC_ID[0])))
			{
				//cout<<"Atom reached core for ID "<<c_KMC_ID[0]<<" at time step "<<MCStep<<" at real time "<<KMC_time<<"s"<<endl;
				c_ID[c_KMC_ID[0]][3]=0;
				vec_update();
				kappa_tot_actualization();
			}
			else
				c_update.push_back(c_KMC_ID[0]);
			for(int i=0; i<KMC_ID.size(); i++)
				if(dist(c_x(KMC_ID.at(i)),c_x(c_KMC_ID[0]))<10 && c_ID[KMC_ID.at(i)][3]==1 && KMC_ID.at(i)!=c_KMC_ID[0])
					c_update.push_back(KMC_ID.at(i));

			if(dt_accum>1.0/(4.0*Diff_coeff))
				if(!outside_mapping())
					return false;
			return true;
		}

		bool outside_mapping()
		{
			dt = vf_MCM*sqrt(dt_accum);
			for(int i=0; i<NC; i++)
				if(c_ID[i][3]==2)
				{
					c_x(i)[0] += dt*cos(rand());
					c_x(i)[1] += dt*sin(rand());
					c_x(i)[0] += -lx(0)*floor((c_x(i)[0]-boxlo[0])/lx(0));
					c_x(i)[1] += -lx(1)*floor((c_x(i)[1]-boxlo[1])/lx(1));
					if(cdist(i)<147.5)
					{
						//cout<<"Atom "<<c_ID[i][0]<<" has reached the cylinder and will be put around Fe ";
						for(bool criterion=true; criterion; )
						{
							criterion=false;
							c_ID[i][1]=Fe_skin.at(rand()%Fe_skin.size());
							c_ID[i][2]=rand()%3;
							for(int j=0; j<3; j++)
								c_x(i)[j] = f_x(c_ID[i][1])[j]+0.5*lll*A(c_ID[i][2],j);
							c_ID[i][3]=1;
							if(cdist(i)>150)
								criterion=true;
							for(int j=0; j<NC && !criterion; j++)
								if(dist(c_x(i),c_x(j))<5.0 && j!=i)
									criterion=true;
						}
						//cout<<c_ID[i][1]<<" at "<<f_x(c_ID[i][1])[0]<<" "<<f_x(c_ID[i][1])[1]<<" "<<f_x(c_ID[i][1])[2]<<endl;
						vec_update();
						c_update.push_back(i);
					}
				}
			dt_accum=0;
			kappa_tot_actualization();
			return true;
		}
		bool LinCoSS(int ID)
		{
			double xn[18], Sn[24];	//stress values of neighboring atoms.
			//double S_err_temp = 0, S_ave_temp=0;

			for(int i=0; i<3; i++)
				xn[i]=f_x(c_ID[ID][1])[i];
			for(int i=0; i<5; i++)
				for(int j=0; j<3; j++)
					xn[3*(i+1)+j]=f_x(Fe_index(ID,i))[j];
			for(int i=0; i<4; i++)
				for(int j=0; j<6; j++)
					Sn[6*i+j] = 0.25*(S[6*Fe_index(ID,i)+j]+S[6*Fe_index(ID,(i+3)%4)+j]);

			for(int i=0; i<4; i++)
				for(int j=0; j<6; j++)
					Sn[6*i+j] += 0.25*(S[6*c_ID[ID][1]+j]+S[6*Fe_index(ID,4)+j]);

			trafo(Sn, xn, ID);

			return true;
		}
		void trafo(double * SS, double* xn, int ID) // FUNCTION
		{
			double ex[9], s[6], x_begin[3], x_end[3], s_corr[36], dx[3], dE_corr[4];
			int ix[3], ixx[3], vz[3];
			corr_counter[ID].clear();

			for(int i=0; i<36; i++)
				s_corr[i]=0;
			for(int i=0; i<4; i++)
				dE_corr[i] = 0;
			for(int ii=0; ii<interactable_size(); ii++)
			{
				int i=interactable(ii);
				if(ID!=i && dist(c_x(i), c_x(ID))<10.0 && c_ID[i][3]%2==1)
				{
					dx_PBC(c_x(ID), c_x(i), dx);
					corr_counter[ID].push_back(i);
					if(!chem_en(dE_corr, dx, c_ID[i][2], c_ID[ID][2])){	cout<<ID<<"\\|"<<i<<endl; abort();}
					for(int j=0; j<3; j++)
						ix[(j+2-c_ID[i][2])%3] = (1-2*(dx[j]<0))*int(abs(dx[j]/lll*2.0)+0.5);
					for(int j=0; j<4; j++)
						for(int k=0; k<6; k++)
							s_corr[6*(1+(j%2)+(2*(j>1)-1)*(1-2*(ix[j%2]<0)))+k] += st_corr[abs(ix[0])][abs(ix[1])][abs(ix[2])][j][k];
				}
			}
			dx_PBC(&(xn[15]), &(xn[0]), ex);
			normalize(ex);
			x_mid_PBC(&(xn[15]), &(xn[0]), x_begin);

			for(int i=0; i<4; i++)
			{
				dx_PBC(&(xn[3*(i+1)]), &(xn[3*((i+3)%4+1)]), &(ex[3]));
				normalize(&(ex[3]));
				x_mid_PBC(&(xn[3*(i+1)]), &(xn[3*((i+3)%4+1)]), x_end);
				dx_PBC(x_end, x_begin, &(ex[6]));
				normalize(&(ex[6]));
				for(int j=0; j<6; j++)
					s[j]=0;
				for(int j=0; j<3; j++)
					for(int k=0; k<=j; k++)
						for(int l=0; l<3; l++)
							for(int m=0; m<3; m++)
								s[j+(k+2)*(j!=k)]+=ex[3*k+m]*ex[3*j+l]*SS[6*i+m+(l+2)*(m!=l)];
				//if(s_corr[0]+s_corr[1]+s_corr[2]!=0)
					//cout<<"Plain: "<<calc_energy(s);
				for(int j=0; j<6; j++)
					s[j] += s_corr[6*i+j];
				c_dE[ID][i]=calc_energy(s);
				//if(s_corr[0]+s_corr[1]+s_corr[2]!=0)
					//cout<<" Elastic correction: "<<c_dE[ID][i];
				c_dE[ID][i] += dE_corr[i];
				//if(s_corr[0]+s_corr[1]+s_corr[2]!=0)
					//cout<<" Chemical correction: "<<c_dE[ID][i]<<endl;
				if(c_dE[ID][i]<=-1 || c_dE[ID][i]>500)
				{
					cout<<"Error with energy "<<c_dE[ID][i]<<" for "<<ID<<" "<<i<<" "<<c_ID[ID][1]<<" "<<Fe_index(ID,0)<<" "<<corr_counter[ID].size()<<endl;
					//cout<<"Error with energy "<<c_dE[ID][i]<<" for "<<ID<<" "<<i<<" "<<c_ID[ID][1]<<" "<<Fe_index(ID,0)<<" "<<S_ave<<" "<<S_err<<" "<<corr_counter<<endl;
					for(int j=0; j<3; j++)
						cout<<" "<<c_x(ID)[j];
					cout<<endl;
					for(int j=0; j<6; j++)
					{
						for(int k=0; k<6; k++)
							cout<<" "<<s_corr[6*j+k];
						cout<<endl;
					}
					cout<<endl;
					abort();
				}
				//if(c_ID[ID][4]==2)
				//{
					//for(int j=0; j<3; j++)
						//ex[j+3]=x_begin[j]+0.5*lll*ex[j];
					//if(core_region(&(ex[3])))
						//c_dE[ID][i]=100;
					//cout<<"Preventing the C "<<ID<<" from entering "<<ex[3]<<" "<<ex[4]<<" "<<ex[5]<<endl;
				//}
			}
		}

		double KMC_choice()
		{

			/*****************MONTE CARLO TIME**********************/
			double r=rand()/(double)RAND_MAX;
			double deltat;
			if(kappa_tot>0)
				deltat=-log(r)/kappa_tot;
			else
			{
				kappa_tot_actualization();
				if(kappa_tot>0){cout<<"ERROR kappa_tot<=0 "<<kappa_tot<<endl; abort();}
			}

			/*****************MONTE CARLO ENERGY********************/
			double s=rand()/(double)RAND_MAX;
			s*=kappa_tot*0.5;
			while(s>=0)
			{
				c_KMC_ID[0]=KMC_ID.at(rand()%KMC_ID.size());
				c_KMC_ID[1]=rand()%4;
				s-=c_expdE[c_KMC_ID[0]][c_KMC_ID[1]];
				if(s<=0)
					return deltat;
			}
			
			cout<<"Error: KMC atom could not be chosen"<<endl;
			return 0;
		}
		void KMC_output(int MCStep)
		{
			positionen<<"ITEM: TIMESTEP\n"<<MCStep<<"\nITEM: NUMBER OF ATOMS\n"<<NC<<"\nITEM: BOX BOUNDS pp pp pp"<<endl;
			for(int i=0; i<3; i++)
				positionen<<boxlo[i]<<" "<<boxhi[i]<<endl;
			positionen<<"ITEM: ATOMS id type core x y z"<<endl;
			for(int i=0; i<NC; i++)
				positionen<<i+1<<" "<<c_ID[i][3]<<" "<<core_region(c_x(i))<<" "<<c_x(i)[0]<<" "<<c_x(i)[1]<<" "<<c_x(i)[2]<<endl;
		}
};

int main(int narg, char **arg)
{
	if(narg<8 || narg>9)
	{
		cout<<"Syntax: executable filename r_c_before r_c Nmax concentration disloc_density (seed)"<<endl;
		abort();
	}
	if(narg==9)
		srand(time(NULL)*time(NULL)*(atoi(arg[8])+1)*(atoi(arg[8])+1));
	else
		srand(time(NULL)*time(NULL));
	chem corr;
	cout<<rand()<<endl;
	C_concentration=atof(arg[5]);
	disloc_density=atof(arg[6]);
	T=atof(arg[7]);

	C c(arg[1], atof(arg[2]), atof(arg[3]), atoi(arg[4]));

	for(int MCStep=0; finalizer; MCStep++)
	{
		//cout<<"                                                 "<<endl;
		//cout<<"  88      a8P  88b           d88    ,ad8888ba,   "<<endl;
		//cout<<"  88    ,88'   888b         d888   d8''    `'8b  "<<endl;
		//cout<<"  88  ,88'     88`8b       d8'88  d8'            "<<endl;
		//cout<<"  88,d88'      88 `8b     d8' 88  88             "<<endl;
		//cout<<"  8888'88,     88  `8b   d8'  88  88             "<<endl;
		//cout<<"  88P   Y8b    88   `8b d8'   88  Y8,            "<<endl;
		//cout<<"  88     '88,  88    `888'    88   Y8a.    .a8P  "<<endl;
		//cout<<"  88       Y8  88     `8'     88    `'Y8888Y''   "<<endl;
		//cout<<"                                                 "<<endl;
		c.KMC(MCStep);
		if(MCStep%1000000==0)
			c.KMC_output(MCStep);
	}

	return 0;
}

double calc_energy(double *sss)
{
	double deltaE;
	//double a[]= {1.4e-11, 9.4e-13, 7.7e-12, 2.3e-11, 3.6e-11, 3.3e-12, 1.77e-11};
	//double b[]= {23.5e-7, 23.9e-7, 5.7e-7, 1.5e-6, 1.6e-6};
	double a[] = {9.0e-12, 6.0e-12, 8.0e-12, 2.2e-11, 3.5e-11, 3.0e-12, 1.7e-11};
	double b[] = {24.9e-7, 24.3e-7, 5.3e-7, 14.5e-7, 15.5e-7};
	double dE0= 0.8153;
	double c=0.829;
	double sdash=1.08e4;
	deltaE =-a[0]*sss[0]*sss[0]+b[0]*sss[0];
	deltaE+= a[1]*sss[1]*sss[1]-b[1]*sss[1];
	deltaE+=-a[2]*sss[2]*sss[2]-b[2]*sss[2];
	deltaE+=dE0;
	if(sss[4]<-sdash)
		deltaE+=-a[5]*sss[4]*sss[4]+b[4]*sss[4]+c-dE0;
	else if(sss[4]>sdash)
		deltaE+=-a[5]*sss[4]*sss[4]-b[4]*sss[4]+c-dE0;
	else
		deltaE+=-a[6]*sss[4]*sss[4];
	if(sss[5]<-sdash)
		deltaE+=-a[3]*sss[5]*sss[5]+b[3]*sss[5]+c-dE0;
	else if(sss[5]>sdash)
		deltaE+=-a[3]*sss[5]*sss[5]-b[3]*sss[5]+c-dE0;
	else
		deltaE+=-a[4]*sss[5]*sss[5];

	return deltaE;
}

void ID_update(int jump_ID, int* c_ID, int* Fe_ID)
{
	if(jump_ID==0)
	{
		c_ID[1] = Fe_ID[3];
		c_ID[2] = (c_ID[2]+1)%3;
	}
	else if(jump_ID==1)
	{
		c_ID[1] = Fe_ID[1];
		c_ID[2] = (c_ID[2]+2)%3;
	}
	else if(jump_ID==2)
	{
		c_ID[1] = Fe_ID[2];
		c_ID[2] = (c_ID[2]+1)%3;
	}
	else if(jump_ID==3)
	{
		c_ID[1] = Fe_ID[2];
		c_ID[2] = (c_ID[2]+2)%3;
	}
}

char const* FPN(int n, double x)
{
	int result;
		result=n-int(log(abs(x)+(abs(x)<1))/log(10.0))-(x<0);
	return strdup(string(result*(0<result),32).c_str());
}

bool chem_en(double* dE, double * dx, int ID, int typ)
{
	int n[3], N[3], vz[3];
	double xx[3];
	for(int i=0; i<3; i++)
	{
		xx[i] = dx[0]*A(i,0)+dx[1]*A(i,1)+dx[2]*A(i,2);
		n[(i+2-ID)%3] = int(abs(2.0*xx[i]/lll)+0.5);
		vz[(i+2-ID)%3] = 1-2*(2.0*xx[i]/lll<-0.25);
	}
	for(int i=0; i<3; i++)
		dx[i]=xx[i];
	for(int i=0; i<3; i++)
		N[i]=n[i];
	N[2]++;
	N[(typ-ID+2)%3]++;
	if((N[0]%2+N[1]%2+N[2]%2)%3!=0)		// Assuming distance between two Fe atoms being i*ex+j*ey+k*ez+l*(ex+ey+ez)/2
		{cout<<"Error: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<xx[0]<<" "<<xx[1]<<" "<<xx[2]<<" "<<typ<<" "<<ID<<endl; return false;}

	for(int i=0; i<4; i++)
	{
		//cout<<dE[i]<<"\t"<<chem_E[n[0]][n[1]][n[2]][1+i%2+(2*(i>1)-1)*vz[(ID+1+i%2)%3]]<<endl;
		//dE[i] += chem_E[n[0]][n[1]][n[2]][1+i%2+(2*(i>1)-1)*vz[(ID+1+i%2)%3]]+chem_stop_E[n[0]][n[1]][n[2]][1+i%2+(2*(i>1)-1)*vz[(ID+1+i%2)%3]];
		for(int j=0; j<3; j++)
			N[j] = n[j]+vz[j]*(((i==0)-(i==2))*(j==2)+((i==1)-(i==3))*(j==1));
		if(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]<=8 && N[0]*N[0]+N[1]*N[1]+N[2]*N[2]<n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
			dE[i] = 9;
		else
			dE[i] += chem_E[n[0]][n[1]][n[2]][1+i%2+(2*(i>1)-1)*vz[(ID+1+i%2)%3]];
	}
	//if(dE[0]>=100 && dE[1]>=100 && dE[2]>=100 && dE[3]>=100)
		//{cout<<"Error EN: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<xx[0]<<" "<<xx[1]<<" "<<xx[2]<<" "<<typ<<" "<<ID<<" "<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl; return false;}
	return true;
}

