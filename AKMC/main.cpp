#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <string>
#include <vector>
#include <time.h>
#include <numeric>
#include <sstream>
#include "cxxopts.hpp"
#include "iofstream.cpp"
#include "../TREE/tree.h"
#include <cassert>
#include <algorithm>

// #define NDEBUG

using namespace std;

void error_exit(string str)
{
	cout<<"ERROR: "<<str<<endl;
	exit(EXIT_FAILURE);
}

void warning(bool condition, string message)
{
	if(condition)
		cout<<"WARNING: "<<message<<endl;
}

class Energy{
	private:
		double a[7], b[5], dE0, c, sdash, kBT, logtau;
	public:
		Energy() : 
		a{9.0e-12, 6.0e-12, 8.0e-12, 2.2e-11, 3.5e-11, 3.0e-12, 1.7e-11},
		b{24.9e-7, 24.3e-7, 5.3e-7, 14.5e-7, 15.5e-7},
		dE0(0.8153), c(0.829), sdash(1.08e4), kBT(0), logtau(0)
		{}

		double calc_energy(double *sss)
		{
			double deltaE=0;
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

		void set_temperature(double temperature, double tau=1.0e-13)
		{
			assert(("temperature must be a positive float", temperature>0));
			assert(("tau must be a positive float", tau>0));
			kBT = 8.617e-5*temperature;
			logtau = log(tau);
		}

		double return_kappa(double *sss)
		{
			assert(("Temperature and tau have not been set", kBT>0 && logtau>0));
			return exp(-calc_energy(sss)/kBT-logtau);
		}
}energy;

class Atom{
	protected:
		Eigen::Vector3d x, *cell;
	public:
		Atom(){
		}
		void set_cell(Eigen::Vector3d *cell_in)
		{
			cell = cell_in;
		}
		void set_position(double *xx)
		{
			x = Eigen::Map<Eigen::Vector3d>(xx, 3);
		}
		void set_position(Eigen::Vector3d xx)
		{
			x = xx;
		}
		Eigen::Vector3d get_position()
		{
			return x;
		}
		Eigen::Vector3d operator-(const Atom& atom)
		{
			auto v = x.array()-atom.x.array();
			return (v-(*cell).array()*(2*v/(*cell).array()).cast<int>().cast<double>()).matrix();
		}
		Eigen::Vector3d operator+(const Atom& atom)
		{
			auto v = 0.5*((*this)-atom).array()+atom.x.array();
			return (v-(*cell).array()*((v/(*cell).array()-0.5).cast<int>().cast<double>())).matrix();
		}
};

class Octa : public Atom{
	private:
		double ** S;
		int * id_NN, count;
		double * x;
	public:
		Octa() : count(0)
		{
			S = new double *[4];
			for(int i=0; i<4; i++)
				S[i] = new double[6];
			id_NN = new int [4];
			x = new double[3];
		}
		void set_neighbor(int id_in, double *S_in){
			assert(("number of neighbors exceeding 4", count<4));
			assert(("problem with the neighbor value", id_in>=0));
			id_NN[count] = id_in;
			for(int i=0; i<6; i++)
				S[count][i] = S_in[i];
			count++;
		}
		~Octa(){
			for(int i=0; i<4; i++)
				delete [] S[i];
			delete [] S;
			delete [] id_NN;
			delete [] x;
		}
};

class CAtom : public Atom{
	private:
		double *kappa;
		Octa * octa;
		int id_octa;
		int zone_id;
	public:
		CAtom(){
			kappa = new double [4];
		}
		void initialize(Octa * octa_in, int id_octa_in)
		{
			octa = octa_in;
			update_id(id_octa_in);
		}
		void update_id(int id_octa_in)
		{
			id_octa = id_octa_in;
			set_position(octa[id_octa].get_position());
		}
};

class AKMC{
	private:
		Octa * octa;
		CAtom * Catom;
		Eigen::Vector3d cell;
		double a_0;
		int N_octa, N_C;
		Eigen::Matrix3d tmat;
		Node *head = new Node;
	public:
		AKMC(double dislocation_density, string input_file, double temperature_in,
			 double box_height, double number_of_C_atoms, vector<double>core_position,
			 double r_core, double r_kmc) : a_0(2.855312531), N_octa(0)
		{
			cout<<"AKMC initialization"<<endl;
			tmat << 1.0/sqrt(3), -1.0/sqrt(2), 1.0/sqrt(6),
					1.0/sqrt(3), 0, -2.0/sqrt(6),
					1.0/sqrt(3), 1.0/sqrt(2),  1.0/sqrt(6);
			energy.set_temperature(temperature_in);
			set_box(box_height, dislocation_density);
			set_number_of_C_atoms(number_of_C_atoms);
			read_octa(input_file, core_position, r_core, r_kmc);
			initialize_C();
		}

		void set_number_of_C_atoms(double nn)
		{
			assert(("illegal number of C atoms", nn>0));
			int N_Fe = int(2.0*cell.prod()/exp(3.0*log(a_0)));
			if(nn<1)
				nn = nn/(1.0-nn)*N_Fe;
			warning(nn>0.1*N_Fe, "Number of C atoms unrealistically high.");
			N_C = int(nn);
		}

		double read_octa(string input_file, vector<double> core_position, double r_core, double r_kmc)
		{
			assert(("core_position has to be x,y-coordinates", core_position.size()!=2));
			warning(abs(core_position.at(0))+abs(core_position.at(1))!=0, "core_position shift not implemented");
			assert(("r_core must be smaller than or equal to r_kmc", r_core<=r_kmc));
			fstream eingabe;
			eingabe.open(input_file, ios::in);
			if(!eingabe)
				error_exit("input file "+input_file+" does not exist or is empty");
			string line;
			int n_in;
			while(eingabe>>n_in)
			{
				if (n_in>=N_octa)
					N_octa = n_in+1;
				getline(eingabe, line);
			}
			cout<<"Number of octahedral interstitial sites: "<<N_octa<<endl;
			eingabe.close();
			octa = new Octa[N_octa];
			eingabe.open(input_file, ios::in);
			double x[3], S[6];
			int id_NN, id_self;
			for(int i=0; i<N_octa && getline(eingabe, line); i++)
			{
				stringstream ss(line);
				ss>>id_self;
				for(int j=0; j<3 && ss>>x[j]; j++);
				ss>>id_NN;
				for(int j=0; j<6 && ss>>S[j]; j++);
				assert(("ID of neighbor is the ID of itself", id_self!=id_NN));
				assert(("neighbor ID larger than number of sites", max(id_self, id_NN)<N_octa));
				assert(("negative ID makes no sense", min(id_NN, id_self)>=0));
				octa[id_self].set_position(x);
				octa[id_self].set_neighbor(id_NN, S);
				octa[id_self].set_cell(&cell);
			}
			eingabe.close();
			cout<<"Octahedral interstitial site configuration import: complete"<<endl;
		}

		void set_box(double box_height, double dislocation_density)
		{
			assert(("Problem with the box height", box_height>0));
			assert(("Dislocation density has to be a positive float", dislocation_density>0));
			cell<<1.0/sqrt(dislocation_density), 1.0/sqrt(dislocation_density), box_height;
			cout<<"Box size: "<<cell.prod()<<endl;
		}

		void initialize_C()
		{
			Catom = new CAtom[N_C];
			for(int i_c=0; i_c<N_C; i_c++)
				Catom[i_c].initialize(octa, rand()%N_octa);
		}

		void run(int number_of_steps=1)
		{
			for(int i_step=0; i_step<number_of_steps; i_step++)
			{
			}
		}
};

int main(int arg, char **name){
	cxxopts::Options options(name[0], "AKMC code");
	int rand_seed=0, nsteps=10000000;
	double temperature=300, dislocation_density=1.0e-6, box_height=27.9762, number_of_C_atoms=1, r_core=10, r_kmc=120;
	vector<double> core_position (2, 0);
	string octa_config="octa.dat", output_file="akmc.log";
	options.add_options()
		("T, temperature", "Temperature K (default: 300 K)", cxxopts::value<double>(temperature))
		("s, seed", "Seed for random number (default: 0)", cxxopts::value<int>(rand_seed))
		("d, density", "Dislocation density in A^2 (default: 10^-6)", cxxopts::value<double>(dislocation_density))
		("n, number", "Number of C atoms or fraction if smaller than 1 (default: 1)", cxxopts::value<double>(number_of_C_atoms))
		("o, output", "Output file name (default: akmc.log)", cxxopts::value<string>(output_file))
		("i, input", "Input file name (default: octa.dat)", cxxopts::value<string>(octa_config))
		("z, height", "Box height (default: 27.9762)", cxxopts::value<double>(box_height))
		("c, core", "Core position in x,y (default: 0, 0)", cxxopts::value<vector<double> >(core_position))
		("r, radius_core", "Radius of the core region in A (default: 10)", cxxopts::value<double>(r_core))
		("R, radius_kmc", "Radius of the KMC region in A (default: 120)", cxxopts::value<double>(r_kmc))
		("h, help", "Print help")
	;
	auto results = options.parse(arg, name);
	srand(rand_seed);
	if (results.count("help"))
	{
		cout<<options.help({"", "Group"})<<endl;
		return 0;
	}
	AKMC akmc = AKMC(dislocation_density, octa_config, temperature, box_height, number_of_C_atoms, core_position, r_core, r_kmc);
}
