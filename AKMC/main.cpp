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

using namespace std;

void error_exit(string str)
{
    cout<<str<<endl;
    exit(EXIT_FAILURE);
}

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

class FeAtom : public Atom{
    private:
        Eigen::Matrix3d S;
    public:
        FeAtom(){
        }
        void set_stress(double *SS)
        {
            S << SS[0], SS[3], SS[4],
                 SS[3], SS[1], SS[5],
                 SS[4], SS[5], SS[2];
        }
};

class CAtom : public Atom{
    private:
        double * E;
    public:
        CAtom(){
            E = new double [4];
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
        void set_position(double *xx){
            for(int i=0; i<3; i++)
                x[i] = xx[i];
        }
        void set_neighbor(int id_in, double *S_in){
            if(count==4)
                error_exit("ERROR: number of neighbors exceeding 4");
            if(id_in<0)
                error_exit("ERROR: problem with the neighbor value");
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

class AKMC{
    private:
        FeAtom * Featom;
        CAtom * Catom;
        Eigen::Vector3d cell;
        double temperature, a_0;
        int N_Fe, N_C;
        Eigen::Matrix3d tmat;
    public:
        AKMC(double dislocation_density, string input_file, double temperature_in, double box_height) : a_0(2.855312531)
        {
            cout<<"AKMC initialization"<<endl;
            tmat << 1.0/1.7320508, -1.0/1.4142135, 1.0/2.449489,
		            1.0/1.7320508, 0, -2.0/2.449489,
		            1.0/1.7320508, 1.0/1.4142135,  1.0/2.449489;
            temperature = temperature_in;
            if (temperature<=0)
                error_exit("ERROR: temperature must be a positive float");
            set_box(box_height, dislocation_density);
        }

        double read_octa(string input_file)
        {
            fstream eingabe;
            eingabe.open(input_file, ios::in);
            if(!eingabe)
                error_exit("ERROR: input file "+input_file+" does not exist or is empty");
            double box_height=0;
            for(int i=0; i<3; i++)
                eingabe>>box_height;
            if (box_height<=0)
                error_exit("ERROR: problem with the box height: "+to_string(box_height));
            string line;
            for(N_Fe=-1; getline(eingabe, line); N_Fe++);
            cout<<"Number of Fe atoms: "<<N_Fe<<endl;
            eingabe.close();
            Featom = new FeAtom[N_Fe];
            eingabe.open(input_file, ios::in);
            for(int i=0; i<3; i++)
                eingabe>>box_height;
            double x[3], S[6];
            for(int i=0; i<N_Fe && getline(eingabe, line); i++)
            {
                if (line.length()==0){
                    i--;
                    continue;
                }
                stringstream ss(line);
                for(int j=0; j<3 && ss>>x[j]; j++);
                for(int j=0; j<6 && ss>>S[j]; j++);
                Featom[i].set_position(x);
                Featom[i].set_stress(S);
                Featom[i].set_cell(&cell);
            }
            eingabe.close();
            cout<<"Fe configuration import: complete"<<endl;
            return box_height;
        }

        void set_box(double box_height, double dislocation_density)
        {
            if(dislocation_density<=0)
                error_exit("ERROR: Dislocation density has to be a positive float");
            cell<<1.0/sqrt(dislocation_density), 1.0/sqrt(dislocation_density), box_height;
            cout<<"Box size: "<<cell.prod()<<endl;
        }

        void detect_sites()
        {
            clock_t begin = clock();
            int NN_ID[3], NNN_ID[3][4], count_site=0, shift_ID=10;
            vector<int> count_NN, count_NNN, NN_ID_tot;
            double dist;
            Eigen::Vector3d v;
            for(int i_Fe=0; i_Fe<N_Fe; i_Fe++)
            {
                count_NN.assign(3, 0); count_NNN.assign(3, 0);
                for(int j_Fe=((i_Fe-shift_ID+N_Fe)%N_Fe); j_Fe!=((i_Fe-shift_ID+N_Fe-1)%N_Fe); j_Fe=(j_Fe+1)%N_Fe)
                {
                    if (i_Fe==j_Fe)
                        continue;
                    dist = (Featom[j_Fe]-Featom[i_Fe]).norm();
                    if(dist>1.2*a_0)
                        continue;
                    for(int ix=0; ix<3; ix++)
                    {
                        if (count_NN.at(ix)==1 && count_NNN.at(ix)==4)
                            continue;
                        v << 0.5*a_0*(ix==0), 0.5*a_0*(ix==1), 0.5*a_0*(ix==2);
                        v = tmat.transpose()*v;
                        if(((Featom[j_Fe]-Featom[i_Fe])+v).norm()<0.97*a_0)
                        {
                            if (count_NN.at(ix)==0 && ((Featom[j_Fe]-Featom[i_Fe])+2.0*v).norm()<0.05*a_0)
                            {
                                NN_ID[ix] = j_Fe;
                                count_NN.at(ix)++;
                            }
                            else
                            {
                                NNN_ID[ix][count_NNN.at(ix)] = j_Fe;
                                count_NNN.at(ix)++;
                            }
                        }
                    }
                }
                for(int ix=0; ix<3; ix++)
                {
                    if(count_NN.at(ix)!=1 || count_NNN.at(ix)!=4)
                        continue;
                    NN_ID_tot.push_back(i_Fe);
                    NN_ID_tot.push_back(NN_ID[ix]);
                    for(int count=0; count<4; count++)
                        NN_ID_tot.push_back(NNN_ID[ix][count]);
                    count_site++;
                }
                if ((double(clock() - begin) / CLOCKS_PER_SEC) > 10)
                {
                    cout<<double(i_Fe)/N_Fe*100.0<<" percent done; site count: "<<count_site<<endl;
                    begin = clock();
                }
            }
            cout<<"Number of sites: "<<count_site<<endl;
        }
};

int main(int arg, char **name){
	cxxopts::Options options(name[0], "AKMC code");
	int rand_seed=0, nsteps=10000000;
	double temperature=300, dislocation_density=1.0e-6, box_height=27.9762;
    string octa_config="octa.dat", output_file="akmc.log";
	options.add_options()
		("T, temperature", "Temperature K (default: 300 K)", cxxopts::value<double>(temperature))
		("r, randseed", "Seed for random number (default: 0)", cxxopts::value<int>(rand_seed))
		("d, density", "Dislocation density in A^2 (default: 10^-6)", cxxopts::value<double>(dislocation_density))
        ("o, output", "Output file name (default: akmc.log)", cxxopts::value<string>(output_file))
        ("i, input", "Input file name (default: octa.dta)", cxxopts::value<string>(octa_config))
        ("z, height", "Box height (default: 27.9762)", cxxopts::value<double>(box_height))
		("h, help", "Print help")
	;
	auto results = options.parse(arg, name);
	srand(rand_seed);
	if (results.count("help"))
	{
		cout<<options.help({"", "Group"})<<endl;
		exit(EXIT_FAILURE);
	}
    AKMC akmc = AKMC(dislocation_density, octa_config, temperature, box_height);
}
