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

class CAtom : public Atom{
    private:
        double * E, * tau, tau_tot;
        Octa * octa;
        int id_octa;
        int zone_id:
    public:
        CAtom(){
            E = new double [4];
            tau = new double [4];
            tau_tot = 0;
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

class Energy{
    private:
        double a[7], b[5], dE0, c, sdash;
    public:
        Energy() : 
        a{9.0e-12, 6.0e-12, 8.0e-12, 2.2e-11, 3.5e-11, 3.0e-12, 1.7e-11},
        b{24.9e-7, 24.3e-7, 5.3e-7, 14.5e-7, 15.5e-7},
	    dE0(0.8153), c(0.829), sdash(1.08e4)
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
};

class AKMC{
    private:
        Octa * octa;
        CAtom * Catom;
        Eigen::Vector3d cell;
        double kBT, a_0;
        int N_octa, N_C;
        Eigen::Matrix3d tmat;
        Energy energy;
    public:
        AKMC(double dislocation_density, string input_file, double temperature_in, double box_height, double number_of_C_atoms) : a_0(2.855312531), N_octa(0)
        {
            cout<<"AKMC initialization"<<endl;
            tmat << 1.0/1.7320508, -1.0/1.4142135, 1.0/2.449489,
		            1.0/1.7320508, 0, -2.0/2.449489,
		            1.0/1.7320508, 1.0/1.4142135,  1.0/2.449489;
            kBT = 8.617e-5*temperature_in;
            if (kBT<=0)
                error_exit("ERROR: temperature must be a positive float");
            set_box(box_height, dislocation_density);
            set_number_of_C_atoms(number_of_C_atoms);
            read_octa(input_file);
            initialize_C();
        }

        void set_number_of_C_atoms(double nn)
        {
            if(nn<=0)
                error_exit("ERROR: Illegal number of C atoms");
            int N_Fe = int(2.0*cell.prod()/exp(3.0*log(a_0)));
            if(nn<1)
                nn = nn/(1.0-nn)*N_Fe;
            if(nn>0.1*N_Fe)
                cout<<"WARNING: Number of C atoms unrealistically high."<<endl;
            N_C = int(nn);
        }

        double read_octa(string input_file)
        {
            fstream eingabe;
            eingabe.open(input_file, ios::in);
            if(!eingabe)
                error_exit("ERROR: input file "+input_file+" does not exist or is empty");
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
                if(id_self==id_NN)
                    error_exit("ERROR: ID of neighbor is the ID of itself");
                if(id_NN>=N_octa)
                    error_exit("ERROR: neighbor ID larger than number of sites");
                octa[id_self].set_position(x);
                octa[id_self].set_neighbor(id_NN, S);
                octa[id_self].set_cell(&cell);
            }
            eingabe.close();
            cout<<"Octahedral interstitial site configuration import: complete"<<endl;
        }

        void set_box(double box_height, double dislocation_density)
        {
            if(box_height<=0)
                error_exit("ERROR: Problem with the box height");
            if(dislocation_density<=0)
                error_exit("ERROR: Dislocation density has to be a positive float");
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
	double temperature=300, dislocation_density=1.0e-6, box_height=27.9762, number_of_C_atoms=1;
    string octa_config="octa.dat", output_file="akmc.log";
	options.add_options()
		("T, temperature", "Temperature K (default: 300 K)", cxxopts::value<double>(temperature))
		("s, seed", "Seed for random number (default: 0)", cxxopts::value<int>(rand_seed))
		("d, density", "Dislocation density in A^2 (default: 10^-6)", cxxopts::value<double>(dislocation_density))
		("n, number", "Number of C atoms or fraction if smaller than 1 (default: 1)", cxxopts::value<double>(number_of_C_atoms))
        ("o, output", "Output file name (default: akmc.log)", cxxopts::value<string>(output_file))
        ("i, input", "Input file name (default: octa.dat)", cxxopts::value<string>(octa_config))
        ("z, height", "Box height (default: 27.9762)", cxxopts::value<double>(box_height))
        ("r, radius_core", "Radius of the core region in A (default: 10)", cxxopts::value<double>(r_core))
        ("R, radius_kmc", "Radius of the KMC region in A (default: 120)", cxxopts::value<double>(r_kmc))
		("h, help", "Print help")
	;
	auto results = options.parse(arg, name);
	srand(rand_seed);
	if (results.count("help"))
	{
		cout<<options.help({"", "Group"})<<endl;
		exit(EXIT_FAILURE);
	}
    AKMC akmc = AKMC(dislocation_density, octa_config, temperature, box_height, number_of_C_atoms);
}
