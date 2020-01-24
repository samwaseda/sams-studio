#ifndef MC_H
#define MC_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<string>
#include<iomanip>
#include<unordered_map>
#include<ctime>
#include<random>
#include<vector>
#include "cxxopts.hpp"
#include "tree.cpp"
#include<algorithm>

using namespace std;

double zufall(){
	return 1.0-2.0*(rand()/(double)RAND_MAX);
}

double square(double xxx){
	return xxx*xxx;
}

double quartic(double xxx){
	return square(xxx)*square(xxx);
}

namespace MC{
	class Atom{
		private:
			double mabs, mabs_old, theta, theta_old, phi, phi_old, *m_old, *J, A, B, **m_n, E_current, slope, m_0;
			int n_neigh, n_max, acc, count, N_spread;
			string type;
			bool E_uptodate; // This does not work when neighbors change their m
		public:
			double *m;
			float  *x;
			Atom();
			void set_num_neighbors(int num_neighbors=96);
			bool set_type(string type_in, bool restart=false);
			string get_type();
			float acceptance_ratio();
			double E(bool force_compute=false);
			double E_harmonic();
			double dE_harmonic();
			double dE();
			void set_m(double mabs_new, double theta_new, double phi_new);
			void revoke();
			void flip_z();
			int nx(int ix);
			void modify_AB(double A_in, double B_in);
			void set_AB(double A_in, double B_in);
			bool check_saddle();
			void set_neighbor(double* mm, double JJ);
			bool modify_neighbor(double* mm, double JJ);
			void propose_new_state();
	};

	class Coeff{
		private:
			double ***JJ, **AA, **BB;
			unordered_map<string, int> ID_type = {{"Fe", 0}, {"Mn", 1}, {"Fes", 2}};
			unordered_map<string, int> defect_type = {{"b", 0}, {"v", 1}, {"s", 2}};
			int shell_max;
		public:
			Coeff(string input_file, fstream &ausgabe);
			void initialize_coeff();
			double A(string elem, string def="b");
			double B(string elem, string def="b");
			double J(string elem_one, string elem_two, int shell_in, string def_one="b", string def_two="b");
			int shell();
			int Elmax();
			int Defmax();
			int El(string str);
			int Def(string str);
	};

	class Shell{
		private:
			float *dist_sb, *dist_bb, *dist_tmp;
			int shell_max;
			string bravais;
		public:
			Shell(int shell_max_in, string bravais);
			int get_shell(bool saddle, float distance)
			~Shell();
	};

	class Energy{
		private:
			int acc, MC_count, N_tot, N_rep, shell_max;
			clock_t begin = clock();
			bool debug_mode;
			double kB, lambda;
			Atom *atom;
			fstream ausgabe, config_ausgabe;
			default_random_engine generator;
			class average_energy
			{
				private:
					double EE, NN, E_sum, kappa;
				public:
					average_energy(double kappa_in=1.0-1.0e-3);
					void add(double E_in, bool total_energy=false);
					double E();
					void reset();
			} E_tot;
		public:
			Energy(string input_file, string bravais, int N_in, float N_Mn_in, int N_s, int N_v, double lambda_in, int num_neighbors, bool debug_mode_in);

			void create_lattice(Atom *atom, string bravais, int N_s=0, int N_Mn=0);
			void reload_lattice(Atom *atom, string restart_file="config.dat");
			bool check_vacancy(int index, string bravais, kdtree *tree, Shell *shell);
			void initialize_tree(kdtree *tree);
			float dist(Atom *xx, Atom *yy);
			double MC(double T_in);
			void E_min();
			double output(string custom_text="", bool config=false);
			void reset();
	};

	/**
	int main(int arg, char ** name){
		cxxopts::Options options(name[0], "Monte carlo simulation code");
		int N_s=0, N_v=0, rand_seed=0, nsteps=10000000, N=10, num_neighbors=58;
		double T=300, lambda=1;
		float N_Mn=0;
		bool debug_mode = false, antiferro=false;
		string input_file="jij_parameters.dat", bravais_lattice="bcc";
		options.add_options()
			("i, input", "Input file for Jij parameters (default: jij_parameters.dat)", cxxopts::value<string>(input_file))
			("b, bravais", "Bravais lattice, bcc or fcc (default: bcc)", cxxopts::value<string>(bravais_lattice))
			("T, temperature", "Temperature K (default: 300 K)", cxxopts::value<double>(T))
			("L, lambda", "Lambda value for thermodynamic integration (default: 1)", cxxopts::value<double>(lambda))
			("n, number-of-atoms", "if n<1000: number of cells in each direction, otherwise number of atoms (default: 10)", cxxopts::value<int>(N))
			("N, number-of-steps", "Number of Monte Carlo steps (default: 10,000,000)", cxxopts::value<int>(nsteps))
			("s, saddle", "Number of saddle point atoms (default: 0)", cxxopts::value<int>(N_s))
			("v, vacancy", "Number of vacancies (default: 0)", cxxopts::value<int>(N_v))
			("r, randseed", "Seed for random number (default: 0)", cxxopts::value<int>(rand_seed))
			("m, manganese", "Number or fraction of Mn atoms (default: 0)", cxxopts::value<float>(N_Mn))
			("d, debug", "Debug mode (not so helpful for now...)")
			("a, anneal", "Annealing before the beginning of the simulation (default: false)")
			("h, help", "Print help")
		;
		auto results = options.parse(arg, name);
		srand(rand_seed);
		if (results.count("help"))
		{
			cout<<options.help({"", "Group"})<<endl;
			exit(0);
		}
		if (results.count("debug"))
			debug_mode = true;
		if (bravais_lattice=="fcc")
			num_neighbors=96;

		Energy MC(input_file, bravais_lattice, N, N_Mn, N_s, N_v, lambda, num_neighbors, debug_mode);
		if(results.count("anneal"))
		{
			int N_anneal = 10000;
			for(int i=0; i<N_anneal; i++)
				MC.MC(T*(double)i/(double)N_anneal);
		}
		MC.reset();
		for(int i=0; i<=10000000; i++)
		{
			MC.MC(T);
			if(i%1000==0)
			{
				cout<<T<<" "<<i;
				MC.output(to_string(int(T)), i%100000==0);
			}
		}
	}
	**/
}

#endif
