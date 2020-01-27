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

double zufall();

double square(double);

double quartic(double);

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
			bool set_type(string, bool restart=false);
			string get_type();
			float acceptance_ratio();
			double E(bool force_compute=false);
			double E_harmonic();
			double dE_harmonic();
			double dE();
			void set_m(double, double, double);
			void revoke();
			void flip_z();
			int nx(int ix);
			void modify_AB(double, double);
			void set_AB(double, double);
			bool check_saddle();
			void set_neighbor(double*, double);
			bool modify_neighbor(double*, double);
			void propose_new_state();
	};

	class Coeff{
		private:
			double ***JJ, **AA, **BB;
			unordered_map<string, int> ID_type = {{"Fe", 0}, {"Mn", 1}, {"Fes", 2}};
			unordered_map<string, int> defect_type = {{"b", 0}, {"v", 1}, {"s", 2}};
			int shell_max;
		public:
			Coeff(string, fstream);
			void initialize_coeff();
			double A(string, string def="b");
			double B(string, string def="b");
			double J(string, string, int, string def_one="b", string def_two="b");
			int shell();
			int Elmax();
			int Defmax();
			int El(string);
			int Def(string);
	};

	class Shell{
		private:
			float *dist_sb, *dist_bb, *dist_tmp;
			int shell_max;
			string bravais;
		public:
			Shell(int, string);
			int get_shell(bool, float);
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
					void add(double, bool total_energy=false);
					double E();
					void reset();
			} E_tot;
		public:
			Energy(string, string, int, float, int, int, double, int, bool);
			void create_lattice(Atom*, string, int N_s=0, int N_Mn=0);
			void reload_lattice(Atom*, string restart_file="config.dat");
			bool check_vacancy(int, string, kdtree*, Shell*);
			void initialize_tree(kdtree*);
			float dist(Atom*, Atom*);
			double MC(double );
			void E_min();
			double output(string custom_text="", bool config=false);
			void reset();
	};
}

#endif
