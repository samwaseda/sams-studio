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
#include "tree.cpp"
#include<algorithm>

using namespace std;

double zufall();

double square(double);

double quartic(double);

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
		void set_num_neighbors(int);
		bool set_type(string, bool);
		string get_type();
		float acceptance_ratio();
		double E(bool);
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

class Structure{
	private:
		int N_tot, N_rep, shell_max, N_s, N_v;
		string bravais;
		Atom* atom;
		fstream ausgabe;
	public:
		Structure(string, int, bool);
		Atom* get_structure();
		void create_lattice(Atom*, string);
		int get_number();
		void add_vacancy();
		void add_saddle();
		void set_coeff(string, int, bool);
		void set_solute(double);
		void reload_lattice(Atom*, string);
		bool check_vacancy(int, string, kdtree*, Shell*);
		void initialize_tree(kdtree*);
		float dist(Atom*, Atom*);
};

class Coeff{
	private:
		double ***JJ, **AA, **BB;
		unordered_map<string, int> ID_type = {{"Fe", 0}, {"Mn", 1}, {"Fes", 2}};
		unordered_map<string, int> defect_type = {{"b", 0}, {"v", 1}, {"s", 2}};
		int shell_max;
	public:
		Coeff(string, fstream&);
		void initialize_coeff();
		double A(string, string);
		double B(string, string);
		double J(string, string, int, string, string);
		int shell();
		int Elmax();
		int Defmax();
		int El(string);
		int Def(string);
};

class average_energy
{
	private:
		double EE, NN, E_sum, kappa;
	public:
		average_energy();
		void add(double, bool);
		double E();
		void reset();
};

class Energy{
	private:
		int acc, MC_count, N_tot;
		clock_t begin;
		bool debug_mode;
		double kB, lambda;
		Atom *atom;
		fstream ausgabe, config_ausgabe;
		default_random_engine generator;
		average_energy E_tot;
	public:
		Energy(Atom*, int, double, bool);
		double MC(double);
		void E_min();
		double output(string, bool);
		void reset();
};

#endif
