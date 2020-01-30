#ifndef MC_H
#define MC_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<string>
#include<iomanip>
#include<ctime>
#include<random>
#include<vector>
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
		int get_num_neighbors();
		bool modify_neighbor(double*, double);
		void propose_new_state();
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
