#ifndef MC_H
#define MC_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<string>
#include<iomanip>
#include <stdexcept>
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
        double mabs, mabs_old, theta, theta_old, phi, phi_old, *m_old, **J, A[2], B[2], **m_n, E_current, dE_current, dm, dphi, dtheta;
        int n_neigh[2], n_max, acc, count;
        bool E_uptodate[2], dE_uptodate[2], debug; // This does not work when neighbors change their m
        void update_flag(bool);
    public:
        double *m;
        Atom();
        ~Atom();
        void set_num_neighbors(int);
        float acceptance_ratio();
        double E(bool, int);
        double dE(bool, int);
        void set_m(double, double, double);
        void revoke();
        void set_AB(double, double, int);
        void set_neighbor(double*, double, int);
        void activate_debug();
        void propose_new_state();
        void set_magnitude(double, double, double);
};

class average_energy
{
    private:
        double EE, E_sum, EE_sq;
        int NN;
    public:
        average_energy();
        void add(double, bool);
        double E_mean();
        double E_var();
        void reset();
};

class MC{
    private:
        int acc, MC_count, N_tot, thermodynamic_integration_flag;
        clock_t steps_per_second;
        bool debug_mode;
        double kB, lambda;
        Atom *atom;
        default_random_engine generator;
        average_energy E_tot;
        bool thermodynamic_integration();
    public:
        MC();
        ~MC();
        void create_atoms(int, vector<double>, vector<double>, vector<int>, vector<int>, vector<double>);
        void append_parameters(vector<double>, vector<double>, vector<int>, vector<int>, vector<double>);
        void activate_debug();
        double run(double, int);
        void set_lambda(double);
        vector<double> get_magnetic_moments();
        double get_acceptance_ratio();
        double get_energy();
        double get_mean_energy();
        double get_energy_variance();
        double get_steps_per_second();
        void set_magnitude(double, double, double);
        void reset();
};

#endif
