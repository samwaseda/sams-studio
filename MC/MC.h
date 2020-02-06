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
        double mabs, mabs_old, theta, theta_old, phi, phi_old, *m_old, E_current, dE_current, dm, dphi, dtheta;
        vector<double> J[2], A, B;
        vector<double*> m_n[2];
        int acc, count;
        bool E_uptodate[2], dE_uptodate[2], debug; // This does not work when neighbors change their m
        void update_flag(bool);
        void set_m(double, double, double);
    public:
        double *m;
        Atom();
        ~Atom();
        float get_acceptance_ratio();
        double E(bool, int);    // force_compute, index
        double dE(bool, int);   // force_compute, index
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
        double EE[2], E_sum[2], EE_sq[2];
        int NN[2];
    public:
        average_energy();
        void add(double, bool, int);
        double E_mean(int);
        double E_var(int);
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
        void create_atoms(vector<double>, vector<double>, vector<int>, vector<int>, vector<double>);
        void append_parameters(vector<double>, vector<double>, vector<int>, vector<int>, vector<double>);
        void activate_debug();
        void run(double, int);
        void set_lambda(double);
        vector<double> get_magnetic_moments();
        double get_acceptance_ratio();
        double get_energy(int);
        double get_mean_energy(int);
        double get_energy_variance(int);
        double get_steps_per_second();
        void set_magnitude(double, double, double);
        void reset();
};

#endif
