#ifndef MC_H
#define MC_H

#include<random>

using namespace std;

double zufall();

double square(double);
double quartic(double);
double sextic(double);
double octic(double);
double decic(double);
double J_linear(double*, double*, double*);
double J_square(double*, double*, double*);
double J_cross_prod(double*, double*, double*);

class Atom{
    private:
        double mabs, mabs_old, theta, theta_old, phi, phi_old, *m_old, E_current[2], dE_current[2], dm, dphi, dtheta, mmax;
        vector<double> heisen_coeff[2], landau_coeff[2];
        vector<double (*)(double)> landau_func[2];
        vector<double (*)(double*, double*, double*)> heisen_func[2];
        vector<double*> m_n[2];
        int acc, count;
        bool E_uptodate[2], dE_uptodate[2], debug; // This does not work when neighbors change their m
        void update_flag(bool);
        void set_m(double, double, double);
    public:
        double *m;
        Atom();
        ~Atom();
        double get_acceptance_ratio();
        double E(int, bool);    // index, force_compute
        double dE(int, bool);   // index, force_compute
        void revoke();
        void set_landau_coeff(double, int, int);
        void set_heisenberg_coeff(double*, double, int, int);
        void clear_landau_coeff(int);
        void clear_heisenberg_coeff(int);
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
        int acc, MC_count, n_tot;
        clock_t steps_per_second;
        bool debug_mode;
        double kB, lambda, eta;
        Atom *atom;
        default_random_engine generator;
        average_energy E_tot;
        bool thermodynamic_integration();
    public:
        MC();
        ~MC();
        void create_atoms(int);
        void activate_debug();
        void run(double, int);
        bool accept();
        void set_lambda(double);
        void prepare_qmc(double, int);
        vector<double> get_magnetic_moments();
        void set_magnetic_moments(vector<double>);
        void set_landau_coeff(vector<double>, int, int);
        void set_heisenberg_coeff(vector<double>, vector<int>, vector<int>, int, int);
        void clear_landau_coeff(int);
        void clear_heisenberg_coeff(int);
        double get_acceptance_ratio();
        vector<double> get_acceptance_ratios();
        double get_energy(int);
        double get_mean_energy(int);
        double get_energy_variance(int);
        double get_steps_per_second();
        int get_number_of_atoms();
        void set_magnitude(vector<double>, vector<double>, vector<double>);
        void reset();
};

#endif
