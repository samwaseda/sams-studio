#ifndef MC_H
#define MC_H

#include<random>
#include<valarray>

using namespace std;

double zufall();

struct Magnitude{
    virtual double value(double);
    virtual valarray<double> gradient(valarray<double>&);
};

struct Square : Magnitude {
    double value(double);
    valarray<double> gradient(valarray<double>&);
} square;

struct Quartic : Magnitude {
    double value(double);
    valarray<double> gradient(valarray<double>&);
} quartic;

struct Sextic : Magnitude {
    double value(double);
    valarray<double> gradient(valarray<double>&);
} sextic;

struct Octic : Magnitude {
    double value(double);
    valarray<double> gradient(valarray<double>&);
} octic;

struct Decic : Magnitude {
    double value(double);
    valarray<double> gradient(valarray<double>&);
} decic;

struct Bilinear{
    virtual double value(valarray<double>&, valarray<double>&);
    virtual double diff(valarray<double>&, valarray<double>&, valarray<double>&);
    virtual valarray<double> gradient(valarray<double>&, valarray<double>&);
};

struct J_linear : Bilinear {
    double value(valarray<double>&, valarray<double>&);
    double diff(valarray<double>&, valarray<double>&, valarray<double>&);
    valarray<double> gradient(valarray<double>&, valarray<double>&);
} j_linear;

struct J_square : Bilinear {
    double value(valarray<double>&, valarray<double>&);
    double diff(valarray<double>&, valarray<double>&, valarray<double>&);
    valarray<double> gradient(valarray<double>&, valarray<double>&);
} j_square;


// double J_cross_prod(double*, double*, double*);

class Atom{
    private:
        double mabs, mabs_old, theta, theta_old, phi, phi_old, E_current[2], dE_current[2], dm, dphi, dtheta, mmax;
        valarray<double> m_old, gradient;
        vector<double> heisen_coeff[2], landau_coeff[2];
        vector<Magnitude*> landau_func[2];
        vector<Bilinear*> heisen_func[2];
        vector<valarray<double> * > m_n[2];
        int acc, count;
        bool E_uptodate[2], dE_uptodate[2], debug; // This does not work when neighbors change their m
        void update_flag(bool);
        void set_m(double, double, double);
        double get_gradient_residual();
    public:
        valarray<double> m;
        valarray<double> get_gradient(double); // lambda
        Atom();
        ~Atom();
        double get_acceptance_ratio();
        double E(int, bool);    // index, force_compute
        double dE(int, bool);   // index, force_compute
        double run_gradient_descent(double, double, double); // h lambda diff
        void revoke();
        void set_landau_coeff(double, int, int);
        void set_heisenberg_coeff(valarray<double>*, double, int, int);
        void clear_landau_coeff(int);
        void clear_heisenberg_coeff(int);
        void activate_debug();
        void propose_new_state();
        void set_magnitude(double, double, double);
        void update_polar_coordinates();
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
        double kB, lambda, eta, E_min;
        Atom *atom;
        default_random_engine generator;
        average_energy E_tot;
        bool thermodynamic_integration();
        bool accept(int, double, double);
        bool preparing_qmc();
    public:
        MC();
        ~MC();
        void create_atoms(int);
        void prepare_qmc(double, int);
        void activate_debug();
        void run(double, int);
        void set_lambda(double);
        void set_eta(double);
        double get_eta();
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
        double run_gradient_descent(int, double, double);
        void reset();
};

#endif
