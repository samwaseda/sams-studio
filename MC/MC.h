#ifndef MC_H
#define MC_H

#include<cstdlib>
#include<ctime>
#include<valarray>
#include <chrono>

using namespace std;

double zufall();
double power(double, int);

struct Product;
struct Magnitude;;

class Atom{
    private:
        double mabs, mabs_old, E_current[2], dE_current[2], dm, dphi, mmax;
        valarray<double> gradient;
        vector<double> heisen_coeff[2], landau_coeff[2];
        vector<Magnitude*> landau_func[2];
        vector<Product*> heisen_func[2];
        vector<Atom*> neigh[2];
        int acc, count;
        bool E_uptodate[2], dE_uptodate[2], debug, flip; // This does not work when neighbors change their m
        void update_flag(bool ff=false);
        friend Product;
    public:
        valarray<double> m, m_old;
        valarray<double> get_gradient(double); // lambda
        Atom();
        ~Atom();
        void set_m(valarray<double>&, bool diff=false);
        double get_gradient_residual();
        double get_acceptance_ratio();
        double get_magnitude(int exponent=1, bool old=false);
        double E(int index=0, bool force_compute=false);
        double dE(int index=0, bool force_compute=false);
        double run_gradient_descent(double, double); // h lambda diff
        void revoke();
        void set_landau_coeff(double, int, int);
        void set_heisenberg_coeff(Atom&, double, int deg=1, int index=0);
        void clear_landau_coeff(int);
        void clear_heisenberg_coeff(int);
        void activate_debug();
        void propose_new_state();
        void set_magnitude(double, double, bool flip_in=true);
        void check_consistency();
};

class average_energy
{
    private:
        double EE[2], E_sum[2], EE_sq[2];
        int NN[2];
    public:
        average_energy();
        void add(double, bool total_energy=false, int index=0);
        double E_mean(int);
        double E_var(int);
        void reset();
};

class MC{
    private:
        long long int acc, MC_count;
        int n_tot;
        double steps_per_second;
        bool debug_mode;
        double kB, lambda, eta, E_min;
        Atom *atom;
        average_energy E_tot;
        bool thermodynamic_integration();
        bool accept(int, double, double);
        bool bose_einstein();
        vector<int> selectable_ID;
    public:
        MC();
        ~MC();
        void create_atoms(int);
        void activate_debug();
        void run(double, int number_of_iterations=1);
        void set_lambda(double);
        void set_eta(double);
        double get_eta();
        vector<double> get_magnetic_moments();
        vector<double> get_magnetic_gradients();
        void set_magnetic_moments(vector<double>);
        void set_landau_coeff(vector<double>, int, int index=0);
        void set_heisenberg_coeff(vector<double>, vector<int>, vector<int>, int, int index=0);
        void clear_landau_coeff(int index=0);
        void clear_heisenberg_coeff(int index=0);
        double get_acceptance_ratio();
        vector<double> get_acceptance_ratios();
        double get_energy(int);
        double get_mean_energy(int index=0);
        double get_energy_variance(int index=0);
        double get_ground_state_energy();
        double get_steps_per_second();
        int get_number_of_atoms();
        void run_debug();
        void set_magnitude(vector<double>, vector<double>, vector<int>);
        double run_gradient_descent(int, double step_size=1, double decrement=0.001, double diff=1.0e-8);
        void select_ID(vector<int>);
        void reset();
        void set_metadynamics(double, double, double, int);
};

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

struct Product{
    virtual double value(Atom&, Atom&);
    virtual double diff(Atom&, Atom&);
    virtual valarray<double> gradient(Atom&, Atom&);
};

struct J_lin_lin : Product {
    double value(Atom&, Atom&);
    double diff(Atom&, Atom&);
    valarray<double> gradient(Atom&, Atom&);
} j_lin_lin;

struct J_cub_lin : Product {
    double value(Atom&, Atom&);
    double diff(Atom&, Atom&);
    valarray<double> gradient(Atom&, Atom&);
} j_cub_lin;

struct J_qui_lin : Product {
    double value(Atom&, Atom&);
    double diff(Atom&, Atom&);
    valarray<double> gradient(Atom&, Atom&);
} j_qui_lin;

struct J_cross_forward : Product {
    double value(Atom&, Atom&);
    double diff(Atom&, Atom&);
    valarray<double> gradient(Atom&, Atom&);
} j_cross_forward;

struct J_cross_backward : Product {
    double value(Atom&, Atom&);
    double diff(Atom&, Atom&);
    valarray<double> gradient(Atom&, Atom&);
} j_cross_backward;

#endif
