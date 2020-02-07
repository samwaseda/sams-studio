#ifndef RIDGE_H
#define RIDGE_H

#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<algorithm>
#include<ctime>
#include<iomanip>
#include<Eigen/Dense>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

double square(double);

class Ridge{
    private:
        int n_dim, n_cv;
        double chi_old, val_error, lambda_tol;
        vector<double> yy, increment;
        bool true_error, debug;
        MatrixXd *H;
        VectorXd *hy, lambda, coeff;
        void chi_training(bool);
        int n_set();
        VectorXd dchi();
        MatrixXd ddchi();
        double chi(bool);
    public:
        Ridge();
        ~Ridge();
        void activate_debug();
        void set_number_of_cv_set(int);
        double get_validation_error();
        double get_true_error();
        void initialize_sets(vector<double>, vector<double>, int, int);
        void ridge(double);
        void classic(int);
        void least_square();
        void random(int);
        void raw(int);
        void conjugate_gradient(int);
        void gradient_descent(int, double, double);
        void reset_increment();
        void set_lambda(vector<double>);
        vector<double> get_coeff();
        vector<double> get_lambda();
        vector<double> get_determinant(int);
        vector<double> get_sum_H();
        vector<double> get_derivative();
        vector<double> get_hessian();
};

#endif
