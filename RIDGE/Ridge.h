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
		int n_dim, n_cv, n_set;
		double *increment, chi_old, val_error;
        vector<double> yy;
		bool true_error, debug;
		MatrixXd *H;
		VectorXd *hy, lambda, coeff;
        void chi_training(bool);
    public:
        Ridge();
        ~Ridge();
        void activate_debug();
        void set_number_of_cv_set(int);
        double get_determinant();
        double get_validation_error();
        double get_true_error();
        void initialize_sets(vector<double>, vector<double>, bool, bool);
        void ridge(double);
        void classic(int);
        double chi(bool);
        void least_square();
        void random(int);
        void raw(int);
        vector<double> get_coeff();
        vector<double> get_lambda();
};

#endif
