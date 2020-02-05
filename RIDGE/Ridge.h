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
		bool true_error;
		MatrixXd *H;
		VectorXd *hy, lambda, coeff;
    public:
        Ridge();
        ~Ridge();
        void set_number_of_cv_set(int);
        double get_determinant();
        void initialize_sets(vector<double>, vector<double>, bool);
        void ridge(double);
        void classic(int);
        void chi_training(bool);
        double chi(bool);
        void least_square();
        void random(int);
        void raw(int);
};

#endif
