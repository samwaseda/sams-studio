#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<algorithm>
#include<ctime>
#include<sstream>
#include<iomanip>
#include<Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "cxxopts.hpp"
#include "iofstream.cpp"
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

double square(double);

class Ridge{
	private:
		int N_line, n_dim, n_cv, n_set;
		double *increment, chi_old, val_error;
        vector<double> yy;
		bool true_error;
		MatrixXd *H;
		VectorXd *hy, lambda, coeff;
    public:
        void initialize_sets(vector<double>, vector<double>, bool);
        void ridge(double);
        void classic(int);
        void chi_training(bool);
        void least_square();
        void random(int);
        void raw(int);
};

