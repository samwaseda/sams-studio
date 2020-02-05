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

class Ridge{
	private:
		int N_line, N_dim, n_cv, n_set;
		double *yy, *increment, chi_old, val_error;
		bool true_error, weight, normalize;
		MatrixXd *H;
		VectorXd *hy, Lambda, w, normalizer;
		fstream lambda, coeff;
		OutputAndConsole output;
};

