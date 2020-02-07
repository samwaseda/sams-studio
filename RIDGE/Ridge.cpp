#include "Ridge.h"

double square(double x){
    return x*x;
}

Ridge::Ridge() : n_dim(0), n_cv(10), chi_old(0), val_error(0), lambda_tol(100), debug(false){}

vector<double> Ridge::get_determinant()
{
    vector<double> determinant;
    for(int i=0; i<n_set(); i++)
        determinant.push_back(H[i].determinant());
    return determinant;
}

vector<double> Ridge::get_sum_H()
{
    vector<double> sum_H;
    for(int i=0; i<n_set(); i++)
        sum_H.push_back(H[i].sum());
    return sum_H;
}

void Ridge::activate_debug(){
    debug = true;
}

void Ridge::set_number_of_cv_set(int n_cv_in)
{
    if(n_cv_in<=0)
        throw invalid_argument("Number of cross validation set has to be larger than one");
    if(n_dim!=0)
        throw invalid_argument("Number of cross validation set cannot be changed after initialization");
    n_cv = n_cv_in;
}

void Ridge::chi_training(bool regulation=true){
    if(n_dim<=0)
        throw invalid_argument("data set not initialized");
    val_error = 0;
    for(int tr_set=0; tr_set<n_cv; tr_set++)
    {
        MatrixXd HL_tmp = MatrixXd::Zero(n_dim, n_dim);
        VectorXd hy_tmp = VectorXd::Zero(n_dim);
        for(int j=0; j<n_cv; j++)
            if(j!=tr_set)
            {
                HL_tmp += H[j];
                hy_tmp += hy[j];
            }
        if (regulation)
        {
            for(int j=0; j<n_dim; j++)
                HL_tmp(j,j) += exp(lambda(j));
        }
        coeff = HL_tmp.inverse()*hy_tmp;
        val_error += yy[tr_set]+ (double) (-2.0*coeff.transpose()*hy[tr_set]+(double)(coeff.transpose()*H[tr_set]*coeff));
    }
    val_error = sqrt(val_error/n_cv);
    if(debug && coeff.norm()==0)
        throw invalid_argument("All coefficients are zero");
}

double Ridge::get_validation_error(){
    return val_error;
}

double Ridge::get_true_error(){
    if(true_error)
        return 0;
    return sqrt(yy[n_cv]-2*coeff.transpose()*hy[n_cv]+coeff.transpose()*H[n_cv]*coeff);
}

void Ridge::reset_increment(){
    increment.resize(n_dim, 1.0);
}

int Ridge::n_set(){
    if(true_error)
        return n_cv+1;
    return n_cv;
}

void Ridge::initialize_sets(vector<double> x_in, vector<double> y_in, int zeroth, bool true_error_in)
{
    zeroth = (zeroth>0);
    true_error = true_error_in;
    int n_data = int(y_in.size());
    n_dim = int(int(x_in.size())/int(n_data));
    if(int(x_in.size())!=n_data*n_dim)
        throw invalid_argument("Dimension of x is not wrong");
    if(zeroth>0)
        n_dim++;
    H = new MatrixXd[n_set()];
    hy = new VectorXd[n_set()];
    yy.resize(n_set(), 0);
    double s_in[n_dim];
    lambda = VectorXd::Zero(n_dim);
    coeff = VectorXd::Zero(n_dim);
    reset_increment();
    for(int i=0; i<n_set(); i++)
    {
        H[i] = MatrixXd::Zero(n_dim, n_dim);
        hy[i] = VectorXd::Zero(n_dim);
    }
    if (zeroth>0)
        s_in[0]=1.0;
    for(int line=0; line<n_data; line++)
    {
        int ID_tmp = line%(n_set());
        for(int i=0; i<n_dim; i++)
            s_in[i+int(zeroth)] = x_in[line*(n_dim-int(zeroth))+i];
        for(int i=0; i<n_dim; i++)
        {
            for(int j=0; j<n_dim; j++)
                H[ID_tmp](i,j) += s_in[i]*s_in[j];
            hy[ID_tmp](i) += s_in[i]*y_in.at(line);
        }
        yy.at(ID_tmp) += square(y_in.at(line));
    }
    for(int i=0; debug && i<n_set(); i++)
        if(hy[i].norm()==0)
            throw invalid_argument("There was a problem with initialization");
}

void Ridge::ridge(double lambda_in)
{
    lambda = VectorXd::Constant(n_dim, lambda_in);
    chi_training();
}

void Ridge::classic(int max_step){
    double lambda_old = lambda.sum()/n_dim;
    lambda = VectorXd::Constant(n_dim, lambda_old);
    chi_training();
    double dlambda = 1.0, chi_old=val_error;
    for(int i=0; i<max_step && abs(dlambda)>1.0e-8; i++)
    {
        ridge(lambda_old+dlambda);
        if(val_error>chi_old)
        {
            dlambda *= -0.9;
            lambda = VectorXd::Constant(n_dim, lambda_old);
        }
        else
        {
            chi_old = val_error;
            lambda_old += dlambda;
        }
    }
}

double Ridge::chi(bool ls){
    MatrixXd HL_tmp = MatrixXd::Zero(n_dim, n_dim);
    VectorXd hy_tmp = VectorXd::Zero(n_dim);
    for(int i=0; i<n_set(); i++)
    {
        HL_tmp += H[i];
        hy_tmp += hy[i];
    }
    for(int i=0; i<n_dim; i++)
        HL_tmp(i,i) += exp(lambda(i));
    coeff = HL_tmp.inverse()*hy_tmp;
    if(true_error)
        return sqrt(yy[n_cv]-2*coeff.transpose()*hy[n_cv]+coeff.transpose()*H[n_cv]*coeff);
    else
        return 0;
}

void Ridge::least_square(){
    chi_training(false);
}

void Ridge::random(int max_step){
    if(chi_old == 0)
        chi_training();
    chi_old = val_error;
    int max_number = n_dim;
    for(int i=0; i<max_step; i++)
    {
        for(int j=0; j<max_number; j++)
        {
            int rand_ID = rand()%n_dim;
            double lambda_old = lambda(rand_ID);
            lambda(rand_ID) = lambda_tol*(1.0-2.0*rand()/(double)RAND_MAX);
            chi_training();
            if(val_error<chi_old)
            {
                chi_old = val_error;
                max_number++;
            }
            else
            {
                lambda(rand_ID) = lambda_old;
                max_number -= ((rand()%100)==0);
            }
        }
    }
}

Ridge::~Ridge(){
    delete [] H;
    delete [] hy;
}

void Ridge::raw(int max_step){
    for(int i=0; i<max_step; i++)
    {
        double chi_init = val_error;
        for(int j=0; j<n_dim; j++)
        {
            int ID = rand()%n_dim;
            bool flag_stop = true;
            for(int l=0; l<i+1 && flag_stop && l<10; l++)
            {
                double lambda_old = lambda(ID);
                if(chi_old != 0)
                    lambda(ID) += increment[ID];
                chi_training();
                if(chi_old == 0)
                {
                    chi_old = val_error;
                    chi_init = chi_old;
                    l--;
                }
                else if(val_error>=chi_old || abs(lambda(ID))>lambda_tol)
                {
                    lambda(ID) = lambda_old;
                    increment[ID] *= -0.9;
                    flag_stop = false;
                }
                else
                {
                    chi_old = val_error;
                    increment[ID] *= 1.001;
                }
            }
        }
        if(i>10 && abs(chi_init-val_error)<0.00001*val_error)
            break;
    }
}

void Ridge::set_lambda(vector<double> lambda_in){
    if (int(lambda_in.size())!=n_dim)
        throw invalid_argument("Length of lambda incorrect");
    for(int i=0; i<n_dim; i++)
        lambda(i) = lambda_in.at(i);
    chi_training();
}

VectorXd Ridge::dchi(){
    VectorXd dchi_dlambda = VectorXd::Zero(n_dim);
    val_error = 0;
    for(int tr_set=0; tr_set<n_cv; tr_set++)
    {
        MatrixXd HL_tmp = MatrixXd::Zero(n_dim, n_dim);
        VectorXd hy_tmp = VectorXd::Zero(n_dim);
        VectorXd dchi_1 = VectorXd::Zero(n_dim);
        VectorXd dchi_2 = VectorXd::Zero(n_dim);
        for(int j=0; j<n_cv; j++)
            if(j!=tr_set)
            {
                HL_tmp += H[j];
                hy_tmp += hy[j];
            }
        for(int j=0; j<n_dim; j++)
            HL_tmp(j,j) += exp(lambda(j));
        //w = HL_tmp.inverse()*hy_tmp;
        //val_error += yy[tr_set]+ (double) (-2.0*w.transpose()*hy[tr_set]+(double)(w.transpose()*H[tr_set]*w));
        MatrixXd HL_inv = HL_tmp.inverse();
        dchi_1 = HL_inv*hy_tmp;
        dchi_2 = HL_inv*(-hy[tr_set]+H[tr_set]*HL_inv*hy_tmp);
        for(int j=0; j<n_dim; j++)
            dchi_dlambda(j) += -2*dchi_1(j)*dchi_2(j);
    }
    //val_error = sqrt(val_error/N_line/n_cv*(1.0+n_cv));
    for(int i=0; i<n_dim; i++)
        dchi_dlambda(i) *= exp(lambda(i)-lambda.mean());
    return dchi_dlambda;
}

MatrixXd Ridge::ddchi(){
    MatrixXd ddchi_ddlambda = MatrixXd::Zero(n_dim, n_dim);
    val_error = 0;
    for(int tr_set=0; tr_set<n_cv; tr_set++)
    {
        MatrixXd HL_tmp = MatrixXd::Zero(n_dim, n_dim);
        VectorXd hy_tmp = VectorXd::Zero(n_dim);
        MatrixXd ddchi_1 = MatrixXd::Zero(n_dim, n_dim);
        MatrixXd ddchi_21 = MatrixXd::Zero(n_dim, n_dim);
        MatrixXd ddchi_22 = MatrixXd::Zero(n_dim, n_dim);
        for(int i=0; i<n_cv; i++)
            if(i!=tr_set)
            {
                HL_tmp += H[i];
                hy_tmp += hy[i];
            }
        for(int i=0; i<n_dim; i++)
            HL_tmp(i,i) += exp(lambda(i));
        //w = HL_tmp.inverse()*hy_tmp;
        //val_error += yy[tr_set]+ (double) (-2.0*w.transpose()*hy[tr_set]+(double)(w.transpose()*H[tr_set]*w));
        MatrixXd HL_inv = HL_tmp.inverse();
        ddchi_1 = 2*HL_inv*(-hy[tr_set]+H[tr_set]*HL_inv*hy_tmp)*hy_tmp.transpose()*HL_inv;
        //ddchi_1 += HL_inv*hy_tmp*(-hy[tr_set].transpose()+hy_tmp.transpose()*HL_inv*H[tr_set])*HL_inv;
        ddchi_21 = HL_inv*H[tr_set]*HL_inv;
        ddchi_22 = HL_inv*hy_tmp*hy_tmp.transpose()*HL_inv;
        for(int j=0; j<n_dim; j++)
            for(int k=0; k<n_dim; k++)
            {
                ddchi_1(j,k) *= HL_inv(j,k);
                ddchi_21(j,k) *= ddchi_22(j,k);
            }
        ddchi_ddlambda += 2*(ddchi_1+ddchi_21);
    }
    //val_error = sqrt(val_error/N_line/n_cv*(1.0+n_cv));
    for(int i=0; i<n_dim; i++)
        for(int j=0; j<n_dim; j++)
            ddchi_ddlambda(i,j) *= exp(lambda(i)+lambda(j)-lambda.mean());
    return ddchi_ddlambda;
}

void Ridge::gradient_descent(int max_cycle, double prefactor, double damper=1.0)
{
    for(int i_cycle=0; i_cycle<max_cycle; i_cycle++)
    {
        lambda -= prefactor*dchi();
        chi_training();
        prefactor *= damper;
        if(lambda.maxCoeff()>lambda_tol || lambda.minCoeff()<-lambda_tol)
            i_cycle = max_cycle;
    }
}

void Ridge::conjugate_gradient(int max_cycle){
    VectorXd r_cg = VectorXd::Zero(n_dim);
    VectorXd p_cg = VectorXd::Zero(n_dim);
    double alpha = 0, beta = 0;
    VectorXd dchi_dlambda = VectorXd::Zero(n_dim);
    MatrixXd ddchi_ddlambda = MatrixXd::Zero(n_dim, n_dim);
    for(int i=0; i<max_cycle+1; i++)
    {
        dchi_dlambda = dchi();
        ddchi_ddlambda = ddchi();
        if(i==0)
        {
            //r_cg = -dchi_dlambda-ddchi_ddlambda*lambda;
            r_cg = -ddchi_ddlambda*lambda;
            p_cg = r_cg;
        }
        alpha = (double) (r_cg.transpose()*r_cg)/(p_cg.transpose()*ddchi_ddlambda*p_cg);
        if(abs(alpha)<1.0e-20 || !isnormal(alpha))
            i=max_cycle;
        else
        {
            lambda += alpha*p_cg;
            beta = (double) (r_cg.transpose()*r_cg);
            r_cg += -alpha*ddchi_ddlambda*p_cg;
            beta = (double) (r_cg.transpose()*r_cg)/beta;
            p_cg = r_cg+beta*p_cg;
            Eigen::EigenSolver<MatrixXd> es(ddchi_ddlambda,false);
        }
    }
}


vector<double> Ridge::get_coeff()
{
        vector<double> vec(coeff.begin(), coeff.end());
        return vec;
}

vector<double> Ridge::get_lambda()
{
        vector<double> vec(lambda.begin(), lambda.end());
        return vec;
}

vector<double> Ridge::get_derivative(){
    VectorXd d = dchi();
    vector<double> v(d.begin(), d.end());
    return v;
}

vector<double> Ridge::get_hessian(){
    VectorXd ddchi_ddlambda = ddchi();
    vector<double> v(ddchi_ddlambda.data(), ddchi_ddlambda.data() + ddchi_ddlambda.rows() * ddchi_ddlambda.cols());
    return v;
}
