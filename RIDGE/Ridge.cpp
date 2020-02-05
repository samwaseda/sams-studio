#include "Ridge.h"

double square(double x){
    return x*x;
}

Ridge::Ridge() : n_dim(0), n_cv(10), chi_old(0), val_error(0), debug(false){}

double Ridge::get_determinant()
{
    double dd_min=0, dd_max, dd=0, dd_sum, dd_var;
    for(int i=0; i<n_set(); i++)
    {
        dd = H[i].determinant();
        if(i==0 || dd_min>dd)
            dd_min = dd;
        if(i==0 || dd_max<dd)
            dd_max = dd;
        dd_sum += dd;
        dd_var += dd*dd;
    }
    return dd_min;
    //output<<"Determinants: min = "<<dd_min<<", max = "<<dd_max<<", sum = "<<dd_sum<<" += "<<sqrt(dd_var-dd_sum*dd_sum/n_set)<<"\n";
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
                HL_tmp(j,j) += exp(-lambda(j));
        }
        coeff = HL_tmp.inverse()*hy_tmp;
        val_error += yy[tr_set]+ (double) (-2.0*coeff.transpose()*hy[tr_set]+(double)(coeff.transpose()*H[tr_set]*coeff));
    }
    val_error = sqrt(val_error/n_cv);
    if(coeff.norm()==0)
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

void Ridge::initialize_sets(vector<double> x_in, vector<double> y_in, bool zeroth, bool true_error_in)
{
    true_error = true_error_in;
    int n_data = int(y_in.size());
    n_dim = int(int(x_in.size())/int(n_data));
    if(zeroth)
        n_dim++;
    H = new MatrixXd[n_set()];
    hy = new VectorXd[n_set()];
    yy.resize(n_set(), 0);
    double s_in[n_dim];
    lambda = VectorXd::Zero(n_dim);
    reset_increment();
    for(int i=0; i<n_set(); i++)
    {
        H[i] = MatrixXd::Zero(n_dim, n_dim);
        hy[i] = VectorXd::Zero(n_dim);
    }
    coeff = VectorXd::Zero(n_dim);
    for(int line=0; line<n_data; line++)
    {
        int ID_tmp = line%(n_set());
        if (zeroth)
            s_in[0]=1.0;
        for(int i=int(zeroth); i<n_dim; i++)
            s_in[i] = x_in[line*n_dim+i];
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
        HL_tmp(i,i) += exp(-lambda(i));
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
            lambda(rand_ID) = 100.0*(1.0-2.0*rand()/(double)RAND_MAX);
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
                else if(val_error>=chi_old || abs(lambda(ID))>100)
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
