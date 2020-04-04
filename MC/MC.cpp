#include "MC.h"

using namespace std;

double zufall(){
    return 1.0-2.0*(rand()/(double)RAND_MAX);
}

double sq_mag(valarray<double> &m){
    return (m*m).sum();
}

double Magnitude::value(double xxx){return 0;}
double Square::value(double xxx){ return xxx*xxx; }
double Quartic::value(double xxx){ return square.value(xxx)*square.value(xxx); }
double Sextic::value(double xxx){ return quartic.value(xxx)*square.value(xxx); }
double Octic::value(double xxx){ return quartic.value(xxx)*quartic.value(xxx); }
double Decic::value(double xxx){ return sextic.value(xxx)*quartic.value(xxx); }

valarray<double> Magnitude::gradient(valarray<double> &m){
    return 0*m;
}

valarray<double> Square::gradient(valarray<double> &m){
    return 2*m;
}

valarray<double> Quartic::gradient(valarray<double> &m){
    return 4*m.apply([](double x){return x*x;}).sum()*m;
}

valarray<double> Sextic::gradient(valarray<double> &m){
    return 6*m.apply([](double x){return x*x*x*x;}).sum()*m;
}

valarray<double> Octic::gradient(valarray<double> &m){
    return 8*m.apply([](double x){return x*x*x*x*x*x;}).sum()*m;
}

valarray<double> Decic::gradient(valarray<double> &m){
    return 10*m.apply([](double x){return x*x*x*x*x*x*x*x;}).sum()*m;
}

double Product::value(valarray<double> &m_neigh, valarray<double> &m_me){
    return 0;
}

double Product::diff(valarray<double> &m_neigh, valarray<double> &m_new, valarray<double> &m_old){
    return 0;
}

valarray<double> Product::gradient(valarray<double> &m_neigh, valarray<double> &m_me){
    return 0*m_neigh;
}

double J_lin_lin::value(valarray<double> &m_neigh, valarray<double> &m_me){
    return (m_neigh*m_me).sum();
}

double J_lin_lin::diff(valarray<double> &m_neigh, valarray<double> &m_new, valarray<double> &m_old){
    return (m_neigh*(m_new-m_old)).sum();
}

valarray<double> J_lin_lin::gradient(valarray<double> &m_neigh, valarray<double> &m_me){
    return m_neigh;
}

double J_cub_lin::value(valarray<double> &m_neigh, valarray<double> &m_me){
    return 0.5*(sq_mag(m_neigh)+sq_mag(m_me))*(m_neigh*m_me).sum();
}

double J_cub_lin::diff(valarray<double> &m_neigh, valarray<double> &m_new, valarray<double> &m_old){
    return 0.5*(sq_mag(m_new)*(m_new*m_neigh).sum()-sq_mag(m_old)*(m_old*m_neigh).sum()+sq_mag(m_neigh)*((m_new-m_old)*m_neigh).sum());
}

valarray<double> J_cub_lin::gradient(valarray<double> &m_neigh, valarray<double> &m_me){
    return m_me*((m_me*m_neigh).sum()+0.5*(sq_mag(m_me)+sq_mag(m_neigh)));
}

double J_qui_lin::value(valarray<double> &m_neigh, valarray<double> &m_me){
    return 0.5*(square.value(sq_mag(m_neigh))+square.value(sq_mag(m_me)))*(m_neigh*m_me).sum();
}

double J_qui_lin::diff(valarray<double> &m_neigh, valarray<double> &m_new, valarray<double> &m_old){
    return 0.5*(square.value(sq_mag(m_new))*(m_new*m_neigh).sum()-square.value(sq_mag(m_old))*(m_old*m_neigh).sum()+square.value(sq_mag(m_neigh))*((m_new-m_old)*m_neigh).sum());
}

valarray<double> J_qui_lin::gradient(valarray<double> &m_neigh, valarray<double> &m_me){
    return m_me*(2*(m_me*m_neigh).sum()+sq_mag(m_me)+sq_mag(m_neigh));
}

Atom::Atom() : mmax(10), acc(0), count(0), debug(false)
{
    m.resize(3);
    m_old.resize(3);
    gradient.resize(3);
    mabs = 1;
    phi = 0;
    theta = 0;
    set_magnitude(0, 0, 0);
    m[0] = mabs*cos(phi)*sin(theta);
    m[1] = mabs*sin(phi)*sin(theta);
    m[2] = mabs*cos(theta);
    update_flag(false);
}

void Atom::update_flag(bool ff=false){
    for(int i=0; i<2; i++)
    {
        E_uptodate[i] = ff;
        dE_uptodate[i] = ff;
    }
}

void Atom::activate_debug(){
    debug = true;
}

double Atom::get_acceptance_ratio(){
    if(count!=0)
        return acc/(double) count;
    return 0;
}

double Atom::E(int index=0, bool force_compute=false){
    if(E_uptodate[index] && !force_compute)
        return E_current[index];
    E_current[index] = 0;
    for(int i_atom=0; i_atom<int(m_n[index].size()); i_atom++)
        E_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->value(*m_n[index].at(i_atom), m);
    E_current[index] *= 0.5;
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        E_current[index] += landau_coeff[index].at(i)*landau_func[index].at(i)->value(mabs);
    if(!debug)
        E_uptodate[index] = true;
    return E_current[index];
}

double Atom::dE(int index=0, bool force_compute=false){
    if(dE_uptodate[index] && !force_compute)
        return dE_current[index];
    dE_current[index] = 0;
    count++;
    acc++;
    for(int i_atom=0; i_atom<int(m_n[index].size()); i_atom++)
        dE_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->diff(*m_n[index].at(i_atom), m, m_old);
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        dE_current[index] += landau_coeff[index].at(i)*(landau_func[index].at(i)->value(mabs)-landau_func[index].at(i)->value(mabs_old));
    if(!debug)
        dE_uptodate[index] = true;
    return dE_current[index];
}

void Atom::set_m(double mabs_new, double theta_new, double phi_new, bool diff){
    update_flag(false);
    mabs_old = mabs;
    theta_old = theta;
    phi_old = phi;
    m_old = m;
    if(diff){
        mabs += mabs_new;
        theta += theta_new;
        phi += phi_new;
    }
    else{
        mabs = mabs_new;
        theta = theta_new;
        phi = phi_new;
    }
    if(abs(mabs)>mmax)
        throw invalid_argument("Magnetic moment exploding");
    m[0] = mabs*cos(phi)*sin(theta);
    m[1] = mabs*sin(phi)*sin(theta);
    m[2] = mabs*cos(theta);
}

void Atom::update_polar_coordinates(){
    update_flag(false);
    mabs = sqrt((m*m).sum());
    theta = atan2(sqrt(m[0]*m[0]+m[1]*m[1]), m[2]);
    phi = atan2(m[1], m[2]);
}

void Atom::revoke(){
    acc--;
    set_m(mabs_old, theta_old, phi_old);
}

void Atom::set_landau_coeff(double value, int deg, int index=0){
    if(value==0)
        return;
    update_flag(false);
    if(dm==0)
        dm = 0.1;
    landau_coeff[index].push_back(value);
    switch(deg){
        case 2:
            landau_func[index].push_back(&square);
            break;
        case 4:
            landau_func[index].push_back(&quartic);
            break;
        case 6:
            landau_func[index].push_back(&sextic);
            break;
        case 8:
            landau_func[index].push_back(&octic);
            break;
        case 10:
            landau_func[index].push_back(&decic);
            break;
        default:
            throw invalid_argument("Longitudinal function not found");
    }
}

void Atom::set_heisenberg_coeff(valarray<double>* mm, double JJ, int deg=1, int index=0){
    if(JJ==0)
        return;
    update_flag(false);
    m_n[index].push_back(mm);
    heisen_coeff[index].push_back(JJ);
    if(dphi==0 && dtheta==0)
    {
        dphi = 0.1*2.0*M_PI;
        dtheta = 0.1;
    }
    switch(deg){
        case 1:
            heisen_func[index].push_back(&j_lin_lin);
            break;
        default:
            throw invalid_argument("Pairwise interaction not found");
    }
}

void Atom::clear_heisenberg_coeff(int index){
    heisen_coeff[index].clear();
    heisen_func[index].clear();
    m_n[index].clear();
}

void Atom::clear_landau_coeff(int index){
    landau_coeff[index].clear();
    landau_func[index].clear();
}

void Atom::propose_new_state(){
    double mabs_new = mabs+dm*zufall();
    double theta_new = cos(theta)+dtheta*zufall();
    double phi_new = phi+dphi*zufall();
    while(abs(theta_new)>1)
        theta_new = cos(theta)+dtheta*zufall();
    theta_new = acos(theta_new);
    set_m(mabs_new, theta_new, phi_new);
}

void Atom::set_magnitude(double ddm, double ddphi, double ddtheta)
{
    if(ddm<0 || ddphi<0 || ddtheta<0)
        throw invalid_argument( "Magnitude cannot be a negative value" );
    if(int(landau_coeff[0].size())+int(landau_coeff[1].size())==0 && ddm!=0)
        throw invalid_argument("You cannot change moment magnitude without defining Landau coefficients");
    dm = ddm;
    dphi = ddphi*2.0*M_PI;
    dtheta = ddtheta;
}

valarray<double> Atom::get_gradient(double lambda){
    valarray<double> grad(3);
    for(int index=0; index<2; index++)
    {
        if(lambda==1 && index==0)
            continue;
        for(int i_atom=0; i_atom<int(m_n[index].size()); i_atom++)
            grad -= (1-index-lambda*(1-2*index))*heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->gradient(*m_n[index].at(i_atom), m);
        for(int i=0; i<int(landau_coeff[index].size()); i++)
            grad += (1-index-lambda*(1-2*index))*landau_coeff[index].at(i)*(landau_func[index].at(i)->gradient(m));
        if(lambda==0)
            break;
    }
    return grad;
}

double Atom::get_gradient_residual(){
    return sqrt((gradient*gradient).sum());
}

double Atom::run_gradient_descent(double h, double lambda){
    valarray<double> grad(3);
    grad = get_gradient(lambda);
    valarray<double> mr = m/mabs;
    valarray<double> mphi = {-m[1], m[0], 0};
    valarray<double> mtheta = {m[0]*m[2], m[1]*m[2], -(m[0]*m[0]+m[1]*m[1])};
    if((mphi*mphi).sum()>1.0e-8)
        mphi = mphi/sqrt((mphi*mphi).sum());
    if((mtheta*mtheta).sum()>1.0e-8)
        mtheta = mtheta/sqrt((mtheta*mtheta).sum());
    set_m(-h*dm*(grad*mr).sum(), -h*dtheta*(grad*mtheta).sum(), -h*dphi*(grad*mphi).sum(), true);
    grad = m-m_old;
    double return_value = (gradient*grad).sum();
    gradient = grad;
    return return_value;
}

Atom::~Atom(){
    for(int i=0; i<2; i++)
    {
        clear_landau_coeff(i);
        clear_heisenberg_coeff(i);
    }
}

average_energy::average_energy(){ reset();}

void average_energy::add(double E_in, bool total_energy=false, int index=0)
{
    if (E_in==0)
        return;
    if (total_energy)
        E_sum[index] = E_in;
    else
        E_sum[index] += E_in;
    EE[index] += E_sum[index];
    EE_sq[index] += square.value(E_sum[index]);
    NN[index] += 1;
}

double average_energy::E_mean(int index=0){
    if(NN[index]>0)
        return EE[index]/(double)NN[index];
    else
        return 0;
}

double average_energy::E_var(int index=0){
    if(NN[index]>0)
        return (EE_sq[index]-square.value(EE[index])/(double)NN[index])/(double)NN[index];
    else
        return 0;
}

void average_energy::reset()
{
    for(int i=0; i<2; i++)
    {
        EE[i] = 0;
        NN[i] = 0;
        E_sum[i] = 0;
        EE_sq[i] = 0;
    }
}

MC::MC(): n_tot(0), debug_mode(false), kB(8.6173305e-5), lambda(-1), eta(1), E_min(0)
{
    reset();
}

void MC::set_lambda(double lambda_in)
{
    if(lambda_in<0 || lambda_in>1)
        throw invalid_argument( "Lambda must be between 0 and 1" );
    lambda = lambda_in;
}

void MC::activate_debug()
{
    debug_mode = true;
}

void MC::set_landau_coeff(vector<double> coeff, int deg, int index=0)
{
    if(int(coeff.size())!=n_tot)
        throw invalid_argument("Number of coefficients is not the same as the number of atoms");
    for(int i=0; i<n_tot; i++)
        atom[i].set_landau_coeff(coeff[i], deg, index);
}

void MC::set_heisenberg_coeff(vector<double> coeff, vector<int> me, vector<int> neigh, int deg, int index=0)
{
    if(int(coeff.size())!=int(me.size()) || int(me.size())!=int(neigh.size()))
        throw invalid_argument("Number of coefficients is not the same as the indices");
    for(int i=0; i<int(coeff.size()); i++)
        atom[me.at(i)].set_heisenberg_coeff(&atom[neigh.at(i)].m, coeff.at(i), deg, index);
}

void MC::create_atoms(int number_of_atoms)
{

    if(number_of_atoms<=0)
        throw invalid_argument("Number of atoms has to be a positive integer");
    if(n_tot!=0)
        throw invalid_argument("You cannot change the number of atoms during the simulation");
    n_tot = number_of_atoms;
    atom = new Atom[n_tot];
}

int MC::get_number_of_atoms(){
    if (n_tot==0)
        throw invalid_argument("Atoms not created yet");
    return n_tot;
}

void MC::clear_landau_coeff(int index=0)
{
    for(int i=0; i<n_tot; i++)
        atom[i].clear_landau_coeff(index);
}

void MC::clear_heisenberg_coeff(int index=0)
{
    for(int i=0; i<n_tot; i++)
        atom[i].clear_heisenberg_coeff(index);
}

bool MC::thermodynamic_integration(){
    if(lambda>=0)
        return true;
    return false;
}

bool MC::accept(int ID_rand, double kBT, double E_current){
    atom[ID_rand].propose_new_state();
    double dE = atom[ID_rand].dE();
    if(thermodynamic_integration())
        dE = (1-lambda)*dE+lambda*atom[ID_rand].dE(1);
    if(dE<=0)
        return true;
    else if(kBT==0)
        return false;
    if(bose_einstein())
    {
        if((exp((E_current-get_ground_state_energy())/kBT/n_tot)-1)/(exp(((E_current-get_ground_state_energy())/n_tot+dE)/kBT)-1)>rand()/(double)RAND_MAX)
            return true;
        else
            return false;
    }
    if(exp(-dE/(kBT*eta))>rand()/(double)RAND_MAX)
        return true;
    return false;
}

bool MC::bose_einstein()
{
    if(eta>0)
        return false;
    return true;
}

double MC::get_energy(int index=0){
    double EE=0;
    for(int i=0; i<n_tot; i++)
        EE += atom[i].E(index, true);
    return EE;
}

double MC::run_gradient_descent(int max_iter, double step_size=1, double decrement=0.001, double diff = 1.0e-8)
{
    reset();
    double residual = 0, residual_max = 0, dot_product = 0;
    for(int iter=0; iter<max_iter; iter++)
    {
        dot_product = 0;
        for(int i_atom=0; i_atom<n_tot; i_atom++)
        {
            dot_product += atom[i_atom].run_gradient_descent(step_size, lambda*thermodynamic_integration());
            residual = atom[i_atom].get_gradient_residual();
            if(i_atom==0 || residual_max<residual)
                residual_max = residual;
        }
        if(iter>0 && residual_max<diff)
            iter = max_iter;
        if(dot_product>0)
            step_size *= 1+decrement;
        else if (dot_product<0)
            step_size *= 1-decrement;
    }
    double E_min_tmp = get_energy();
    if(E_min>E_min_tmp)
        E_min = E_min_tmp;
    return residual_max;
}

double MC::get_ground_state_energy(){
    if(E_min==0)
    {
        vector<double> m = get_magnetic_moments();
        run_gradient_descent(n_tot*n_tot);
        set_magnetic_moments(m);
    }
    return E_min;
}

double MC::get_eta(){
    return eta;
}

void MC::set_eta(double eta_in){
    if(eta_in<0)
        throw invalid_argument("Invalid eta value");
    eta = eta_in;
}

void MC::run(double T_in, int number_of_iterations=1){
    double kBT = kB*T_in, dEE_tot[2], EE_tot[2];
    auto begin = std::chrono::high_resolution_clock::now();
    int ID_rand;
    for (int i=0; i<2; i++)
    {
        EE_tot[i] = get_energy(i);
        E_tot.add(EE_tot[i], true, i);
        if(!thermodynamic_integration())
            break;
    }
    for(int iter=0; iter<number_of_iterations; iter++)
    {
        for(int i=0; i<2; i++)
        {
            dEE_tot[i] = 0;
            if (debug_mode)
                EE_tot[i] = -get_energy(i);
        }
        for(int i=0; i<n_tot; i++)
        {
            MC_count++;
            ID_rand = rand()%n_tot;
            if(accept(ID_rand, kBT, EE_tot[0]+dEE_tot[0]))
            {
                acc++;
                dEE_tot[0] += atom[ID_rand].dE();
                if(thermodynamic_integration())
                    dEE_tot[1] += atom[ID_rand].dE(1);
            }
            else
                atom[ID_rand].revoke();
        }
        for(int i=0; debug_mode && i<2; i++)
        {
            EE_tot[i] += get_energy(i);
            if(abs(EE_tot[i]-dEE_tot[i])>1.0e-6*n_tot)
                throw invalid_argument( "Problem with the energy difference "+to_string(EE_tot[i])+" "+to_string(dEE_tot[i]) );
            if(!thermodynamic_integration())
                break;
        }
        E_tot.add(dEE_tot[0]);
        if(thermodynamic_integration())
            E_tot.add(dEE_tot[1], false, 1);
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - begin);
    steps_per_second = n_tot*number_of_iterations/double(duration.count())/1.0e6;
}

double MC::get_steps_per_second(){
    return (double)steps_per_second;
}

vector<double> MC::get_magnetic_moments(){
    vector<double> m(n_tot*3);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        for(int ix=0; ix<3; ix++)
            m.at(i_atom*3+ix) = atom[i_atom].m[ix];
    return m;
}

vector<double> MC::get_magnetic_gradients(){
    vector<double> m(n_tot*3);
    valarray<double> grad(3);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
    {
        grad = atom[i_atom].get_gradient(lambda);
        for(int ix=0; ix<3; ix++)
            m.at(i_atom*3+ix) = grad[ix];
    }
    return m;
}

void MC::set_magnetic_moments(vector<double> m_in)
{
    if(int(m_in.size())!=3*n_tot)
        throw invalid_argument("Length of magnetic moments not correct");
    for(int i_atom=0; i_atom<n_tot; i_atom++)
    {
        for(int ix=0; ix<3; ix++)
             atom[i_atom].m[ix] = m_in.at(i_atom*3+ix);
        atom[i_atom].update_polar_coordinates();
    }
    reset();
}

double MC::get_mean_energy(int index=0){
    return E_tot.E_mean(index);
}

double MC::get_energy_variance(int index=0){
    return E_tot.E_var(index);
}

double MC::get_acceptance_ratio(){
    if(MC_count==0)
        return 0;
    return acc/(double)MC_count;
} 

vector<double> MC::get_acceptance_ratios(){
    vector<double> v(n_tot);
    for(int i=0; i<n_tot; i++)
        v.at(i) = atom[i].get_acceptance_ratio();
    return v;
}

void MC::set_magnitude(vector<double> dm, vector<double> dphi, vector<double> dtheta)
{
    if(int(dm.size())!=int(dphi.size()) || int(dphi.size())!=int(dtheta.size()) || n_tot!=int(dm.size()))
        throw invalid_argument("Length of vectors not consistent");
    for(int i=0; i<n_tot; i++)
        atom[i].set_magnitude(dm[i], dphi[i], dtheta[i]);
}

void MC::reset()
{
    acc = 0;
    MC_count = 0;
    E_tot.reset();
}

MC::~MC()
{
    delete [] atom;
}

