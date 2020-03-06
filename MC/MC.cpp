#include "MC.h"

using namespace std;

double zufall(){
    return 1.0-2.0*(rand()/(double)RAND_MAX);
}

double square(double xxx){ return xxx*xxx; }
double quartic(double xxx){ return square(xxx)*square(xxx); }
double sextic(double xxx){ return quartic(xxx)*square(xxx); }
double octic(double xxx){ return quartic(xxx)*quartic(xxx); }
double decic(double xxx){ return sextic(xxx)*quartic(xxx); }

double J_linear(double *m_one, double *m_two, double *m_three=NULL){
    if (m_three!=NULL)
        return m_one[0]*(m_two[0]-m_three[0])+m_one[1]*(m_two[1]-m_three[1])+m_one[2]*(m_two[2]-m_three[2]);
    else
        return m_one[0]*m_two[0]+m_one[1]*m_two[1]+m_one[2]*m_two[2];
}

double J_square(double *m_one, double *m_two, double *m_three=NULL){
    if (m_three!=NULL)
        return square(J_linear(m_one, m_two))-square(J_linear(m_one, m_three));
    else
        return square(J_linear(m_one, m_two));
}

double cross_prod_sq(double *m_one, double *m_two){
    return (square(m_one[1]*m_two[2]-m_one[2]*m_two[1])
           +square(m_one[2]*m_two[0]-m_one[0]*m_two[2])
           +square(m_one[0]*m_two[1]-m_one[1]*m_two[0]));
}

double J_cross_prod(double *m_one, double *m_two, double *m_three=NULL){
    if (m_three==NULL)
        return cross_prod_sq(m_one, m_two);
    return cross_prod_sq(m_one, m_two)-cross_prod_sq(m_one, m_three);
}

Atom::Atom() : mmax(10), acc(0), count(0), debug(false)
{
    m = new double[3];
    m_old = new double[3];
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
        E_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)(m_n[index].at(i_atom), m, NULL);
    E_current[index] *= 0.5;
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        E_current[index] += landau_coeff[index].at(i)*landau_func[index].at(i)(mabs);
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
        dE_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)(m_n[index].at(i_atom), m, m_old);
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        dE_current[index] += landau_coeff[index].at(i)*(landau_func[index].at(i)(mabs)-landau_func[index].at(i)(mabs_old));
    if(!debug)
        dE_uptodate[index] = true;
    return dE_current[index];
}

void Atom::set_m(double mabs_new, double theta_new, double phi_new){
    update_flag(false);
    if(abs(mabs)>mmax)
        throw invalid_argument("Magnetic moment exploding");
    mabs_old = mabs;
    theta_old = theta;
    phi_old = phi;
    for(int i=0; i<3; i++)
        m_old[i] = m[i];
    mabs = mabs_new;
    theta = theta_new;
    phi = phi_new;
    m[0] = mabs*cos(phi)*sin(theta);
    m[1] = mabs*sin(phi)*sin(theta);
    m[2] = mabs*cos(theta);
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
            landau_func[index].push_back(square);
            break;
        case 4:
            landau_func[index].push_back(quartic);
            break;
        case 6:
            landau_func[index].push_back(sextic);
            break;
        case 8:
            landau_func[index].push_back(octic);
            break;
        case 10:
            landau_func[index].push_back(decic);
            break;
        default:
            throw invalid_argument("Longitudinal function not found");
    }
}

void Atom::set_heisenberg_coeff(double* mm, double JJ, int deg=1, int index=0){
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
            heisen_func[index].push_back(J_linear);
            break;
        case 2:
            heisen_func[index].push_back(J_square);
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
        theta_new = cos(theta)+0.2*zufall();
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

Atom::~Atom(){
    delete[] m_old;
    delete[] m;
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
    EE_sq[index] += square(E_sum[index]);
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
        return (EE_sq[index]-square(EE[index])/(double)NN[index])/(double)NN[index];
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

MC::MC(): n_tot(0), debug_mode(false), kB(8.6173305e-5), lambda(-1), eta(1)
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
        atom[me.at(i)].set_heisenberg_coeff(atom[neigh.at(i)].m, coeff.at(i), deg, index);
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

bool MC::accept(int ID_rand, double kBT){
    atom[ID_rand].propose_new_state();
    double dE = atom[ID_rand].dE();
    if(thermodynamic_integration())
        dE = (1-lambda)*dE+lambda*atom[ID_rand].dE(1);
    if(dE<=0)
        return true;
    else if(kBT==0)
        return false;
    if(preparing_qmc())
    {
        double E_old = atom[ID_rand].E()-dE;
        if((exp(E_old/kBT)-1)/(exp((E_old+dE)/kBT)-1)>rand()/(double)RAND_MAX)
            return true;
        else
            return false;
    }
    if(exp(-dE/(kBT*eta))>rand()/(double)RAND_MAX)
        return true;
    return false;
}

bool MC::preparing_qmc()
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

void MC::prepare_qmc(double T_in, int number_of_iterations){
    if(thermodynamic_integration())
        throw invalid_argument("QMC+Thermodynamic integration now allowed");
    eta = 0;
    vector<double> m = get_magnetic_moments();
    run(0, number_of_iterations);
    double E_current = get_energy();
    set_magnetic_moments(m);
    run(T_in, number_of_iterations);
    reset();
    run(T_in, number_of_iterations);
    eta = (E_tot.E_mean()-E_current)/n_tot/(kB*T_in);
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
    clock_t begin = clock();
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
            if(accept(ID_rand, kBT))
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
    steps_per_second = n_tot*number_of_iterations/(double)(clock()-begin)*CLOCKS_PER_SEC;
}

double MC::get_steps_per_second(){
    return (double)steps_per_second;
}

vector<double> MC::get_magnetic_moments(){
    vector<double> m;
    m.resize(n_tot*3);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        for(int ix=0; ix<3; ix++)
            m.at(i_atom*3+ix) = atom[i_atom].m[ix];
    return m;
}

void MC::set_magnetic_moments(vector<double> m_in)
{
    if(int(m_in.size())!=3*n_tot)
        throw invalid_argument("Length of magnetic moments not correct");
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        for(int ix=0; ix<3; ix++)
             atom[i_atom].m[ix] = m_in.at(i_atom*3+ix);
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

