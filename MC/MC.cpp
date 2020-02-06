#include "MC.h"

using namespace std;

double zufall(){
    return 1.0-2.0*(rand()/(double)RAND_MAX);
}

double square(double xxx){
    return xxx*xxx;
}

double quartic(double xxx){
    return square(xxx)*square(xxx);
}

Atom::Atom() : acc(0), count(0), debug(false)
{
    m = new double[3];
    m_old = new double[3];
    A.resize(2, 0);
    B.resize(2, 0);
    mabs = 1;
    phi = 0;
    theta = 0;
    set_magnitude(0.1, 0.1, 0.1);
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

float Atom::get_acceptance_ratio(){
    if(count!=0)
        return acc/(float) count;
    return 0;
}

double Atom::E(bool force_compute=false, int index=0){
    if(E_uptodate[index] && !force_compute)
        return E_current;
    E_current = 0;
    for(int i_atom=0; i_atom<int(m_n[index].size()); i_atom++)
        E_current -= J[index].at(i_atom)*(m_n[index].at(i_atom)[0]*m[0]
                                          +m_n[index].at(i_atom)[1]*m[1]
                                          +m_n[index].at(i_atom)[2]*m[2]);
    E_current *= 0.5;
    E_current += A[index]*square(mabs)+B[index]*quartic(mabs);
    if(!debug)
        E_uptodate[index] = true;
    return E_current;
}

double Atom::dE(bool force_compute=false, int index=0){
    if(dE_uptodate[index] && !force_compute)
        return dE_current;
    dE_current = 0;
    count++;
    acc++;
    for(int i_atom=0; i_atom<int(m_n[index].size()); i_atom++)
        dE_current -= J[index].at(i_atom)*(m_n[index].at(i_atom)[0]*(m[0]-m_old[0])
                                           +m_n[index].at(i_atom)[1]*(m[1]-m_old[1])
                                           +m_n[index].at(i_atom)[2]*(m[2]-m_old[2]));
    dE_current += A[index]*(square(mabs)-square(mabs_old))+B[index]*(quartic(mabs)-quartic(mabs_old));
    if(!debug)
        dE_uptodate[index] = true;
    return dE_current;
}

void Atom::set_m(double mabs_new, double theta_new, double phi_new){
    update_flag(false);
    if(abs(mabs)>5)
        cout<<"WARNING: Magnetic moment value = "<<mabs<<endl;
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

void Atom::set_AB(double A_in, double B_in, int index=0){
    update_flag(false);
    if(A[index]!=0 || B[index]!=0)
        cout<<"WARNING: A and B have already been set"<<endl;
    if(A_in==0 || B_in==0)
        cout<<"WARNING: Setting A=0 or B=0"<<endl;
    if(B_in<0)
        cout<<"WARNING: Negative B value will make it diverge"<<endl;
    A[index] = A_in;
    B[index] = B_in;
}

void Atom::set_neighbor(double* mm, double JJ, int index=0){
    update_flag(false);
    m_n[index].push_back(mm);
    J[index].push_back(JJ);
}

void Atom::propose_new_state(){
    double mabs_new = abs(mabs+dm*zufall());
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

MC::MC(): thermodynamic_integration_flag(0), kB(8.6173305e-5)
{
    reset();
}

void MC::set_lambda(double lambda_in)
{
    if(lambda_in<0 || lambda_in>1)
        throw invalid_argument( "Lambda must be between 0 and 1" );
    if(thermodynamic_integration_flag%2==0)
        thermodynamic_integration_flag++;
    lambda = lambda_in;
}

void MC::activate_debug()
{
    debug_mode = true;
}

void MC::create_atoms(vector<double> A, vector<double> B, vector<int> me, vector<int> neigh, vector<double> J)
{

    N_tot = int(A.size());
    atom = new Atom[N_tot];
    for(int i=0; i<N_tot; i++)
        atom[i].set_AB(A[i], B[i]);
    for(int i=0; i<int(J.size()); i++)
        atom[me.at(i)].set_neighbor(atom[neigh.at(i)].m, J.at(i));
}

void MC::append_parameters(vector<double> A, vector<double> B, vector<int> me, vector<int> neigh, vector<double> J)
{
    for(int i=0; i<N_tot; i++)
        atom[i].set_AB(A[i], B[i], 1);
    for(int i=0; i<int(J.size()); i++)
        atom[me.at(i)].set_neighbor(atom[neigh.at(i)].m, J.at(i), 1);
    if(thermodynamic_integration_flag<2)
        thermodynamic_integration_flag += 2;
}

bool MC::thermodynamic_integration(){
    switch (thermodynamic_integration_flag){
        case 0:
            return false;
        case 1:
            throw invalid_argument("Parameters not set for lambda=1");
        case 2:
            throw invalid_argument("Lambda parameter not set");
        case 3:
            return true;
        default:
            throw invalid_argument("Something went wrong");
    }
}

void MC::run(double T_in, int number_of_iterations=1){
    double kBT = kB*T_in, dEE_tot[2], dE, EE_tot[2];
    clock_t begin = clock();
    int ID_rand;
    EE_tot[0] = 0;
    EE_tot[1] = 0;
    for(int i=0; i<N_tot; i++)
        EE_tot[0] += atom[i].E(true);
    E_tot.add(EE_tot[0], true, 0);
    if(thermodynamic_integration())
        for(int i=0; i<N_tot; i++)
            EE_tot[1] += atom[i].E(true, 1);
    E_tot.add(EE_tot[1], true, 1);
    for(int iter=0; iter<number_of_iterations; iter++)
    {
        for(int i=0; i<2; i++)
        {
            dEE_tot[i] = 0;
            EE_tot[i] = 0;
        }
        for(int i=0; debug_mode && i<N_tot; i++)
            EE_tot[0] -= atom[i].E(true);
        for(int i=0; i<N_tot; i++)
        {
            MC_count++;
            ID_rand = rand()%N_tot;
            atom[ID_rand].propose_new_state();
            dE = atom[ID_rand].dE();
            if(thermodynamic_integration())
                dE = (1-lambda)*dE+lambda*atom[ID_rand].dE(false, 1);
            if(dE<0 || (kBT>0)*exp(-dE/kBT)>rand()/(double)RAND_MAX)
            {
                acc++;
                dEE_tot[0] += atom[ID_rand].dE(false, 0);
                if(thermodynamic_integration())
                    dEE_tot[1] += atom[ID_rand].dE(false, 1);
            }
            else
                atom[ID_rand].revoke();
        }
        if(debug_mode)
        {
            for(int i=0; i<N_tot; i++)
                EE_tot[0] += atom[i].E(true);
            if(abs(EE_tot[0]-dEE_tot[0])>1.0e-6*N_tot)
                throw invalid_argument( "Problem with the energy difference "+to_string(EE_tot[0])+" "+to_string(dEE_tot[0]) );
        }
        E_tot.add(dEE_tot[0]);
        if(thermodynamic_integration())
            E_tot.add(dEE_tot[1], false, 1);
    }
    steps_per_second = N_tot*number_of_iterations/(double)(clock()-begin)*CLOCKS_PER_SEC;
}

double MC::get_steps_per_second(){
    return (double)steps_per_second;
}

vector<double> MC::get_magnetic_moments(){
    vector<double> m;
    m.resize(N_tot*3);
    for(int i_atom=0; i_atom<N_tot; i_atom++)
        for(int ix=0; ix<3; ix++)
            m.at(i_atom*3+ix) = atom[i_atom].m[ix];
    return m;
}

double MC::get_energy(int index=0){
    double EE=0;
    for(int i=0; i<N_tot; i++)
        EE += atom[i].E(true, index);
    return EE;
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

void MC::set_magnitude(double dm, double dphi, double dtheta)
{
    for(int i=0; i<N_tot; i++)
        atom[i].set_magnitude(dm, dphi, dtheta);
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

