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

Atom::Atom() : A(0), B(0), n_neigh(0), n_max(0), acc(0), count(0), E_uptodate(false), debug(false)
{
    m = new double[3];
    m_old = new double[3];
    mabs = 1;
    phi = 0;
    theta = 0;
    set_magnitude(0.1, 0.1, 0.1);
    m[0] = mabs*cos(phi)*sin(theta);
    m[1] = mabs*sin(phi)*sin(theta);
    m[2] = mabs*cos(theta);
}

void Atom::set_num_neighbors(int num_neighbors=96){
    n_max = num_neighbors;
    m_n = new double*[n_max];
    J = new double[n_max];
    for(int i=0; i<n_max; i++)
        J[i] = 0;
}

void Atom::activate_debug(){
    debug = true;
}

float Atom::acceptance_ratio(){
    if(count!=0)
        return acc/(float) count;
    return 0;
}

double Atom::E(bool force_compute=false){
    if(E_uptodate && !force_compute)
        return E_current;
    E_current = 0;
    for(int i_atom=0; i_atom<n_neigh; i_atom++)
        E_current -= J[i_atom]*(m_n[i_atom][0]*m[0]
                                +m_n[i_atom][1]*m[1]
                                +m_n[i_atom][2]*m[2]);
    E_current *= 0.5;
    E_current += A*square(mabs)+B*quartic(mabs);
    if(!debug)
        E_uptodate = true;
    return E_current;
}

double Atom::dE(){
    double EE = 0;
    count++;
    acc++;
    for(int i_atom=0; i_atom<n_neigh; i_atom++)
        EE -= J[i_atom]*(m_n[i_atom][0]*(m[0]-m_old[0])
                         +m_n[i_atom][1]*(m[1]-m_old[1])
                         +m_n[i_atom][2]*(m[2]-m_old[2]));
    return A*(square(mabs)-square(mabs_old))+B*(quartic(mabs)-quartic(mabs_old))+EE;
}

void Atom::set_m(double mabs_new, double theta_new, double phi_new){
    E_uptodate = false;
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

void Atom::set_AB(double A_in, double B_in){
    E_uptodate = false;
    if(A!=0 || B!=0)
        cout<<"WARNING: A and B have already been set"<<endl;
    if(A_in==0 || B_in==0)
        cout<<"WARNING: Setting A=0 or B=0"<<endl;
    if(B_in<0)
        cout<<"WARNING: Negative B value will make it diverge"<<endl;
    A = A_in;
    B = B_in;
}

void Atom::set_neighbor(double* mm, double JJ){
    if(n_neigh>=n_max)
        throw invalid_argument("number of neighbors too small");
    E_uptodate = false;
    m_n[n_neigh] = mm;
    J[n_neigh] = JJ;
    n_neigh++;
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
    delete m_old;
    delete J;
    delete m_n;
    delete m;
}

average_energy::average_energy(){ reset();}

void average_energy::add(double E_in, bool total_energy=false)
{
    if (E_in==0)
        return;
    if (total_energy)
        E_sum = E_in;
    else
        E_sum += E_in;
    EE += E_sum;
    EE_sq += square(E_sum);
    NN += 1;
}

double average_energy::E_mean(){
    if(NN>0)
        return EE/(double)NN;
    else
        return 0;
}

double average_energy::E_var(){
    if(NN>0)
        return (EE_sq-square(EE)/(double)NN)/(double)NN;
    else
        return 0;
}

void average_energy::reset()
{
    EE = 0;
    NN = 0;
    E_sum = 0;
    EE_sq = 0;
}

MC::MC(): thermodynamic_integration(false), kB(8.6173305e-5)
{
    reset();
}

void MC::set_lambda(double lambda_in)
{
    if(lambda_in<0 || lambda_in>1)
        throw invalid_argument( "Lambda must be between 0 and 1" );
    thermodynamic_integration = true;
    lambda = lambda_in;
}

void MC::activate_debug()
{
    debug_mode = true;
}

void MC::create_atoms(int num_neigh, vector<double> A, vector<double> B, vector<int> me, vector<int> neigh, vector<double> J)
{

    N_tot = int(A.size());
    atom = new Atom[N_tot];
    for(int i=0; i<N_tot; i++)
    {
        atom[i].set_num_neighbors(num_neigh);
        atom[i].set_AB(A[i], B[i]);
    }
    for(int i=0; i<int(J.size()); i++)
        atom[me.at(i)].set_neighbor(atom[neigh.at(i)].m, J.at(i));
}

double MC::run(double T_in, int number_of_iterations=1){
    double kBT = kB*T_in, dEE_tot = 0, dE, EE_tot=0;
    clock_t begin = clock();
    int ID_rand;
    for(int i=0; i<N_tot; i++)
        EE_tot += atom[i].E(true);
    E_tot.add(EE_tot, true);
    for(int iter=0; iter<number_of_iterations; iter++)
    {
        dEE_tot = 0, EE_tot = 0;
        for(int i=0; debug_mode && i<N_tot; i++)
            EE_tot -= atom[i].E(true);
        for(int i=0; i<N_tot; i++)
        {
            MC_count++;
            ID_rand = rand()%N_tot;
            atom[ID_rand].propose_new_state();
            dE = atom[ID_rand].dE();
            if(dE<0 || (kBT>0)*exp(-dE/kBT)>rand()/(double)RAND_MAX)
            {
                acc++;
                dEE_tot += dE;
            }
            else
                atom[ID_rand].revoke();
        }
        if(debug_mode)
        {
            for(int i=0; i<N_tot; i++)
                EE_tot += atom[i].E(true);
            if(abs(EE_tot-dEE_tot)>1.0e-6*N_tot)
                throw invalid_argument( "Problem with the energy difference "+to_string(EE_tot)+" "+to_string(dEE_tot) );
        }
        E_tot.add(dEE_tot);
    }
    steps_per_second = N_tot*number_of_iterations/(double)(clock()-begin)*CLOCKS_PER_SEC;
    return dEE_tot;
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

double MC::get_energy(){
    double EE=0;
    for(int i=0; i<N_tot; i++)
        EE += atom[i].E(true);
    return EE;
}

double MC::get_mean_energy(){
    return E_tot.E_mean();
}

double MC::get_energy_variance(){
    return E_tot.E_var();
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
    delete atom;
}

