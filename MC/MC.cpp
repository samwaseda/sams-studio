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

Atom::Atom() : n_max(0), n_neigh(0), A(0), B(0), E_uptodate(false), acc(0), count(0), type("Fe"), N_spread(4), m_0(0)
{
    x = new float[3];
    m = new double[3];
    m_old = new double[3];
    mabs = 1;
    phi = 0;
    theta = 0;
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

bool Atom::set_type(string type_in, bool restart=false){
    E_uptodate = false;
    if(type==type_in && !restart)
        cout<<"New type the same as the old one"<<endl;
    else
    {
        type = type_in;
        return true;
    }
    return false;
}

string Atom::get_type(){
    return type;
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
    E_uptodate = true;
    return E_current;
}

double Atom::E_harmonic(){
    return slope*square(mabs-m_0);
}

double Atom::dE_harmonic(){
    return slope*(square(mabs-m_0)-square(mabs_old-m_0));
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

void Atom::flip_z(){
    set_m(mabs, -theta, -phi);
}

int Atom::nx(int ix){
    if(ix<0 || ix>3)
        cout<<"ERROR: x[3]"<<endl;
    return round(N_spread*x[ix]);
}

void Atom::modify_AB(double A_in, double B_in){
    E_uptodate = false;
    if(A==0 && B==0)
        cout<<"WARNING: A and B seem not to have been set"<<endl;
    A += A_in;
    B += B_in;
}

void Atom::set_AB(double A_in, double B_in){
    E_uptodate = false;
    if(A!=0 || B!=0)
        cout<<"WARNING: A and B have already been set"<<endl;
    if(A_in==0 || B_in==0)
        cout<<"WARNING: Setting A=0 or B=0"<<endl;
    A = A_in;
    B = B_in;
    slope = 2.0*abs(A);
    if (A<0)
        m_0 = sqrt(-0.5*A/B);
}

bool Atom::check_saddle(){
    if(nx(0)%2==1 || nx(1)%2==1 || nx(2)%2==1)
        return true;
    return false;
}

void Atom::set_neighbor(double* mm, double JJ){
    if(n_neigh>=n_max)
    {
        cout<<"ERROR: n_max too small: "<<n_neigh<<" "<<n_max<<endl;
        exit(EXIT_FAILURE);
    }
    E_uptodate = false;
    m_n[n_neigh] = mm;
    J[n_neigh] = JJ;
    n_neigh++;
}

int Atom::get_num_neighbors(){
    return n_neigh;
}

bool Atom::modify_neighbor(double* mm, double JJ){
    for (int i=0; i<n_neigh; i++)
        if(m_n[i] == mm)
        {
            J[i] += JJ;
            E_uptodate = false;
            return true;
        }
    cout<<"ERROR: neighbor not found"<<endl;
    return false;
}

void Atom::propose_new_state(){
    double mabs_new = abs(mabs+0.1*zufall());
    double theta_new = cos(theta)+0.2*zufall();
    double phi_new = phi+0.4*M_PI*zufall();
    while(abs(theta_new)>1)
        theta_new = cos(theta)+0.2*zufall();
    theta_new = acos(theta_new);
    set_m(mabs_new, theta_new, phi_new);
}

average_energy::average_energy(): kappa(1.0-1.0e-3){ reset();}

void average_energy::add(double E_in, bool total_energy=false)
{
    if (total_energy)
        E_sum = E_in;
    else if (E_sum!=0)
        E_sum += E_in;
    else
        return;
    EE = E_sum+EE*kappa;
    NN += 1;
}

double average_energy::E(){
    if(NN>0)
        return EE*(1-kappa)/(1-exp(log(kappa)*NN));
    else
        return 0;
}

void average_energy::reset()
{
    EE = 0;
    NN = 0;
    E_sum = 0;
}

Energy::Energy(double lambda_in=1, bool debug_mode_in=true): acc(0), MC_count(0), kB(8.6173305e-5), lambda(lambda_in)
{
    cout<<"Starting initialization"<<endl;
    lambda = lambda_in;
    ausgabe.open("MC.log", ios::out);
    ausgabe<<"# Starting a Metropolis Monte Carlo simulation"<<endl;

    debug_mode = debug_mode_in;
    if(debug_mode)
        ausgabe<<"## Debug mode activated"<<endl;
}

void Energy::create_atoms(int n_tot, int num_neigh, vector<double> A, vector<double> B)
{
    N_tot = n_tot;
    atom = new Atom[N_tot];
    for(int i=0; i<N_tot; i++)
    {
        atom[i].set_num_neighbors(num_neigh);
        atom[i].set_AB(A[i], B[i]);
    }

    double EE = 0, E = 0, E_max, E_min, E_harm = 0;
    for(int i=0; i<N_tot; i++)
    {
        E += atom[i].E(true);
        EE += square(atom[i].E());
        E_harm += atom[i].E_harmonic();
        if(i==0 || E_max<atom[i].E())
            E_max = atom[i].E();
        if(i==0 || E_min>atom[i].E())
            E_min = atom[i].E();
    }
    E_tot.add(E-E_harm, true);
    ausgabe<<"# Initial total energy: "<<E<<" per atom: "<<E/N_tot<<"+-"<<sqrt(EE*N_tot-square(E))/(double)N_tot<<" Emax: "<<E_max<<" Emin: "<<E_min<<endl;
    if(E/N_tot<-1)
        cout<<"WARNING: Energy per atom "<<E/N_tot<<endl;
    begin = clock();
}

double Energy::MC(double T_in){
    double kBT = kB*T_in, dEE_tot = 0, dE, EE_tot=0;
    int ID_rand;
    for(int i=0; debug_mode && i<N_tot; i++)
    {
        if(lambda<1)
            EE_tot -= atom[i].E(true)-atom[i].E_harmonic();
        else
            EE_tot -= atom[i].E(true);
    }
    for(int i=0; i<N_tot; i++)
    {
        MC_count++;
        ID_rand = rand()%N_tot;
        atom[ID_rand].propose_new_state();
        dE = atom[ID_rand].dE();
        if(lambda<1)
            dE = lambda*dE+(1.0-lambda)*atom[ID_rand].dE_harmonic();
        if(dE<0 || (kBT>0)*exp(-dE/kBT)>rand()/(double)RAND_MAX)
        {
            acc++;
            dEE_tot += atom[ID_rand].dE()-atom[ID_rand].dE_harmonic();
        }
        else
            atom[ID_rand].revoke();
    }
    if(debug_mode)
    {
        for(int i=0; i<N_tot; i++)
            EE_tot += atom[i].E(true)-atom[i].E_harmonic();
        if(abs(EE_tot-dEE_tot)>1.0e-6*N_tot)
        {
            cout<<"ERROR: Problem with the energy difference: "<<EE_tot<<" "<<dEE_tot<<endl;
            exit(EXIT_FAILURE);
        }
    }
    E_tot.add(dEE_tot);
    return dEE_tot;
}

double Energy::output(string custom_text="", bool config=false){
    double m_ave[3], mm_ave[3], E=0, m_min=0, m_max=0, EE=0, E_min=0, E_max=0, m_tmp, E_harm=0;
    for(int i=0; i<3; i++)
    {
        m_ave[i] = 0;
        mm_ave[i] = 0;
    }
    for(int i=0; i<N_tot; i++)
    {
        for(int j=0; j<3; j++)
        {
            m_ave[j] += atom[i].m[j];
            mm_ave[j] += atom[i].m[j]*atom[i].m[j];
        }
        m_tmp = sqrt(square(atom[i].m[0])+square(atom[i].m[1])+square(atom[i].m[2]));
        if(m_min>m_tmp || i==0)
            m_min = m_tmp;
        if(m_max<m_tmp || i==0)
            m_max = m_tmp;
        E += atom[i].E(true);
        E_harm += atom[i].E_harmonic();
        EE += square(atom[i].E());
        if(E_min>atom[i].E() || i==0)
            E_min = atom[i].E();
        if(E_max>atom[i].E() || i==0)
            E_max = atom[i].E();
    }
    E_tot.add(E-E_harm, true);
    if(custom_text.size()>0)
        custom_text += " ";
    ausgabe<<int((double)(clock()-begin)/CLOCKS_PER_SEC+0.5)<<" "<<custom_text
        <<m_ave[0]/N_tot<<" "<<m_ave[1]/N_tot<<" "<<m_ave[2]/N_tot<<" "
        <<sqrt(mm_ave[0]*N_tot-m_ave[0]*m_ave[0])/N_tot<<" "<<sqrt(mm_ave[1]*N_tot-m_ave[1]*m_ave[1])/N_tot<<" "<<sqrt(mm_ave[2]*N_tot-m_ave[2]*m_ave[2])/N_tot<<" "
        <<(double)acc/(double)MC_count<<" "<<fixed<<setprecision(6)<<E<<" "<<fixed<<setprecision(6)<<E_harm<<" "<<fixed<<setprecision(6)<<E_tot.E()<<" "<<sqrt(EE*N_tot-E*E)/N_tot<<endl;
    m_tmp = m_ave[0]*m_ave[0]+m_ave[1]*m_ave[1]+m_ave[2]*m_ave[2];
    cout<<setw(7)<<int((double)(clock()-begin)/CLOCKS_PER_SEC+0.5)
        <<setw(10)<<sqrt(m_tmp)/N_tot
        <<setw(11)<<sqrt((mm_ave[0]+mm_ave[1]+mm_ave[2])*N_tot-m_tmp)/N_tot
        <<setw(10)<<(double)acc/(double)MC_count<<setw(10)<<E_tot.E()
        <<setw(10)<<E_harm
        <<" "<<E/N_tot<<" "<<sqrt(EE*N_tot-E*E)/N_tot<<endl;
    if(config)
    {
        config_ausgabe.open("config.dat", ios::out);
        for(int i=0; i<N_tot; i++)
            config_ausgabe<<i<<" "<<atom[i].get_type()<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<" "
                <<atom[i].m[0]<<" "<<atom[i].m[1]<<" "<<atom[i].m[2]<<" "
                <<atom[i].E()<<" "<<atom[i].acceptance_ratio()<<endl;
        config_ausgabe.close();
    }
    return E;
}

void Energy::E_min(){
    cout<<"Starting to calculate the energy minimum"<<endl;
    //for(int i=0; i<N_tot; i++)
    //	m[i] = 2.2;
    double dE = 1;
    do
    {
        dE = 0;
        for(int i=0; i<1000; i++)
            dE += MC(0);
        output();
    } while(abs(dE)>1.0e-4);
}

void Energy::reset()
{
    acc = 0;
    MC_count = 0;
    E_tot.reset();
}

