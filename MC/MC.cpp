#include "MC.h"

double zufall(){
    return 1.0-2.0*((double)rand()/(double)RAND_MAX);
}

double m_norm(valarray<double> mm){
    return sqrt((mm*mm).sum());
}

double power(double x, int exponent){
    switch(exponent){
        case 1:
            return x;
        case 2:
            return square.value(x);
        case 4:
            return quartic.value(x);
        case 6:
            return sextic.value(x);
        case 8:
            return octic.value(x);
        case 10:
            return decic.value(x);
        default:
            return exp(exponent*log(x));
    }
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

double Product::value(Atom &neigh, Atom &me){
    return 0;
}

double Product::diff(Atom &neigh, Atom &me){
    return 0;
}

valarray<double> Product::gradient(Atom &neigh, Atom &me){
    return 0*neigh.m;
}

double J_lin_lin::value(Atom &neigh, Atom &me){
    return (neigh.m*me.m).sum();
}

double J_lin_lin::diff(Atom &neigh, Atom &me){
    return (neigh.m*(me.m-me.m_old)).sum();

}

valarray<double> J_lin_lin::gradient(Atom &neigh, Atom &me){
    return neigh.m;
}

double J_cub_lin::value(Atom &neigh, Atom &me){
    return 0.5*(neigh.get_magnitude(2)+me.get_magnitude(2))*(neigh.m*me.m).sum();
}

double J_cub_lin::diff(Atom &neigh, Atom &me){
    return 0.5*(me.get_magnitude(2)*(me.m*neigh.m).sum()-me.get_magnitude(2, true)*(me.m_old*neigh.m).sum()+neigh.get_magnitude(2)*((me.m-me.m_old)*neigh.m).sum());
}

valarray<double> J_cub_lin::gradient(Atom &neigh, Atom &me){
    return (me.m*(me.m*neigh.m).sum()+0.5*neigh.m*(me.get_magnitude(2)+neigh.get_magnitude(2)));
}

double J_qui_lin::value(Atom &neigh, Atom &me){
    return 0.5*(neigh.get_magnitude(4)+me.get_magnitude(4))*(neigh.m*me.m).sum();
}

double J_qui_lin::diff(Atom &neigh, Atom &me){
    return 0.5*(me.get_magnitude(4)*(me.m*neigh.m).sum()-me.get_magnitude(4, true)*(me.m_old*neigh.m).sum()+neigh.get_magnitude(4)*((me.m-me.m_old)*neigh.m).sum());
}

valarray<double> J_qui_lin::gradient(Atom &neigh, Atom &me){
    return (2*me.m*me.get_magnitude(2)*(me.m*neigh.m).sum()+0.5*neigh.m*(me.get_magnitude(4)+neigh.get_magnitude(4)));
}

double J_cross_forward::value(Atom &neigh, Atom &me){
    return neigh.get_magnitude(2)-square.value((neigh.m*me.m).sum())/me.get_magnitude(2);
}

double J_cross_forward::diff(Atom &neigh, Atom &me){
    return square.value((neigh.m*me.m).sum())/me.get_magnitude(2)-square.value((neigh.m*me.m_old).sum())/me.get_magnitude(2, true);
}

valarray<double> J_cross_forward::gradient(Atom &neigh, Atom &me){
    return 2*(-me.m*(me.m*neigh.m).sum()/me.get_magnitude(2)+neigh.m*square.value((me.m*neigh.m).sum())/me.get_magnitude(4));
}

double J_cross_backward::value(Atom &neigh, Atom &me){
    return me.get_magnitude(2)-square.value((neigh.m*me.m).sum())/neigh.get_magnitude(2);
}

double J_cross_backward::diff(Atom &neigh, Atom &me){
    return me.get_magnitude(2)-square.value((neigh.m*me.m).sum())/neigh.get_magnitude(2)-(me.get_magnitude(2, true)-square.value((neigh.m_old*me.m).sum())/neigh.get_magnitude(2));
}

valarray<double> J_cross_backward::gradient(Atom &neigh, Atom &me){
    return 2*(me.m*me.get_magnitude()-neigh.m*(me.m*neigh.m).sum()/neigh.get_magnitude(2));
}

Atom::Atom() : mmax(100), acc(0), count(0), debug(false)
{
    m.resize(3);
    m_old.resize(3);
    gradient.resize(3);
    mabs = 1;
    m[0] = mabs;
    set_magnitude(0, 0);
    update_flag(false);
}

void Atom::update_flag(bool ff){
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

double Atom::get_magnitude(int exponent, bool old)
{
    if(old)
        return power(mabs_old, exponent);
    else
        return power(mabs, exponent);
}

double Atom::E(int index, bool force_compute){
    if(E_uptodate[index] && !force_compute)
        return E_current[index];
    E_current[index] = 0;
    for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
        E_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->value(*neigh[index].at(i_atom), *this);
    E_current[index] *= 0.5;
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        E_current[index] += landau_coeff[index].at(i)*landau_func[index].at(i)->value(mabs);
    if(!debug)
        E_uptodate[index] = true;
    return E_current[index];
}

double Atom::dE(int index, bool force_compute){
    if(dE_uptodate[index] && !force_compute)
        return dE_current[index];
    dE_current[index] = 0;
    count++;
    acc++;
    for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
        dE_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->diff(*neigh[index].at(i_atom), *this);
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        dE_current[index] += landau_coeff[index].at(i)*(landau_func[index].at(i)->value(mabs)-landau_func[index].at(i)->value(mabs_old));
    if(!debug)
        dE_uptodate[index] = true;
    return dE_current[index];
}

void Atom::set_m(valarray<double>& m_new, bool diff){
    update_flag(false);
    mabs_old = mabs;
    m_old = m;
    if(diff){
        m += dphi*m_new;
        m *= sqrt(((m_old+dm*m_new)*(m_old+dm*m_new)).sum()/(m*m).sum());
    }
    else{
        m = m_new;
    }
    mabs = sqrt((m*m).sum());
    if(mabs>mmax)
        throw invalid_argument("Magnetic moment exploding");
}

void Atom::check_consistency() {
    if ( abs(sqrt((m*m).sum())-abs(mabs))>1.0e-8 )
        throw invalid_argument("mabs: "+to_string(sqrt((m*m).sum()))+" vs. "+to_string(abs(mabs)));
}

void Atom::revoke(){
    acc--;
    update_flag(false);
    m = m_old;
    mabs = mabs_old;
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
            return landau_func[index].push_back(&square);
        case 4:
            return landau_func[index].push_back(&quartic);
        case 6:
            return landau_func[index].push_back(&sextic);
        case 8:
            return landau_func[index].push_back(&octic);
        case 10:
            return landau_func[index].push_back(&decic);
        default:
            throw invalid_argument("Longitudinal function not found");
    }
}

void Atom::set_heisenberg_coeff(Atom &neigh_in, double JJ, int deg, int index){
    if(JJ==0)
        return;
    update_flag(false);
    neigh[index].push_back(&neigh_in);
    heisen_coeff[index].push_back(JJ);
    if(dphi==0)
        dphi = 0.1;
    if (deg!=1 && dm==0)
        dm = 0.1;
    switch(deg){
        case 1:
            return heisen_func[index].push_back(&j_lin_lin);
        case 3:
            return heisen_func[index].push_back(&j_cub_lin);
        case 5:
            return heisen_func[index].push_back(&j_qui_lin);
        case 11:
            return heisen_func[index].push_back(&j_cross_forward);
        case 13:
            return heisen_func[index].push_back(&j_cross_backward);
        default:
            throw invalid_argument("Pairwise interaction not found");
    }
}

void Atom::clear_heisenberg_coeff(int index){
    heisen_coeff[index].clear();
    heisen_func[index].clear();
    neigh[index].clear();
}

void Atom::clear_landau_coeff(int index){
    landau_coeff[index].clear();
    landau_func[index].clear();
}

void Atom::propose_new_state(){
    valarray<double> m_new(3);
    for(int i=0; i<3; i++)
        m_new[i] = zufall();
    m_new *= zufall()/sqrt((m_new*m_new).sum());
    set_m(m_new, true);
    if(flip && rand()%2==1)
        m *= -1;
}

valarray<double> Atom::delta_m(){
    return m-m_old;
}

void Atom::set_magnitude(double ddm, double ddphi, bool flip_in)
{
    if(ddm<0 || ddphi<0)
        throw invalid_argument( "Magnitude cannot be a negative value" );
    if(int(landau_coeff[0].size())+int(landau_coeff[1].size())==0 && ddm!=0)
        throw invalid_argument("You cannot change moment magnitude without defining Landau coefficients");
    dm = ddm;
    dphi = ddphi;
    flip = flip_in;
}

valarray<double> Atom::get_gradient(double lambda){
    valarray<double> grad(3);
    for(int index=0; index<2; index++)
    {
        if(lambda==1 && index==0)
            continue;
        for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
            grad -= (1-index-lambda*(1-2*index))*heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->gradient(*neigh[index].at(i_atom), *this);
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
    grad = -get_gradient(lambda);
    set_m(grad, true);
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

void average_energy::add(double E_in, bool total_energy, int index)
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

MC::MC(): n_tot(0), debug_mode(false), kB(8.617333262145e-5), lambda(-1), eta(1), E_min(0)
{
    srand (time(NULL));
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

void MC::set_landau_coeff(vector<double> coeff, int deg, int index)
{
    if(int(coeff.size())!=n_tot)
        throw invalid_argument("Number of coefficients is not the same as the number of atoms");
    for(int i=0; i<n_tot; i++)
        atom[i].set_landau_coeff(coeff[i], deg, index);
}

void MC::set_heisenberg_coeff(vector<double> coeff, vector<int> me, vector<int> neigh, int deg, int index)
{
    if(int(coeff.size())!=int(me.size()) || int(me.size())!=int(neigh.size()))
        throw invalid_argument("Number of coefficients is not the same as the indices");
    for(int i=0; i<int(coeff.size()); i++)
        atom[me.at(i)].set_heisenberg_coeff(atom[neigh.at(i)], coeff.at(i), deg, index);
}

void MC::create_atoms(int number_of_atoms)
{

    if(number_of_atoms<=0)
        throw invalid_argument("Number of atoms has to be a positive integer");
    if(n_tot!=0)
        throw invalid_argument("You cannot change the number of atoms during the simulation");
    n_tot = number_of_atoms;
    atom = new Atom[n_tot];
    selectable_ID.resize(n_tot);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        selectable_ID.at(i_atom) = i_atom;
}

void MC::select_ID(vector<int> select_ID_in)
{
    if(int(select_ID_in.size())>n_tot)
        throw invalid_argument("select_ID longer than the number of atoms");
    selectable_ID = select_ID_in;
}

int MC::get_number_of_atoms(){
    if (n_tot==0)
        throw invalid_argument("Atoms not created yet");
    return n_tot;
}

void MC::clear_landau_coeff(int index)
{
    for(int i=0; i<n_tot; i++)
        atom[i].clear_landau_coeff(index);
}

void MC::clear_heisenberg_coeff(int index)
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
    double dEE = atom[ID_rand].dE();
    if (meta.initialized)
        dEE += meta.get_biased_energy(
            m_norm(magnetization+atom[ID_rand].delta_m()/n_tot), m_norm(magnetization));
    update_magnetization(ID_rand);
    if(thermodynamic_integration())
        dEE = (1-lambda)*dEE+lambda*atom[ID_rand].dE(1);
    if(dEE<=0)
        return true;
    else if(kBT==0)
        return false;
    if(bose_einstein())
    {
        if((exp((E_current-get_ground_state_energy())/kBT/n_tot)-1)/(exp(((E_current-get_ground_state_energy())/n_tot+dEE)/kBT)-1)>rand()/(double)RAND_MAX)
            return true;
        else
            return false;
    }
    if(exp(-dEE/(kBT*eta))>(double)(rand())/(double)RAND_MAX)
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

double MC::run_gradient_descent(int max_iter, double step_size, double decrement, double diff)
{
    reset();
    double residual = 0, residual_max = 0, dot_product = 0;
    for(int iter=0; iter<max_iter; iter++)
    {
        dot_product = 0;
        for(int i_atom=0; i_atom<n_tot; i_atom++)
        {
            dot_product += atom[i_atom].run_gradient_descent(
                step_size, lambda*thermodynamic_integration());
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

void MC::run_debug(){
    vector<double> m_tmp(3);
    double E = get_energy();
    for(int i=0; i<n_tot; i++)
    {
        for(int j=0; j<3; j++)
            m_tmp.at(j) = atom[i].m[j];
        if(abs(atom[i].E(0, true)-atom[i].E())>1.0e-8)
            throw invalid_argument(
                "force compute not working: "
                +to_string(atom[i].E(0, true))+" "+to_string(atom[i].E()));
        double E_before = get_energy();
        accept(i, 0, 0);
        if(abs(get_energy()-E_before-atom[i].dE())>1.0e-8)
            throw invalid_argument(
                "energy difference wrong: "
                +to_string(atom[i].E()-E_before)+" "+to_string(atom[i].dE()));
        atom[i].revoke();
        for(int j=0; j<3; j++)
            if(abs(atom[i].m[j]-m_tmp.at(j))>1.0e-8)
                throw invalid_argument("revoke not worked properly for magnetic moments");
    }
    if(abs(E-get_energy())>1.0e-8)
        throw invalid_argument("revoke not worked properly for total energy");
    for(int i=0; i<n_tot; i++)
        atom[i].check_consistency();
    reset();
    run(1000, 100);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        for(int ix=0; ix<3; ix++)
            magnetization[ix] -= atom[i_atom].m[ix]/n_tot;
    if (abs(magnetization).sum()>1.0e-6)
        throw invalid_argument("Magnetization wrong: "+to_string(abs(magnetization).sum()));
    reset();
}

void MC::run(double T_in, int number_of_iterations){
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
            ID_rand = selectable_ID.at(rand()%selectable_ID.size());
            if(accept(ID_rand, kBT, EE_tot[0]+dEE_tot[0]))
            {
                acc++;
                dEE_tot[0] += atom[ID_rand].dE();
                if(thermodynamic_integration())
                    dEE_tot[1] += atom[ID_rand].dE(1);
            }
            else
            {
                update_magnetization(ID_rand, true);
                atom[ID_rand].revoke();
            }
        }
        for(int i=0; debug_mode && i<2; i++)
        {
            EE_tot[i] += get_energy(i);
            if(abs(EE_tot[i]-dEE_tot[i])>1.0e-6*n_tot)
                throw invalid_argument(
                    "Problem with the energy difference "
                    +to_string(EE_tot[i])+" "+to_string(dEE_tot[i]) );
            if(!thermodynamic_integration())
                break;
            for(int i=0; i<n_tot; i++)
                atom[i].check_consistency();
        }
        magnetization_hist.push_back(m_norm(magnetization));
        if (meta.initialized)
            meta.append_value(m_norm(magnetization));
        E_tot.add(dEE_tot[0]);
        if(thermodynamic_integration())
            E_tot.add(dEE_tot[1], false, 1);
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - begin);
    steps_per_second = n_tot*number_of_iterations/double(duration.count())*1.0e6;
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
        valarray<double> mm(&(m_in.at(i_atom*3)), 3);
        atom[i_atom].set_m(mm);
    }
    reset();
}

double MC::get_mean_energy(int index){
    return E_tot.E_mean(index);
}

double MC::get_energy_variance(int index){
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

void MC::set_magnitude(vector<double> dm, vector<double> dphi, vector<int> flip)
{
    if(int(dm.size())!=int(dphi.size()) || n_tot!=int(dm.size()))
        throw invalid_argument("Length of vectors not consistent");
    for(int i=0; i<n_tot; i++)
        atom[i].set_magnitude(dm[i], dphi[i], (flip[i]>0));
}

void MC::set_metadynamics(double aa, double bb, double cc, int dd, double ee, int ff)
{
    meta.set_metadynamics(aa, bb, cc, dd, ee ,ff);
}

void MC::update_magnetization(int mc_id, bool backward)
{
    if (backward)
        magnetization -= (atom[mc_id].m-atom[mc_id].m_old)/n_tot;
    else
        magnetization += (atom[mc_id].m-atom[mc_id].m_old)/n_tot;
}

vector<double> MC::get_magnetization(){
    return magnetization_hist;
}

void MC::reset()
{
    acc = 0;
    MC_count = 0;
    E_tot.reset();
    reset_magnetization();
    magnetization_hist.clear();
}

void MC::reset_magnetization()
{
    magnetization.resize(3);
    magnetization = 0*magnetization;
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        magnetization += atom[i_atom].m;
    magnetization /= n_tot;
}

MC::~MC()
{
    delete [] atom;
}

Metadynamics::Metadynamics() : initialized(false){};

void Metadynamics::set_metadynamics(
    double max_range_in,
    double energy_increment_in,
    double length_scale_in,
    int bins,
    double cutoff_in,
    int derivative)
{
    if (max_range_in<=0)
        throw invalid_argument("max_range must be a positive float");
    if (energy_increment_in<0)
        throw invalid_argument("energy_increment must be a positive float");
    if (length_scale_in<=0)
        throw invalid_argument("length_scale must be a positive float");
    if (bins<=0)
        throw invalid_argument("bins must be a positive integer");
    if (cutoff_in<=0)
        throw invalid_argument("cutoff_in must be a positive float");
    initialized = true;
    max_range = max_range_in;
    energy_increment = energy_increment_in;
    denominator = length_scale_in*length_scale_in*2;
    hist.assign(bins, 0);
    cutoff = cutoff_in*length_scale_in;
    use_derivative = false;
    if (derivative != 0)
        use_derivative = true;
}

double Metadynamics::get_biased_energy(double m_new, double m_old){
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    double mass = max_range/hist.size();
    if (m_new>=max_range || m_old>=max_range)
        return 0;
    if (use_derivative)
        return (m_new-m_old)*hist.at(int((m_new+m_old)*0.5/mass));
    return hist.at(int(m_new/mass))-hist.at(int(m_old/mass));
}

void Metadynamics::append_value(double m){
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    int i_min = max(0, int(hist.size()*(m-cutoff)/max_range));
    int i_max = min(int(hist.size()), int(hist.size()*(m+cutoff)/max_range));
    for (int i=i_min && !use_derivative; i<i_max; i++)
        hist.at(i) += energy_increment*exp(
            -square.value(m-max_range*i/hist.size())/denominator);
    for (int i=i_min && use_derivative; i<i_max; i++)
    {
        double m_tmp = m-max_range*i/hist.size();
        hist.at(i) += m_tmp*energy_increment*exp(-square.value(m_tmp)/denominator);
    }
}

