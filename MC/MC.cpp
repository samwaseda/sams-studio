#include "MC.h"

double RandomNumberFactory::uniform(bool symmetric, double max_value){
    if (symmetric)
        return max_value*(1.0-2.0*((double)rand()/(double)RAND_MAX));
    else
        return max_value*((double)rand()/(double)RAND_MAX);
}

valarray<double> RandomNumberFactory::on_sphere(int size){
    valarray<double> m_new(size);
    for(int i=0; i<size; i++)
        m_new[i] = uniform();
    m_new *= uniform()/sqrt((m_new*m_new).sum());
    return m_new;
}

double RandomNumberFactory::normal(){
    return distribution(generator);
}

valarray<double> RandomNumberFactory::n_on_sphere(int size){
    return normal()*on_sphere(size);
}

double m_norm(valarray<double> mm){
    return sqrt((mm*mm).sum());
}

valarray<double> m_cross(valarray<double>& m_one, valarray<double> m_two){
    return m_one.cshift(1)*m_two.cshift(2)-m_one.cshift(2)*m_two.cshift(1);
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
    return (neigh.m*(me.m-me.m_tmp)).sum();
}

valarray<double> J_lin_lin::gradient(Atom &neigh, Atom &me){
    return neigh.m;
}

double J_cub_lin::value(Atom &neigh, Atom &me){
    return 0.5*(neigh.get_magnitude(2)+me.get_magnitude(2))*(neigh.m*me.m).sum();
}

double J_cub_lin::diff(Atom &neigh, Atom &me){
    return 0.5*(me.get_magnitude(2)*(me.m*neigh.m).sum()
                -me.get_magnitude(2, true)*(me.m_tmp*neigh.m).sum()
                +neigh.get_magnitude(2)*((me.m-me.m_tmp)*neigh.m).sum());
}

valarray<double> J_cub_lin::gradient(Atom &neigh, Atom &me){
    return (me.m*(me.m*neigh.m).sum()+0.5*neigh.m*(me.get_magnitude(2)+neigh.get_magnitude(2)));
}

double J_qui_lin::value(Atom &neigh, Atom &me){
    return 0.5*(neigh.get_magnitude(4)+me.get_magnitude(4))*(neigh.m*me.m).sum();
}

double J_qui_lin::diff(Atom &neigh, Atom &me){
    return 0.5*(me.get_magnitude(4)*(me.m*neigh.m).sum()
                -me.get_magnitude(4, true)*(me.m_tmp*neigh.m).sum()
                +neigh.get_magnitude(4)*((me.m-me.m_tmp)*neigh.m).sum());
}

valarray<double> J_qui_lin::gradient(Atom &neigh, Atom &me){
    return (2*me.m*me.get_magnitude(2)*(me.m*neigh.m).sum()
            +0.5*neigh.m*(me.get_magnitude(4)+neigh.get_magnitude(4)));
}

Atom::Atom() : mabs(1), mmax(100), acc(0), count(0), debug(false)
{
    m.resize(3);
    m_tmp.resize(3);
    gradient.resize(3);
    m[0] = mabs;
    set_magnitude(0, 0);
    update_flag(false);
}

void Atom::update_flag(bool ff){
    up_to_date.E.assign(2, ff);
    up_to_date.dE.assign(2, ff);
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
        return power(mabs_tmp, exponent);
    else
        return power(mabs, exponent);
}

double Atom::E(int index, bool force_compute){
    if(up_to_date.E.at(index) && !force_compute)
        return E_current[index];
    E_current[index] = 0;
    for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
        E_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->value(
            *neigh[index].at(i_atom), *this);
    E_current[index] *= 0.5;
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        E_current[index] += landau_coeff[index].at(i)*landau_func[index].at(i)->value(mabs);
    if(!debug)
        up_to_date.E.at(index) = true;
    return E_current[index];
}

double Atom::dE(int index, bool force_compute){
    if(up_to_date.dE.at(index) && !force_compute)
        return dE_current[index];
    dE_current[index] = 0;
    count++;
    acc++;
    for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
        dE_current[index] -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->diff(
            *neigh[index].at(i_atom), *this);
    for(int i=0; i<int(landau_coeff[index].size()); i++)
        dE_current[index] += landau_coeff[index].at(i)*(
            landau_func[index].at(i)->value(mabs)-landau_func[index].at(i)->value(mabs_tmp));
    if(!debug)
        up_to_date.dE.at(index) = true;
    return dE_current[index];
}

void Atom::set_m(valarray<double> m_new, bool diff){
    update_flag(false);
    mabs_tmp = mabs;
    m_tmp = m;
    if(diff && abs(dm-1)+abs(dphi-1)==0)
        m += m_new;
    else if(diff){
        m += dphi*m_new;
        m *= m_norm(m_tmp+dm*m_new)/sqrt((m*m).sum());
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
        throw invalid_argument(
            "mabs: "+to_string(sqrt((m*m).sum()))+" vs. "+to_string(abs(mabs)));
}

void Atom::revoke(){
    acc--;
    update_flag(false);
    m = m_tmp;
    mabs = mabs_tmp;
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
    set_m(rand_generator.on_sphere(3), true);
    if(flip && rand()%2==1)
        m *= -1;
}

valarray<double> Atom::delta_m(){
    return m-m_tmp;
}

void Atom::rescale_magnitude(double rescale_m, double rescale_phi){
    if (rescale_m==1 && dm>0)
        dm = 1;
    else
        dm *= rescale_m;
    if (rescale_phi==1 && dphi>0)
        dphi = 1;
    else
        dphi *= rescale_phi;
}

void Atom::set_magnitude(double ddm, double ddphi, bool flip_in)
{
    if(ddm<0 || ddphi<0)
        throw invalid_argument( "Magnitude cannot be a negative value" );
    if(int(landau_coeff[0].size())+int(landau_coeff[1].size())==0 && ddm!=0)
        throw invalid_argument(
            "You cannot change moment magnitude without defining Landau coefficients");
    dm = ddm;
    dphi = ddphi;
    flip = flip_in;
}

valarray<double> Atom::get_gradient(double lambda){
    valarray<double> grad(3), grad_tot(3);
    if (lambda<0)
        lambda = 0;
    for(int index=0; index<2; index++)
    {
        if(lambda==1 && index==0)
            continue;
        for(int i_atom=0; i_atom<int(neigh[index].size()); i_atom++)
            grad -= heisen_coeff[index].at(i_atom)*heisen_func[index].at(i_atom)->gradient(
                *neigh[index].at(i_atom), *this);
        for(int i=0; i<int(landau_coeff[index].size()); i++)
            grad += landau_coeff[index].at(i)*(landau_func[index].at(i)->gradient(m));
        if(lambda==0) {
            grad_tot = grad;
            break;
        }
        else
            grad_tot += (1-index-lambda*(1-2*index))*grad;
    }
    return grad_tot;
}

double Atom::get_gradient_residual(){
    return sqrt((gradient*gradient).sum());
}

double Atom::run_gradient_descent(double h, double lambda){
    valarray<double> grad(3);
    grad = -h*get_gradient(lambda);
    set_m(grad, true);
    grad = m-m_tmp;
    double return_value = (gradient*grad).sum();
    gradient = grad;
    return return_value;
}

void Atom::calc_spin_dynamics(
    double damping_parameter, double delta_t, double mu_s, double lambda){
    valarray<double> h_field = -get_gradient(lambda);
    valarray<double> h_stochastic = mu_s*rand_generator.n_on_sphere();
    if (dm==0) {
        m_tmp = m_cross(m, h_field+h_stochastic);
        m_tmp -= damping_parameter*m_cross(m, m_cross(m, h_field));
    }
    else {
        m_tmp = m_cross(m, h_field);
        m_tmp += h_stochastic;
        m_tmp += damping_parameter*h_field;
    }
    m_tmp *= delta_t/constants.hbar;
}

void Atom::update_spin_dynamics(){
    set_m(m_tmp, true);
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
        E_sum.at(index) = E_in;
    else
        E_sum.at(index) += E_in;
    EE.at(index) += E_sum.at(index);
    EE_sq.at(index) += square.value(E_sum.at(index));
    NN += 1;
}

double average_energy::E_mean(int index=0){
    return EE.at(index)/(double)NN;
}

double average_energy::E_var(int index=0){
    return (EE_sq.at(index)-square.value(EE.at(index))/(double)NN)/(double)NN;
}

void average_energy::reset()
{
    EE.assign(2, 0);
    E_sum.assign(2, 0);
    EE_sq.assign(2, 0);
    NN = 0;
}

MC::MC(): n_tot(0), lambda(-1), debug_mode(false), spin_dynamics_flag(false)
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

void MC::set_heisenberg_coeff(
    vector<double> coeff, vector<int> me, vector<int> neigh, int deg, int index)
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
    selectable_id.resize(n_tot);
    for(int i_atom=0; i_atom<n_tot; i_atom++)
        selectable_id.at(i_atom) = i_atom;
}

void MC::select_id(vector<int> select_id_in)
{
    if(int(select_id_in.size())>n_tot)
        throw invalid_argument("select_id longer than the number of atoms");
    selectable_id = select_id_in;
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

void MC::run_spin_dynamics(double kBT, int threads){
    double mu_s = sqrt(2*constants.damping_parameter*constants.hbar*kBT/constants.delta_t);
    #pragma omp parallel num_threads(threads)
    {
        #pragma omp for
        for (int i=0; i<n_tot; i++)
            atom[i].calc_spin_dynamics(
                constants.damping_parameter, constants.delta_t, mu_s, lambda);
        #pragma omp for
        for (int i=0; i<n_tot; i++)
            atom[i].update_spin_dynamics();
    }
    E_tot.add(get_energy(0), true);
    if(thermodynamic_integration())
        E_tot.add(get_energy(1), true);
    reset_magnetization();
}

void MC::run_mc(double kBT){
    int id_rand;
    double EE_tot[2];
    vector<double> dEE_tot (2,0);
    for(int i=0; i<2 && debug_mode; i++)
        EE_tot[i] = -get_energy(i);
    for(int i=0; i<n_tot; i++)
    {
        MC_count++;
        id_rand = selectable_id.at(rand()%selectable_id.size());
        atom[id_rand].propose_new_state();
        double dEE = atom[id_rand].dE();
        if (meta.initialized)
            dEE += meta.get_biased_energy(
                m_norm(magnetization+atom[id_rand].delta_m()/n_tot),
                sqrt((magnetization*magnetization).sum()));
        update_magnetization(id_rand);
        if(thermodynamic_integration())
            dEE = (1-lambda)*dEE+lambda*atom[id_rand].dE(1);
        if(metropolis(kBT, dEE))
        {
            acc++;
            dEE_tot.at(0) += atom[id_rand].dE();
            if(thermodynamic_integration())
                dEE_tot.at(1) += atom[id_rand].dE(1);
        }
        else
        {
            update_magnetization(id_rand, true);
            atom[id_rand].revoke();
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
    E_tot.add(dEE_tot[0]);
    if(thermodynamic_integration())
        E_tot.add(dEE_tot[1], false, 1);
}

bool MC::metropolis(double kBT, double energy_difference){
    if(energy_difference<=0)
        return true;
    else if(kBT==0)
        return false;
    if(exp(-energy_difference/kBT)>(double)(rand())/(double)RAND_MAX)
        return true;
    return false;
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
    return residual_max;
}

void MC::run(double T_in, int number_of_iterations, int threads){
    double kBT = constants.kB*T_in;
    vector<double> dEE_tot;
    auto begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<2; i++)
    {
        E_tot.add(get_energy(i), true, i);
        if(!thermodynamic_integration())
            break;
    }
    for(int iter=0; iter<number_of_iterations; iter++)
    {
        if (spin_dynamics_flag)
            run_spin_dynamics(kBT, threads);
        else
            run_mc(kBT);
        magnetization_hist.push_back(m_norm(magnetization));
        if (meta.initialized)
            meta.append_value(magnetization_hist.back());
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

void MC::switch_spin_dynamics(bool on, double damping_parameter, double delta_t, bool rescale_mag)
{
    if (spin_dynamics_flag != on && rescale_mag){
        for (int i=0; i<n_tot && on; i++)
            atom[i].rescale_magnitude(1*(damping_parameter>0), 1);
        for (int i=0; i<n_tot && !on; i++)
            atom[i].rescale_magnitude(0.1*(damping_parameter>0), 0.1);
    }
    spin_dynamics_flag = on;
    if (damping_parameter < 0)
        throw invalid_argument("damping_parameter must be a positive float");
    constants.damping_parameter = damping_parameter;
    if (delta_t<=0)
        throw invalid_argument("delta_t must be a positive float");
    constants.delta_t = delta_t;
}

void MC::set_metadynamics(double aa, double bb, double cc, int dd, double ee, int ff)
{
    meta.set_metadynamics(aa, bb, cc, dd, ee ,ff);
}

void MC::update_magnetization(int mc_id, bool backward)
{
    if (backward)
        magnetization -= atom[mc_id].delta_m()/n_tot;
    else
        magnetization += atom[mc_id].delta_m()/n_tot;
}

vector<double> MC::get_magnetization(){
    return magnetization_hist;
}

vector<double> MC::get_histogram(int derivative){
    return meta.get_histogram(magnetization_hist, derivative);
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

double Metadynamics::get_biased_gradient(double m){
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    if (!use_derivative)
        throw invalid_argument("use_derivative not chosen");
    double mass = max_range/hist.size();
    if (m>=max_range)
        return 0;
    return hist.at(int(m*0.5/mass));
}

double Metadynamics::get_biased_energy(double m_new, double m_tmp){
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    double mass = max_range/hist.size();
    if (m_new>=max_range || m_tmp>=max_range)
        return 0;
    if (use_derivative)
        return (m_new-m_tmp)*hist.at(int((m_new+m_tmp)*0.5/mass));
    return hist.at(int(m_new/mass))-hist.at(int(m_tmp/mass));
}

double Metadynamics::gauss_exp(double m, int i)
{
    return energy_increment*exp(-square.value(m-max_range*(i+0.5)/hist.size())/denominator);
}

void Metadynamics::append_value(double m)
{
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    for (int i=i_min(m); !use_derivative && i<i_max(m); i++)
        hist.at(i) += gauss_exp(m, i);
    for (int i=i_min(m); use_derivative && i<i_max(m); i++)
        hist.at(i) += 2*(m-max_range*(0.5+i)/hist.size())/denominator*gauss_exp(m, i);
}

int Metadynamics::i_min(double m){return max(0, int(hist.size()*(m-cutoff)/max_range));}
int Metadynamics::i_max(double m){
    return min(int(hist.size()), int(hist.size()*(m+cutoff)/max_range));
}

vector<double> Metadynamics::get_histogram(vector<double>& magnetization, int derivative)
{
    if (!initialized)
        throw invalid_argument("metadynamics not initialized yet");
    if (derivative!=0 && !use_derivative)
        throw invalid_argument("derivative can be taken only if use_derivative is activated");
    vector<double> m_range(hist.size());
    for (int i=0; i<int(m_range.size()); i++)
        m_range.at(i) = max_range*i/m_range.size();
    if (!use_derivative || derivative!=0)
    {
        m_range.insert( m_range.end(), hist.begin(), hist.end() );
        return m_range;
    }
    vector<double> h_tmp (hist.size(), 0);
    for (auto m: magnetization)
        for (int i=i_min(m); i<i_max(m); i++)
            h_tmp.at(i) += gauss_exp(m, i);
    m_range.insert( m_range.end(), h_tmp.begin(), h_tmp.end() );
    return m_range;
}

