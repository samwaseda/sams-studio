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

Coeff::Coeff(string input_file, fstream &ausgabe) : shell_max(6)
{
	initialize_coeff();
	cout<<"Starting to read parameter file"<<endl;
	string line;
	fstream eingabe;
	eingabe.open(input_file, ios::in);
	int dist_in, pot_in;
	string neigh_type_in, elem_one_in, elem_two_in, elem_in;
	double wert_in;
	char Koeff;
	while(eingabe>>Koeff)
	{
		if(Koeff=='J')
		{
			eingabe>>elem_one_in>>elem_two_in>>dist_in>>wert_in;
			ausgabe<<"## Jij parameter: Jij "<<elem_one_in<<"-"<<elem_two_in<<" r_"<<dist_in<<"  = "<<wert_in<<endl;
			if(JJ[El(elem_one_in)+El(elem_two_in)][Def(elem_one_in)+Def(elem_two_in)][dist_in] != 0)
				ausgabe<<"## WARNING: Jij["<<elem_one_in<<","<<elem_two_in<<"]["<<dist_in<<"] set twice"<<endl;
			JJ[El(elem_one_in)+El(elem_two_in)][Def(elem_one_in)+Def(elem_two_in)][dist_in] = wert_in;
			if(dist_in+1>shell_max)
				shell_max = dist_in+1;
		}
		else if (Koeff=='A')
		{
			eingabe>>elem_in>>wert_in;
			ausgabe<<"## Ai parameter: "<<elem_in<<" = "<<wert_in<<endl;
			if(AA[El(elem_in)][Def(elem_in)] != 0)
				ausgabe<<"## WARNING: Ai["<<elem_in<<"] set twice"<<endl;
			AA[El(elem_in)][Def(elem_in)] = wert_in;
		}
		else if (Koeff=='B')
		{
			eingabe>>elem_in>>wert_in;
			ausgabe<<"## Bi parameter: "<<elem_in<<" = "<<wert_in<<endl;
			if(BB[El(elem_in)][Def(elem_in)] != 0)
				ausgabe<<"## WARNING: Bi["<<elem_in<<"] set twice"<<endl;
			BB[El(elem_in)][Def(elem_in)] = wert_in;
		}
		else
			ausgabe<<"## Unknown parameter: "<<Koeff<<endl;
		getline(eingabe, line);
	}
	eingabe.close();
}

void Coeff::initialize_coeff(){
	JJ = new double ** [2*Elmax()];
	for (int i=0; i<2*Elmax(); i++)
		JJ[i] = new double *[2*Defmax()];
	for (int i=0; i<2*Elmax(); i++)
		for (int j=0; j<2*Defmax(); j++)
			JJ[i][j] = new double [shell_max];
	AA = new double *[Elmax()];
	BB = new double *[Elmax()];
	for (int i=0; i<Elmax(); i++){
		AA[i] = new double [Defmax()];
		BB[i] = new double [Defmax()];
	}
	for(int i=0; i<2*Elmax(); i++)
		for(int j=0; j<2*Defmax(); j++)
			for (int k=0; k<shell_max; k++)
				JJ[i][j][k] = 0;
	for(int i=0; i<Elmax(); i++)
		for(int j=0; j<Defmax(); j++){
			AA[i][j] = 0;
			BB[i][j] = 0;
		}
	shell_max = 0;
}

double Coeff::A(string elem, string def="b"){
	return AA[El(elem)][Def(def)];
}

double Coeff::B(string elem, string def="b"){
	return BB[El(elem)][Def(def)];
}

double Coeff::J(string elem_one, string elem_two, int shell_in, string def_one="b", string def_two="b"){
	return JJ[El(elem_one)+El(elem_two)][Def(def_one)+Def(def_two)][shell_in];
}

int Coeff::shell(){ return shell_max; }

int Coeff::Elmax(){return int(pow(2,ID_type.size()))-1;}

int Coeff::Defmax(){return int(pow(2,defect_type.size()))-1;}

int Coeff::El(string str){
	int pos = str.find('_');
	if (pos<=str.size())
		str = str.substr(0, pos);
	auto search = ID_type.find(str);
	if (search != ID_type.end())
		return search->second;
	else
	{
		cout<<"Element not found: "<<str<<endl;
		exit(EXIT_FAILURE);
	}
}

int Coeff::Def(string str){
	int pos = str.find('_');
	if (pos>str.size())
		return 0;
	else
	{
		str = str.substr(pos+1);
		auto search = defect_type.find(str);
		if (search != defect_type.end())
			return search->second;
		else
		{
			cout<<"Defect not found: "<<str<<endl;
			exit(EXIT_FAILURE);
		}
	}
}

Shell::Shell(int shell_max_in, string bravais) : shell_max(shell_max_in+1){
	float dist_fcc_sb[] = {0, 3.0/8.0, 5.0/8.0, 7.0/8.0, 9.0/8.0, 11.0/8.0, 13.0/8.0, 15.0/8.0, 17.0/8.0};
	float dist_fcc_bb[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
	float dist_bcc_sb[] = {0, 11.0/16.0, 19.0/16.0, 27.0/16.0, 35.0/16.0, 43.0/16.0, 51.0/16.0, 59.0/16.0, 67.0/16.0};
	float dist_bcc_bb[] = {0, 3.0/4.0, 1.0, 2.0, 11.0/4.0, 3.0, 4.0, 19.0/4.0, 5.0};
	if(shell_max>8)
	{
		cout<<"Max shell count cannot surpass 8"<<endl;
		abort();
	}
	dist_sb = new float[shell_max];
	dist_bb = new float[shell_max];
	for(int i=0; i<shell_max; i++)
	{
		if(bravais=="bcc")
		{
			dist_sb[i] = 0.5*(sqrt(dist_bcc_sb[i])+sqrt(dist_bcc_sb[i+1]));
			dist_bb[i] = 0.5*(sqrt(dist_bcc_bb[i])+sqrt(dist_bcc_bb[i+1]));
		}
		else
		{
			dist_sb[i] = 0.5*(sqrt(dist_fcc_sb[i])+sqrt(dist_fcc_sb[i+1]));
			dist_bb[i] = 0.5*(sqrt(dist_fcc_bb[i])+sqrt(dist_fcc_bb[i+1]));
		}
	}
}

int Shell::get_shell(bool saddle, float distance)
{
	if(distance>5.1)
		return -1;
	if(!saddle)
		dist_tmp = dist_bb;
	else
		dist_tmp = dist_sb;
	if(distance<dist_tmp[0])
	{
		cout<<"Distance too small: "<<distance<<endl;
		abort();
	}
	for(int i=1; i<shell_max; i++)
		if(distance<dist_tmp[i])
			return i-1;
	return -1;
}

Shell::~Shell(){
	delete[] dist_sb, delete[] dist_bb;
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

Structure::Structure(string bravais, int N_in, bool debug_mode_in) : N_v(0), N_s(0)
{
	if(bravais!="bcc" && bravais!="fcc")
	{
		cout<<"Bravais lattice not recognized"<<endl;
		exit(EXIT_FAILURE);
	}
	N_rep = N_in;
	fstream restart;
	ausgabe.open("structure.dat", ios::out);
	restart.open("config.dat", ios::in);
	if(restart)
	{
		ausgabe<<"## Restarting from previous config.dat"<<endl;
		N_rep = count(istreambuf_iterator<char>(restart), istreambuf_iterator<char>(), '\n');
		ausgabe<<"## "<<N_rep<<" atoms detected"<<endl;
	}
	restart.close();

	if(N_rep>=100)
	{
		if(bravais=="bcc")
			N_rep = round(exp(log(0.5*double(N_rep))/3.0));
		else if(bravais=="fcc")
			N_rep = round(exp(log(0.25*double(N_rep))/3.0));
	}

	N_tot = N_rep*N_rep*N_rep*2;
	if(bravais=="fcc")
		N_tot += N_rep*N_rep*N_rep*2;

	Atom* atom = new Atom[N_tot];
	if(restart)
		reload_lattice(atom, "config.dat");
	else
		create_lattice(atom, bravais);
}

Atom* Structure::get_structure()
{
	return atom;
}

void Structure::set_coeff(string input_file, int num_neighbors, bool debug_mode)
{
	for(int ind=0; ind<N_tot; ind++)
		atom[ind].set_num_neighbors(num_neighbors);
	Coeff cff = Coeff(input_file, ausgabe);
	shell_max = cff.shell();

	ausgabe<<"## Number of shells: "<<shell_max<<" = "<<cff.shell()<<endl;

	vector<int> count_shell(shell_max, 0);
	int A_modif_count = 0;
	double A_modif_value = 0;
	double B_modif_value = 0;
	vector<int> J_modif_count(shell_max, 0);
	vector<double> J_modif_value(shell_max, 0);

	kdtree tree(num_neighbors);
	initialize_tree((&tree));

	Shell shell(shell_max, bravais);

	if(debug_mode)
	{
		int count_vacancy = 0;
		for(int i=0; i<N_tot; i++)
			if(check_vacancy(i, bravais, &tree, &shell))
				count_vacancy++;
		cout<<"Vacancy count: "<<count_vacancy<<endl;
	}

	for(int i=0; i<N_tot; i++)
	{
		atom[i].set_AB(cff.A(atom[i].get_type()), cff.B(atom[i].get_type()));
		if(N_s>0 && shell.get_shell(true, dist(&(atom[i]), &(atom[N_tot-1])))==0)
		{
			A_modif_value += cff.A(atom[i].get_type());
			B_modif_value += cff.B(atom[i].get_type());
			A_modif_count++;
		}
		if(N_v!=0 && check_vacancy(i, bravais, &tree, &shell))
		{
			atom[i].modify_AB(cff.A(atom[i].get_type(), "Fe_v"), cff.B(atom[i].get_type(), "Fe_v"));
			A_modif_value += cff.A(atom[i].get_type(), "Fe_v");
			B_modif_value += cff.B(atom[i].get_type(), "Fe_v");
			A_modif_count++;
		}
		for(int j=0; j<tree.num_neighbors; j++)
		{
			int ind_nn = tree.get_index(i, j);
			//int i_shell = shell.get_shell(atom[i].check_saddle()+atom[ind_nn].check_saddle(), dist(&(atom[i]), &(atom[ind_nn])));
			int i_shell = shell.get_shell(atom[i].check_saddle()+atom[ind_nn].check_saddle(), tree.get_distance(i,j));
			if(i_shell<0)
				continue;
			count_shell.at(i_shell)++;
			atom[i].set_neighbor(atom[ind_nn].m, cff.J(atom[i].get_type(), atom[ind_nn].get_type(), i_shell));
			if (N_s>0 && atom[i].check_saddle()+atom[ind_nn].check_saddle())
			{
				J_modif_count.at(i_shell) += 1;
				J_modif_value.at(i_shell) += cff.J(atom[i].get_type(), atom[ind_nn].get_type(), i_shell);
			}
			if(N_v!=0 && check_vacancy(i, bravais, &tree, &shell) && check_vacancy(ind_nn, bravais, &tree, &shell) && !atom[ind_nn].check_saddle())
			{
				double vv = cff.J(atom[i].get_type(), atom[ind_nn].get_type(), i_shell, "Fe_v", "Fe_v");
				if(!atom[i].modify_neighbor(atom[ind_nn].m, vv))
					exit(EXIT_FAILURE);
				J_modif_count.at(i_shell) += 1;
				J_modif_value.at(i_shell) += vv;
			}
		}
	}

	ausgabe<<"## Shell count statistics:"<<endl<<"##";
	for(int i=0; i<cff.shell(); i++)
		ausgabe<<" "<<count_shell.at(i);
	ausgabe<<endl<<"## Per atom (rounded up):"<<endl<<"##";
	for(int i=0; i<cff.shell(); i++)
		ausgabe<<" "<<round(count_shell.at(i)/(float)N_tot);
	ausgabe<<endl<<"## In perfect crystal:"<<endl<<"##";
	for(int i=0; i<cff.shell(); i++)
		ausgabe<<" "<<N_tot*round(count_shell.at(i)/(float)N_tot);
	ausgabe<<endl<<"## A,B modification count and values: "<<A_modif_count<<" "<<A_modif_value<<" "<<B_modif_value<<endl;
	ausgabe<<"## J modification count:";
	for(int i=0; i<shell_max; i++)
		ausgabe<<" "<<J_modif_count.at(i);
	ausgabe<<endl<<"## J modification values:";
	for(int i=0; i<shell_max; i++)
		ausgabe<<" "<<J_modif_value.at(i);
	ausgabe<<endl;

	ausgabe<<"# Bravais lattice: "<<bravais<<endl;
	ausgabe<<"# Number of lattices in each direction: "<<N_rep<<endl;
	ausgabe<<"# Total number of atoms: "<<N_tot<<endl;
	//ausgabe<<"# Number of vacancies: "<<N_v<<endl;
	//ausgabe<<"# Number of saddle point atoms: "<<N_s<<endl;
	//ausgabe<<"# Number or fraction of Mn atoms: "<<N_Mn<<endl;
}

int Structure::get_number()
{
	return N_tot;
}

void Structure::set_solute(double N_Mn_in)
{
	int N_Mn = 0;
	if(N_Mn_in != 0)
	{
		if(N_Mn_in < 1)
			N_Mn = int(N_Mn_in*N_tot);
		else
			N_Mn = int(N_Mn_in);
	}
	for(int i_Mn=0; i_Mn<N_Mn; i_Mn++)
	{
		int ID_rand = rand()%N_tot;
		if(!atom[ID_rand].set_type("Mn"))
			i_Mn--;
		else
			cout<<"Mn atom ID: "<<ID_rand<<" position: "<<atom[ID_rand].x[0]<<" "<<atom[ID_rand].x[1]<<" "<<atom[ID_rand].x[2]<<endl;
	}
}

void Structure::add_vacancy()
{
	if (N_s>0 || N_v>0)
	{
		cout<<"Only one saddle or one vacancy is permitted"<<endl;
		abort();
	}
	N_tot--;
	N_v++;
}

void Structure::add_saddle()
{
	if (N_s>0 || N_v>0)
	{
		cout<<"Only one saddle or one vacancy is permitted"<<endl;
		abort();
	}
	N_tot--;
	N_s++;
	for(int ix=0; ix<2; ix++)
		atom[N_tot-1].x[ix] += 0.25;
	if(bravais=="bcc")
		atom[N_tot-1].x[2] += 0.25;
	atom[N_tot-1].set_type("Fes");
}

void Structure::create_lattice(Atom *atom, string bravais)
{
	int count = 0;
	for(int ix=0; ix<N_rep; ix++)
		for(int iy=0; iy<N_rep; iy++)
			for(int iz=0; iz<N_rep; iz++)
			{
				if(bravais=="fcc")
				{
					for(int inu=0; inu<4 && count<N_tot; inu++)
					{
						atom[count].x[0] = ix+0.5*((inu==2)+(inu==3));
						atom[count].x[1] = iy+0.5*((inu==1)+(inu==3));
						atom[count].x[2] = iz+0.5*((inu==1)+(inu==2));
						if(iz%2==0)
							atom[count].flip_z();
						count++;
					}
				}
				else if(bravais=="bcc")
				{
					for(int inu=0; inu<2 && count<N_tot; inu++)
					{
						atom[count].x[0] = ix+0.5*inu;
						atom[count].x[1] = iy+0.5*inu;
						atom[count].x[2] = iz+0.5*inu;
						count++;
					}
				}
			}
	if (count!=N_tot)
	{
		cout<<"ERROR: number of atoms not correct"<<endl;
		abort();
	}
}

void Structure::reload_lattice(Atom *atom, string restart_file="config.dat")
{
	int index;
	fstream restart;
	restart.open(restart_file, ios::in);
	string type;
	double x[3], mx, my, mz, mabs, theta, phi;
	string line;
	for(int i_atom=0; i_atom<N_tot; i_atom++)
	{
		restart>>index>>type>>x[0]>>x[1]>>x[2]>>mx>>my>>mz;
		atom[index].set_type(type, true);
		for(int dim=0; dim<3; dim++)
		{
			if(x[dim]>N_rep || x[dim]<0)
			{
				cout<<"Error in position "<<x[dim]<<" of atom "<<index<<endl;
				abort();
			}
			atom[index].x[dim] = x[dim];
		}
		mabs = sqrt(mx*mx+my*my+mz*mz);
		theta = atan2(sqrt(mx*mx+my*my), mz);
		phi = atan2(my, mx);
		atom[index].set_m(mabs, theta, phi);
		if(!getline(restart, line))
		{
			cout<<"Error in config.dat at line "<<i_atom<<endl;
			abort();
		}
	}
	restart.close();
}

bool Structure::check_vacancy(int index, string bravais, kdtree *tree, Shell *shell)
{
	if(atom[index].check_saddle())
		return false;
	if(bravais=="bcc" && shell->get_shell(false, tree->get_distance(index, 7))==0)
		return false;
	if(bravais=="fcc" && shell->get_shell(false, tree->get_distance(index, 11))==0)
		return false;
	return true;
}

void Structure::initialize_tree(kdtree *tree){
	double cell[3], positions[N_tot*3];
	for(int i=0; i<3; i++)
		cell[i] = (double)N_rep;
	for(int n=0; n<N_tot; n++)
		for(int dim=0; dim<3; dim++)
			positions[n*3+dim] = atom[n].x[dim];
	tree->set_cell(cell);
	tree->set_positions(&positions[0], N_tot, 1.75);
}

float Structure::dist(Atom *xx, Atom *yy){
	float dist_tot = 0;
	for(int ix=0; ix<3; ix++)
		dist_tot += square(xx->x[ix]-yy->x[ix]-N_rep*round((xx->x[ix]-yy->x[ix])/N_rep));
	return sqrt(dist_tot);
}


//Energy::Energy(string input_file, string bravais, int N_in, float N_Mn_in, int N_s, int N_v, double lambda_in, int num_neighbors, bool debug_mode_in): acc(0), MC_count(0), kB(8.6173305e-5), N_rep(N_in), lambda(lambda_in)
Energy::Energy(Atom* atom, int n_tot, double lambda_in, bool debug_mode_in): acc(0), MC_count(0), kB(8.6173305e-5), lambda(lambda_in)
{
	cout<<"Starting initialization"<<endl;
	ausgabe.open("MC.log", ios::out);
	ausgabe<<"# Starting a Metropolis Monte Carlo simulation"<<endl;

	debug_mode = debug_mode_in;
	if(debug_mode)
		ausgabe<<"## Debug mode activated"<<endl;

	double EE = 0, E = 0, E_max, E_min, E_harm = 0;
	N_tot = n_tot;
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
		for(int i=0; config && i<N_tot; i++)
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

