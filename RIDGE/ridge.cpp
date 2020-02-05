		polynomial(int n_cv_in, bool true_error_in, bool weight_in, string input_file, bool normalize_in, bool zeroth) : N_dim(1), n_set(10), chi_old(0), val_error(0){
			weight = weight_in;
			normalize = normalize_in;
			if(n_cv_in>0 && n_cv_in<=N_line)
				n_cv = n_cv_in;
			else
				n_cv = N_line;
			true_error = true_error_in;
			n_set = n_cv;
			if(true_error)
				n_set++;
			lambda.open("lambda.dat", ios::out);
			coeff.open("coeff.dat", ios::out);
			output<<N_line<<" lines and "<<N_dim<<" terms detected\n";
			initialize_sets(input_file, skip, zeroth);
			double dd_min, dd_max, dd=0, dd_sum, dd_var;
			for(int i=0; i<n_set; i++)
			{
				dd = H[i].determinant();
				if(i==0 || dd_min>dd)
					dd_min = dd;
				if(i==0 || dd_max<dd)
					dd_max = dd;
				dd_sum += dd;
				dd_var += dd*dd;
			}
			output<<"Determinants: min = "<<dd_min<<", max = "<<dd_max<<", sum = "<<dd_sum<<" += "<<sqrt(dd_var-dd_sum*dd_sum/n_set)<<"\n";
			if (normalize)
				output<<"Normalization vector: "<<normalizer.transpose()<<"\n";
		}

void Ridge::set_number_of_cv_set(int n_set_in)
{
    if(n_set_in<=0)
        throw invalid_argument("Number of cross validation set has to be larger than one");
    n_set = n_set_in;
}

void Ridge::initilize_sets(vector<double> data, bool zeroth)
{
    H = new MatrixXd[n_set];
    hy = new VectorXd[n_set];
    yy = new double [n_set];
    for(int i=0; i<n_set; i++)
        yy[i] = 0;
    fstream eingabe;
    string linie;
    double s_in[N_dim];
    double E_in, E_sum=0, EE_sum=0;
    double g_in=1.0, g_sum[n_set];
    Lambda = VectorXd::Zero(N_dim);
    normalizer  = VectorXd::Zero(N_dim);
    increment = new double[N_dim];
    for(int i=0; i<N_dim; i++)
        increment[i] = 1.0;
    for(int i=0; i<n_set; i++)
    {
        H[i] = MatrixXd::Zero(N_dim, N_dim);
        hy[i] = VectorXd::Zero(N_dim);
        g_sum[i] = 0;
    }
    w = VectorXd::Zero(N_dim);
    eingabe.open(input_file, ios::in);
    for(int line=1; line<=N_line; line++)
    {
        int ID_tmp = line%(n_set);
        if (zeroth)
            s_in[0]=1.0;
        for(int i=int(zeroth); i<N_dim; i++)
            eingabe>>s_in[i];
        eingabe>>E_in;
        if(weight)
            eingabe>>g_in;
        g_sum[ID_tmp] += g_in;
        for(int i=0; i<N_dim; i++)
        {
            for(int j=0; j<N_dim; j++)
                H[ID_tmp](i,j) += s_in[i]*s_in[j]*g_in;
            hy[ID_tmp](i) += s_in[i]*E_in*g_in;
            normalizer(i) += s_in[i]*s_in[i]*g_in;
        }
        yy[ID_tmp] += E_in*E_in*g_in;
        E_sum += E_in;
        EE_sum += E_in*E_in;
        getline(eingabe,linie);
    }
    eingabe.close();
    for(int i=0; i<n_set; i++)
    {
        H[i] /= g_sum[i];
        hy[i] /= g_sum[i];
        yy[i] /= g_sum[i];
    }
    normalizer /= normalizer.sum()/N_dim;
    output<<"Initialization finished. E="<<E_sum/N_line<<"+-"<<sqrt(EE_sum*N_line-E_sum*E_sum)/N_line<<"\n";
    return true;
}

		bool check_input_file(string input_file, bool zeroth)
		{
			bool trenner = false;
			bool skip = false;
			fstream eingabe;
			string linie;
			eingabe.open(input_file, ios::in);
			getline(eingabe,linie);
			for(int i=0; i<linie.length(); i++)
			{
				if(i==0 && linie[i]=='#')
					skip = true;
				else if(linie[i]==' ' || linie[i]=='\t' || linie[i]==',')
				{
					if(trenner)
						N_dim++;
					trenner=false;
				}
				else
					trenner=true;
			}
			if(weight)
				N_dim--;
			if(!zeroth)
				N_dim--;
			for(N_line=1; getline(eingabe,linie); N_line++);
			if(skip)
				N_line--;
			eingabe.close();
			if(N_line==0)
			{
				output<<"ERROR: Empty file\n";
				exit(EXIT_FAILURE);
			}
			return skip;
		}

		~polynomial(){
			coeff.close();
			lambda.close();
		}

		void pearson(){
			fstream eingabe;
			double s_tmp[N_dim], E_tmp, EE_sum=0, E_sum=0;
			double x[N_dim], xx[N_dim], yy[N_dim], xy[N_dim], y[N_dim];
			string line;
			for(int i=0; i<N_dim; i++)
			{
				x[i] = 0;
				xx[i] = 0;
				yy[i] = 0;
				xy[i] = 0;
				y[i] = 0;
			}
			eingabe.open("training_set.dat", ios::in);
			s_tmp[0]=1.0;
			for(int i=0; i<N_line; i++)
			{
				for(int j=1; j<N_dim; j++)
					eingabe>>s_tmp[j];
				eingabe>>E_tmp;
				getline(eingabe, line);
				E_sum += E_tmp;
				EE_sum += E_tmp*E_tmp;
				for(int j=0; j<N_dim; j++)
				{
					x[j] += s_tmp[j];
					xx[j] += s_tmp[j]*s_tmp[j];
					y[j] += E_tmp;
					yy[j] += E_tmp*E_tmp;
					xy[j] += s_tmp[j]*E_tmp;
				}
			}
			for(int i=0; i<N_dim; i++)
			{
				double pearson_r = (N_line*xy[i]-x[i]*y[i])/sqrt((N_line*xx[i]-x[i]*x[i]+(i==0))*(N_line*yy[i]-y[i]*y[i]))+0.5*(i==0);
				Lambda(i) = -10*log(abs(pearson_r));
				output<<i<<setw(13)<<pearson_r<<setw(13)<<Lambda(i)<<"\n";
			}
			eingabe.close();
			output<<"ave(E)="<<E_sum/N_line<<"+-"<<sqrt(EE_sum*N_line-E_sum*E_sum)/N_line<<"\n";
		}

		void ridge(double lambda_in)
		{
			Lambda = VectorXd::Constant(N_dim, lambda_in);
			chi_training();
		}

		void conjugate_gradient(int max_cycle){
			VectorXd r_cg = VectorXd::Zero(N_dim);
			VectorXd p_cg = VectorXd::Zero(N_dim);
			double alpha = 0, beta = 0;
			cout<<"Starting conjugate gradient"<<endl;
			for(int i=0; i<max_cycle+1; i++)
			{
				VectorXd dchi_dlambda = VectorXd::Zero(N_dim);
				MatrixXd ddchi_ddlambda = MatrixXd::Zero(N_dim, N_dim);
				dchi_dlambda = dchi();
				ddchi_ddlambda = ddchi();
				if(i==0)
				{
					//r_cg = -dchi_dlambda-ddchi_ddlambda*Lambda;
					r_cg = -ddchi_ddlambda*Lambda;
					p_cg = r_cg;
				}
				alpha = (double) (r_cg.transpose()*r_cg)/(p_cg.transpose()*ddchi_ddlambda*p_cg);
				if(abs(alpha)<1.0e-20 || !isnormal(alpha))
					i=max_cycle;
				else
				{
					Lambda += alpha*p_cg;
					beta = (double) (r_cg.transpose()*r_cg);
					r_cg += -alpha*ddchi_ddlambda*p_cg;
					beta = (double) (r_cg.transpose()*r_cg)/beta;
					p_cg = r_cg+beta*p_cg;
					Eigen::EigenSolver<MatrixXd> es(ddchi_ddlambda,false);
					output<<setw((int)log(max_cycle)/log(10)+2)<<i*(i>0);
					output<<": alpha="<<setw(13)<<alpha;
					output<<" |dchi|="<<setw(13)<<dchi_dlambda.norm();
					output<<" |chi|="<<setw(13)<<chi(true);
					output<<" min(Eig)="<<setw(13)<<es.eigenvalues().real().minCoeff();
					output<<" |Lambda|="<<setw(13)<<Lambda.norm();
					output<<" |r|="<<setw(13)<<r_cg.norm();
					output<<" |p|="<<setw(13)<<p_cg.norm();
					output<<" |pAp|="<<setw(13)<<(double)(p_cg.transpose()*ddchi_ddlambda*p_cg);
					output<<" beta="<<setw(13)<<beta;
					output<<endl;
				}
			}
		}

		MatrixXd ddchi(){
			MatrixXd ddchi_ddlambda = MatrixXd::Zero(N_dim, N_dim);
			val_error = 0;
			for(int tr_set=0; tr_set<n_cv; tr_set++)
			{
				MatrixXd HL_tmp = MatrixXd::Zero(N_dim, N_dim);
				VectorXd hy_tmp = VectorXd::Zero(N_dim);
				MatrixXd ddchi_1 = MatrixXd::Zero(N_dim, N_dim);
				MatrixXd ddchi_21 = MatrixXd::Zero(N_dim, N_dim);
				MatrixXd ddchi_22 = MatrixXd::Zero(N_dim, N_dim);
				for(int i=0; i<n_cv; i++)
					if(i!=tr_set)
					{
						HL_tmp += H[i];
						hy_tmp += hy[i];
					}
				for(int i=0; i<N_dim; i++)
					HL_tmp(i,i) += exp(Lambda(i));
				//w = HL_tmp.inverse()*hy_tmp;
				//val_error += yy[tr_set]+ (double) (-2.0*w.transpose()*hy[tr_set]+(double)(w.transpose()*H[tr_set]*w));
				MatrixXd HL_inv = HL_tmp.inverse();
				ddchi_1 = 2*HL_inv*(-hy[tr_set]+H[tr_set]*HL_inv*hy_tmp)*hy_tmp.transpose()*HL_inv;
				//ddchi_1 += HL_inv*hy_tmp*(-hy[tr_set].transpose()+hy_tmp.transpose()*HL_inv*H[tr_set])*HL_inv;
				ddchi_21 = HL_inv*H[tr_set]*HL_inv;
				ddchi_22 = HL_inv*hy_tmp*hy_tmp.transpose()*HL_inv;
				for(int j=0; j<N_dim; j++)
					for(int k=0; k<N_dim; k++)
					{
						ddchi_1(j,k) *= HL_inv(j,k);
						ddchi_21(j,k) *= ddchi_22(j,k);
					}
				ddchi_ddlambda += 2*(ddchi_1+ddchi_21);
			}
			//val_error = sqrt(val_error/N_line/n_cv*(1.0+n_cv));
			for(int i=0; i<N_dim; i++)
				for(int j=0; j<N_dim; j++)
					ddchi_ddlambda(i,j) *= exp(Lambda(i)+Lambda(j)-Lambda.mean());
			return ddchi_ddlambda;
		}

		VectorXd dchi(){
			VectorXd dchi_dlambda = VectorXd::Zero(N_dim);
			val_error = 0;
			for(int tr_set=0; tr_set<n_cv; tr_set++)
			{
				MatrixXd HL_tmp = MatrixXd::Zero(N_dim, N_dim);
				VectorXd hy_tmp = VectorXd::Zero(N_dim);
				VectorXd dchi_1 = VectorXd::Zero(N_dim);
				VectorXd dchi_2 = VectorXd::Zero(N_dim);
				for(int j=0; j<n_cv; j++)
					if(j!=tr_set)
					{
						HL_tmp += H[j];
						hy_tmp += hy[j];
					}
				for(int j=0; j<N_dim; j++)
					HL_tmp(j,j) += exp(Lambda(j));
				//w = HL_tmp.inverse()*hy_tmp;
				//val_error += yy[tr_set]+ (double) (-2.0*w.transpose()*hy[tr_set]+(double)(w.transpose()*H[tr_set]*w));
				MatrixXd HL_inv = HL_tmp.inverse();
				dchi_1 = HL_inv*hy_tmp;
				dchi_2 = HL_inv*(-hy[tr_set]+H[tr_set]*HL_inv*hy_tmp);
				for(int j=0; j<N_dim; j++)
					dchi_dlambda(j) += -2*dchi_1(j)*dchi_2(j);
			}
			//val_error = sqrt(val_error/N_line/n_cv*(1.0+n_cv));
			for(int i=0; i<N_dim; i++)
				dchi_dlambda(i) *= exp(Lambda(i)-Lambda.mean());
			return dchi_dlambda;
		}

		void eigen_values(){
			Eigen::EigenSolver<MatrixXd> es(ddchi());
			cout << "The eigenvalues of chi are:" << endl << es.eigenvalues() << endl;
		}

		void classic(int max_step){
			double Lambda_old = Lambda.sum()/N_dim;
			Lambda = VectorXd::Constant(N_dim, Lambda_old);
			chi_training();
			double dLambda = 1.0, chi_old=val_error;
			output<<"Starting Classic calculation with error "<<chi_old<<" for Lambda="<<Lambda_old<<"\n";
			for(int i=0; i<max_step && abs(dLambda)>1.0e-8; i++)
			{
				ridge(Lambda_old+dLambda);
				if(val_error>chi_old)
				{
					dLambda *= -0.9;
					Lambda = VectorXd::Constant(N_dim, Lambda_old);
				}
				else
				{
					chi_old = val_error;
					Lambda_old += dLambda;
				}
			}
			output<<"Ending Classic calculation with CV error "<<chi_old;
			output<<" and true error "<<chi(true);
			output<<" for Lambda="<<Lambda_old<<"\n";
			output<<"Coeff: "<<w.transpose()<<"\n";
		}

		void chi_training(bool regulation=true){
			val_error = 0;
			for(int tr_set=0; tr_set<n_cv; tr_set++)
			{
				MatrixXd HL_tmp = MatrixXd::Zero(N_dim, N_dim);
				VectorXd hy_tmp = VectorXd::Zero(N_dim);
				for(int j=0; j<n_cv; j++)
					if(j!=tr_set)
					{
						HL_tmp += H[j];
						hy_tmp += hy[j];
					}
				if (regulation)
				{
					if (normalize)
					{
						for(int j=0; j<N_dim; j++)
							HL_tmp(j,j) += exp(Lambda(j))*normalizer(j);
					}
					else
					{
						for(int j=0; j<N_dim; j++)
							HL_tmp(j,j) += exp(Lambda(j));
					}
				}
				w = HL_tmp.inverse()*hy_tmp;
				val_error += yy[tr_set]+ (double) (-2.0*w.transpose()*hy[tr_set]+(double)(w.transpose()*H[tr_set]*w));
			}
			val_error = sqrt(val_error/n_cv);
		}

		double chi(bool ls){
			MatrixXd HL_tmp = MatrixXd::Zero(N_dim, N_dim);
			VectorXd hy_tmp = VectorXd::Zero(N_dim);
			int max_cycle = n_set;
			if(true_error)
				max_cycle = n_cv+int(!ls);
			for(int i=0; i<max_cycle; i++)
			{
				HL_tmp += H[i];
				hy_tmp += hy[i];
			}
			if (normalize)
			{
				for(int i=0; i<N_dim; i++)
					HL_tmp(i,i) += exp(Lambda(i))*normalizer(i);
			}
			else
			{
				for(int i=0; i<N_dim; i++)
					HL_tmp(i,i) += exp(Lambda(i));
			}
			w = HL_tmp.inverse()*hy_tmp;
			coeff<<w.transpose()<<"\n";
			lambda<<Lambda.transpose()<<"\n";
			if(true_error)
				return sqrt(yy[n_cv]-2*w.transpose()*hy[n_cv]+w.transpose()*H[n_cv]*w);
			else
				return 0;
		}

		void least_square(){
			chi_training(false);
			output<<"Ending Least Square calculation with CV error "<<val_error;
			output<<" and true error "<<chi(true);
			output<<"\n";
			output<<"Coeff: "<<w.transpose()<<"\n";
		}

		void random(int max_step){
			output<<"Starting a random calculation\n";
			if(chi_old == 0)
				chi_training();
			chi_old = val_error;
			int max_number = N_dim;
			for(int i=0; i<max_step; i++)
			{
				for(int j=0; j<max_number; j++)
				{
					int rand_ID = rand()%N_dim;
					double lambda_old = Lambda(rand_ID);
					Lambda(rand_ID) = 100.0*(1.0-2.0*rand()/(double)RAND_MAX);
					chi_training();
					if(val_error<chi_old)
					{
						chi_old = val_error;
						max_number++;
					}
					else
					{
						Lambda(rand_ID) = lambda_old;
						max_number -= ((rand()%100)==0);
					}
				}
				output<<setw((int)log(max_step)/log(10)+2)<<i*(i>0);
				output<<" |chi_t|="<<setw(13)<<chi(true);
				output<<" |chi_v|="<<setw(13)<<val_error;
				//output<<" |ddchi|="<<setw(13)<<ddchi_ddlambda.determinant();
				output<<" |Lambda|= "<<setw(13)<<Lambda.sum();
				output<<"\n";
			}
			output<<"Coeff: "<<w.transpose()<<"\n";
			output<<"Lambda: "<<Lambda.transpose()<<"\n";
		}

		double raw(int max_step){
			output<<"Starting a raw calculation\n";
			for(int i=0; i<max_step; i++)
			{
				double chi_init = val_error;
				for(int j=0; j<N_dim; j++)
				{
					int ID = rand()%N_dim;
					bool flag_stop = true;
					for(int l=0; l<i+1 && flag_stop && l<10; l++)
					{
						double lambda_old = Lambda(ID);
						if(chi_old != 0)
							Lambda(ID) += increment[ID];
						chi_training();
						if(chi_old == 0)
						{
							chi_old = val_error;
							chi_init = chi_old;
							l--;
						}
						else if(val_error>=chi_old || abs(Lambda(ID))>100)
						{
							Lambda(ID) = lambda_old;
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
				output<<setw((int)log(max_step)/log(10)+2)<<i*(i>0);
				output<<" |chi_t|="<<setw(13)<<chi(true);
				output<<" |chi_v|="<<setw(13)<<val_error;
				//output<<" |ddchi|="<<setw(13)<<ddchi_ddlambda.determinant();
				output<<" |Lambda|= "<<setw(13)<<Lambda.sum();
				output<<"\n";
				if(i>10 && abs(chi_init-val_error)<0.00001*val_error)
				{
					output<<"Calculation converged\n";
					output<<"Coeff: "<<w.transpose()<<"\n";
					output<<"Lambda: "<<Lambda.transpose()<<"\n";
					return chi_old-chi_init;
				}
			}
			return 0;
		}
