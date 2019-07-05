#include<iostream>
#include<cmath>
#include<cstdio>
#include<ctime>
#include<vector> // standard KMC implementation

using namespace std;

double kBT = 300*8.617e-5;
double logtau = -13*log(10);

class Node
{
    private:                // dE : activation energy value (eV)
        double dE, kappa;   // kappa : jump frequency (i.e. [vibration frequency]*exp(-dE/kBT))
        int index;
        bool selected;
        Node *minor, *major, *parent;
    public:
        Node() : dE(0), kappa(0), index(-1), selected(false){
            minor=NULL;     // minor and major are branches conventionally called left and right
            major=NULL;     //
            parent=NULL;    //
        }

        ~Node(){
            if(selected)
            {
                parent->add_kappa(-kappa);
                delete major, delete minor;
            }
        }

        void set_dE(double dE_in, int index_in=-1, double kappa_in=0)
        {
            if (index_in>=0)
                index = index_in;
            dE = dE_in;
            if (kappa_in==0)
                add_kappa(exp(-dE/kBT-logtau));
            else
                kappa = kappa_in;
        }

        void add_kappa(double kappa_in)
        {
            kappa += kappa_in;
            if(parent!=NULL)
                parent->add_kappa(kappa_in);
        }

        double get_kappa()
        {
            if (!isleaf() && kappa==0)
                kappa = minor->get_kappa()+major->get_kappa();
            return kappa;
        }

        void append_dE(double dE_in, int index_in=-1)
        {
            if(dE==0 && isleaf())
            {
                set_dE(dE_in, index_in);
            }
            else
            {
                if(isleaf())
                {
                    minor = create_node(dE, index, kappa);
                    major = create_node(dE_in, index_in);
                    dE = 0;
                }
                else
                {
                    if(minor->get_kappa()<major->get_kappa())
                        minor->append_dE(dE_in, index_in);
                    else
                        major->append_dE(dE_in, index_in);
                }
            }
        }

        Node * create_node(double dE_in, int index_in=-1, double kappa_in=0)
        {
            Node *new_node = new Node;
            new_node->parent = this;
            new_node->set_dE(dE_in, index_in, kappa_in);
            return new_node;
        }

        void get_structure()
        {
            if(isleaf())
            {
                //cout<<"Leaf: "<<index<<" "<<parent->index<<" "<<dE<<" "<<kappa<<endl;
                return;
            }
            //cout<<"Node: "<<index<<" "<<minor->index<<" "<<major->index<<" "<<dE<<" "<<kappa<<endl;
            minor->get_structure();
            major->get_structure();
        }

        bool isleaf()
        {
            if (minor==NULL && major==NULL)
                return true;
            else
                return false;
        }

        int choose_event(double xi)
        {
            if(xi<minor->get_kappa())
            {
                if (minor->isleaf())
                    return copy_node(major, minor);
                else
                    return minor->choose_event(xi);
            }
            else
            {
                if (major->isleaf())
                    return copy_node(minor, major);
                else
                    return major->choose_event(xi-minor->get_kappa());
            }
        }

        int copy_node(Node *node_to_copy, Node *node_to_delete)
        {
            node_to_delete->selected = true;
            int index_to_return = node_to_delete->index;
            //cout<<"Deleting "<<node_to_copy->index<<" and "<<node_to_delete->index<<endl;
            delete node_to_delete;
            dE = node_to_copy->dE;
            kappa = node_to_copy->kappa;
            index = node_to_copy->index;
            Node *major_copy = node_to_copy->major;
            Node *minor_copy = node_to_copy->minor;
            delete node_to_copy;
            major = major_copy;
            minor = minor_copy;
            if (major!=NULL)
                major->parent = this;
            if (minor!=NULL)
                minor->parent = this;
            //if(isleaf())
                //cout<<"New leaf with index="<<index<<" dE="<<dE<<" kappa="<<kappa<<endl;
            //else
                //cout<<"New node with index="<<index<<" dE="<<dE<<" kappa="<<kappa<<endl;
            return index_to_return;
        }
};

class standard{
    private:
        vector<double> E_list, k_list;
        vector<int> index;
    public:
        void append_dE(double E_in, int index_in)
        {
            E_list.push_back(E_in);
            if(k_list.size()>0)
                k_list.push_back(exp(-E_in/kBT-logtau)+k_list.at(k_list.size()-1));
            else
                k_list.push_back(exp(-E_in/kBT-logtau));
            index.push_back(index_in);
        }
        double get_kappa()
        {
            return k_list.at(k_list.size()-1);
        }
        int choose_event(double xi)
        {
            for(int i=0; i<k_list.size(); i++)
            {
                xi -= k_list.at(i);
                if(xi<0)
                {
                    for(int j=i+1; j<k_list.size(); j++)
                        k_list.at(j) -= k_list.at(i);
                    E_list.erase(E_list.begin()+i);
                    k_list.erase(k_list.begin()+i);
                    int index_to_return = index.at(i);
                    index.erase(index.begin()+i);
                    return index_to_return;
                }
            }
        }
};

double get_energy() { return 0.75+0.2*rand()/RAND_MAX; }

int main(int arg, char ** name){
    srand(3);

    int n_events = atoi(name[1]);
    int n_loop = atoi(name[2]);
    cout<<n_events<<" "<<n_loop<<" ";
    double E[n_events+n_loop];
    Node *head = new Node;
    clock_t begin = clock();
    for(int i=0; i<n_events; i++)
    {
        E[i] = get_energy();
        head->append_dE(E[i], i);
    }
    for(int i=0; i<n_loop; i++)
    {
        //cout<<"Cycle: "<<i<<" Total kappa: "<<head->get_kappa()<<endl;
        //head->get_structure();
        double random_number = head->get_kappa()*rand()/RAND_MAX;
        int event = head->choose_event(random_number);
        //cout<<"Selected: "<<event<<" "<<E[event]<<" "<<exp(-E[event]/kBT-logtau)<<" "<<random_number<<endl;
        E[i+n_events] = get_energy();
        head->append_dE(E[i+n_events], n_events+i);
    }
    cout<<double(clock() - begin) / CLOCKS_PER_SEC<<" ";

    standard *kmc = new standard;
    begin = clock();
    for(int i=0; i<n_events; i++)
    {
        E[i] = get_energy();
        kmc->append_dE(E[i], i);
    }
    for(int i=0; i<n_loop; i++)
    {
        //cout<<"Cycle: "<<i<<" Total kappa: "<<kmc->get_kappa()<<endl;
        //kmc->get_structure();
        double random_number = kmc->get_kappa()*rand()/RAND_MAX;
        int event = kmc->choose_event(random_number);
        //cout<<"Selected: "<<event<<" "<<E[event]<<" "<<exp(-E[event]/kBT-logtau)<<" "<<random_number<<endl;
        E[i+n_events] = get_energy();
        kmc->append_dE(E[i+n_events], n_events+i);
    }
    cout<<double(clock() - begin) / CLOCKS_PER_SEC<<endl;
}

