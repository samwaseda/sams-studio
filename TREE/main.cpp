#include<iostream>
#include<cmath>
#include "tree.h"

using namespace std;

double kBT = 300*8.617e-5;
double logtau = -13*log(10);
double get_energy(double value=0)
{
    if(value>0)
        return value;
    else
        return 0.75+0.2*rand()/RAND_MAX;
}

void error_exit(string str)
{
	cout<<"ERROR: "<<str<<endl;
	exit(EXIT_FAILURE);
}

class Atom{
    private:
        double E[4], kappa[4];
        int count;
    public:
        Atom() : count(0)
        {
            displaced();
        }

        int get_count(){
            return count;
        }

        void displaced(double value=0)
        {
            for(int i=0; i<4; i++)
            {
                E[i] = get_energy(value);
                kappa[i] = exp(-E[i]/kBT-logtau);
            }
            count++;
        }

        double * get_kappa()
        {
            return kappa;
        }
};

int main(int arg, char ** name){
    srand(3);

    #ifdef NDEBUG
    cout<<"NDEBUG has to be commented out to test the functionalities of tree.cpp"<<endl;
    return 0;
    #endif

    Node *head = new Node;
    Node *current_node;

    int n_atoms = 10;
    int n_loop = 1e2;
    if(arg==3 && false)
    {
        n_atoms = atoi(name[1]);
        n_loop = atoi(name[2]);
    }
    Atom atom[n_atoms];

    cout<<"Routine 1: One atom with E=0.8153 for all"<<endl;
    atom[0].displaced(0.81);
    head->append_kappa(atom[0].get_kappa(), 0);
    current_node = head->choose_event(0.1);
    cout<<"Check 1: index"<<endl;
    if(current_node->get_index()!=0)
        error_exit(to_string(current_node->get_index())+"!=0");
    cout<<"Check 2: jump_ID chosen"<<endl;
    if(current_node->get_jump_ID()!=0)
        error_exit(to_string(current_node->get_jump_ID())+"!=0");
    current_node->update_kappa();
    cout<<"Check 3: another jump_ID"<<endl;
    current_node = head->choose_event(0.99);
    if(current_node->get_jump_ID()!=3)
        error_exit(to_string(current_node->get_jump_ID())+"!=3");
    cout<<"Check 4: total kappa"<<endl;
    if(abs(head->get_kappa()-4*exp(-0.81/kBT-logtau))>1.0e-4)
        error_exit(to_string(head->get_kappa())+"!="+to_string(4*exp(-0.81/kBT-logtau)));
    cout<<"Check 5: deletion of current node"<<endl;
    delete current_node, head;

    cout<<"Routine 2: Two atoms with different activation energy values"<<endl;
    head = new Node;
    head->append_kappa(atom[0].get_kappa(), 0);
    atom[1].displaced(0.79);
    head->append_kappa(atom[1].get_kappa(), 1);
    cout<<"Check 1: index"<<endl;
    current_node = head->choose_event(0.5);
    if(current_node->get_index()!=1)
        error_exit(to_string(current_node->get_index())+"!=1");
    cout<<"Check 2: total kappa"<<endl;
    if(abs(head->get_kappa()-4*exp(-0.81/kBT-logtau)-4*exp(-0.79/kBT-logtau))>1.0e-4)
        error_exit(to_string(head->get_kappa())+"!="+to_string(4*exp(-0.81/kBT-logtau)+4*exp(-0.79/kBT-logtau)));
    cout<<"Check 3: deletion of current node"<<endl;
    delete current_node;
    cout<<"Check 4: appending of new branch"<<endl;
    head->append_kappa(atom[1].get_kappa(), 1);
    cout<<"Check 5: index"<<endl;
    current_node = head->choose_event(0.5);
    if(current_node->get_index()!=1)
        error_exit(to_string(current_node->get_index())+"!=1");
    cout<<"Check 6: total kappa"<<endl;
    if(abs(head->get_kappa()-4*exp(-0.81/kBT-logtau)-4*exp(-0.79/kBT-logtau))>1.0e-4)
        error_exit(to_string(head->get_kappa())+"!="+to_string(4*exp(-0.81/kBT-logtau)+4*exp(-0.79/kBT-logtau)));
    delete head;
    head = new Node;
    cout<<n_atoms<<" "<<n_loop<<" "<<endl;
    clock_t begin = clock();
    for(int i=0; i<n_atoms; i++)
        head->append_kappa(atom[i].get_kappa(), i);
    int deleted_index;
    double KMC_time = 0;
    for(int i=0; i<n_loop; i++)
    {
        cout<<"Cycle: "<<i<<" Total kappa: "<<head->get_kappa()<<" Total KMC time: "<<KMC_time<<endl;
        //head->get_structure();
        double random_number = rand()/(double)RAND_MAX;
        current_node = head->choose_event(random_number);
        //cout<<"Selected: "<<current_node->get_index()<<" "<<random_number<<endl;
        atom[current_node->get_index()].displaced();
        if(i%10==3)
        {
            deleted_index = current_node->get_index();
            delete current_node;
        }
        else
            current_node->update_kappa();
        if(i%10==7)
            head->append_kappa(atom[deleted_index].get_kappa(), deleted_index);
        KMC_time += -log(rand()/(double)RAND_MAX)/head->get_kappa();
    }
    cout<<double(clock() - begin) / CLOCKS_PER_SEC<<endl;
    for(int i=0; i<n_atoms; i++)
        cout<<atom[i].get_count()<<" ";
    cout<<endl;
}

