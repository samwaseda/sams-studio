#include<iostream>
#include<cmath>
#include<cstdio>
#include<ctime>
#include<vector> // standard KMC implementation
#include "tree.h"

using namespace std;

double kBT = 300*8.617e-5;
double logtau = -13*log(10);
double get_energy() { return 0.75+0.2*rand()/RAND_MAX; }

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

        void displaced()
        {
            for(int i=0; i<4; i++)
            {
                E[i] = get_energy();
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

    int n_atoms = 10;
    int n_loop = 1e2;
    if(arg==3)
    {
        n_atoms = atoi(name[1]);
        n_loop = atoi(name[2]);
    }
    Atom atom[n_atoms];
    //for(int i=0; i<n_atoms; i++)
    //    cout<<atom[i].get_kappa()[0]<<" "<<atom[i].get_kappa()[1]<<" "<<atom[i].get_kappa()[2]<<" "<<atom[i].get_kappa()[3]<<endl;
    cout<<n_atoms<<" "<<n_loop<<" "<<endl;
    Node *head = new Node;
    Node *current_node;
    clock_t begin = clock();
    for(int i=0; i<n_atoms; i++)
        head->append_kappa(atom[i].get_kappa(), i);
    int deleted_index;
    double KMC_time = 0;
    for(int i=0; i<n_loop; i++)
    {
        cout<<"Cycle: "<<i<<" Total kappa: "<<head->get_kappa()<<" Total KMC time: "<<KMC_time<<endl;
        //head->get_structure();
        double random_number = head->get_kappa()*rand()/RAND_MAX;
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

