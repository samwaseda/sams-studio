#ifndef TREE_H
#define TREE_H

#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

class Node
{
    private:
       // kappa : jump frequency (i.e. [vibration frequency]*exp(-dE/kBT))
        vector<double> kappa;
        double kappa_tot;
        int index, jump_ID;
        bool selected;
        Node *minor, *major, *parent;
        double kappa_sum();
        Node * create_node(vector<double>&, int index_in=-1, bool propagate_kappa=false);
        Node* return_chosen_event(double);
        bool isleaf();
        void copy_node(Node *);
    public:
        Node();
        ~Node();
        void update_kappa(vector<double>&);
        void set_kappa(vector<double>&, int index_in=-1, bool propagate_kappa=false);
        void add_kappa(double);
        void remove();
        double get_kappa();
        void append_kappa(vector<double>&, int index_in=-1);
        string get_structure();
        int get_index();
        Node * choose_event(double);
        int get_jump_ID();
};

class Tree{
    private:
        bool selected;
        Node *head, *current_node;

    public:
        Tree();
        void update_kappa(vector<double>&);
        void append(vector<double>&, int);
        int get_index();
        int get_jump_id();
        void remove();
        void choose_event(double);
        double get_kappa();
        string get_structure();
};
#endif
