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
    public:
        Node();
        ~Node();
        void set_kappa(vector<double>&, int index_in=-1, bool propagate_kappa=false);
        void add_kappa(double);
        void update_kappa();
        void remove();
        double get_kappa();
        void append_kappa(vector<double>&, int index_in=-1);
        double kappa_sum();
        Node * create_node(vector<double>&, int index_in=-1, bool propagate_kappa=false);
        void get_structure();
        bool isleaf();
        int get_index();
        Node * choose_event(double);
        Node* return_chosen_event(double);
        void copy_node(Node *);
        int get_jump_ID();
};

class Tree{
    private:
        bool selected;
        Node *head, *current_node;

    public:
        Tree();
        void append(vector<double>&, int);
        int get_index();
        int get_jump_id();
        void remove();
        void choose_event(double);
};
#endif
