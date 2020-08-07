#include<iostream>
#include<cmath>
#include<cstdio>
#include<string>
#include "Tree.h"
// #define NDEBUG

using namespace std;

Node :: Node() : kappa_tot(0), index(-1), selected(false){
    minor=NULL;     // minor and major are branches conventionally called left and right
    major=NULL;     //
    parent=NULL;    //
}

Node :: ~Node(){
    if(selected)
    {
        add_kappa(-kappa_tot); // remove kappa value from all parent nodes
        if(parent!=NULL)
        {
            if(parent->major->selected)             // not-to-be-deleted leaf/node is upbranched
                parent->copy_node(parent->minor);   // to the parent node
            else if (parent->minor->selected)       //
                parent->copy_node(parent->major);   //
        }
        else
            delete parent;
        delete major;
        delete minor;
    }
}

void Node :: remove()
{
    if(parent!=NULL)
        delete this;
    else
    {
        if (major!=NULL || minor!=NULL)
            throw invalid_argument("Leaves are not NULL");
        kappa_tot = 0;
    }
}

void Node :: set_kappa(vector<double> &kappa_in, int index_in, bool propagate_kappa)
{
    if (index_in>=0)
        index = index_in;
    kappa = kappa_in;
    if (propagate_kappa)
        add_kappa(kappa_sum());
    else
        kappa_tot = kappa_sum();
}

double Node :: kappa_sum()
{
    if(!isleaf())
        throw invalid_argument("Not leaf");
    double sum=0;
    for(int i=0; i<int(kappa.size()); i++)
        sum += kappa.at(i);
    if(sum<=0)
        throw invalid_argument("Invalid kappa");
    return sum;
}

int Node :: get_index()
{
    if(index<0)
        throw invalid_argument("index not set");
    return index;
}

void Node :: update_kappa(vector<double> &kappa_in)
{
    kappa = kappa_in;
    if (!isleaf())
        throw invalid_argument("not leaf");
    add_kappa(-kappa_tot);
    kappa_tot = 0;
    add_kappa(kappa_sum());
    selected = false;
    if(kappa_tot!=kappa_sum())
        throw invalid_argument("kappa_tot!=kappa_sum()");
}

void Node :: add_kappa(double kappa_in)
{
    kappa_tot += kappa_in;
    if(parent!=NULL)
        parent->add_kappa(kappa_in);
}

double Node :: get_kappa()
{
    if (!isleaf() && kappa_tot==0)
        kappa_tot = minor->get_kappa()+major->get_kappa();
    return kappa_tot;
}

void Node :: append_kappa(vector<double> &kappa_in, int index_in)
{
    if(kappa_tot==0 && isleaf())
    {                                   // This should be only relevant
        set_kappa(kappa_in, index_in);  // for the very first assignment of kappa
    }                                   // (maybe with assertion?)
    else
    {
        if(isleaf())
        {
            minor = create_node(kappa, index);
            major = create_node(kappa_in, index_in, true);
        }
        else // if node, continue going farther towards a leaf
        {
            if(minor->get_kappa()<major->get_kappa())
                minor->append_kappa(kappa_in, index_in);
            else
                major->append_kappa(kappa_in, index_in);
        }
    }
}

Node * Node :: create_node(vector<double> &kappa_in, int index_in, bool propagate_kappa)
{
    Node *new_node = new Node;
    new_node->parent = this;
    new_node->set_kappa(kappa_in, index_in, propagate_kappa);
    return new_node;
}

string Node :: get_structure() // only relevant for debugging
{
    if(isleaf())
        return "Leaf: "+to_string(index)+" "+to_string(parent->index)+" "+to_string(kappa_tot);
    string to_return="Node: "+to_string(index)+" "+to_string(minor->index)+" "+to_string(major->index)+" "+to_string(kappa_tot)+"\n";
    return minor->get_structure()+major->get_structure()+to_return;
}

bool Node :: isleaf()
{
    if (minor==NULL && major==NULL)
        return true;
    else // only one of them can be NULL
        return false;
}

Node* Node :: return_chosen_event(double xi)
{
    for(jump_ID=0; !selected; jump_ID++)
    {
        xi -= kappa.at(jump_ID);
        if(xi<0)
            selected = true;
    }
    if(jump_ID>int(kappa.size()))
        throw invalid_argument("Jump ID surpassed the maximum number of jump possibilities");
    return this;
}

int Node :: get_jump_ID()
{
    return jump_ID-1;
}

Node* Node :: choose_event(double xi)
{
    if(parent==NULL)
        xi *= get_kappa();
    if(xi>=get_kappa())
        throw invalid_argument("Number too large");
    if(isleaf())
        return return_chosen_event(xi);
    if(abs(kappa_tot-minor->get_kappa()-major->get_kappa())>1.0e-4)
        throw invalid_argument("Kappa value is not consistent");
    //if(minor->selected || major->selected)
    //    throw invalid_argument("Leaves probably not properly updated");
    if(xi<minor->get_kappa())
    {
        if (minor->isleaf())
            return minor->return_chosen_event(xi);
        else // look farther until leaf is reached
            return minor->choose_event(xi);
    }
    else
    {
        if (major->isleaf())
            return major->return_chosen_event(xi-minor->get_kappa());
        else
            return major->choose_event(xi-minor->get_kappa());
    }
}

void Node :: copy_node(Node *node_to_copy) // upbranching when one leaf is deleted
{
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
}

Tree :: Tree() : selected(false){
    head = new Node;
    current_node = new Node;
}

void Tree :: append(vector<double> &kappa_in, int index_in){
    head->append_kappa(kappa_in, index_in);
}

void Tree :: choose_event(double xi){
    current_node = head->choose_event(xi);
    selected = true;
}

int Tree :: get_index(){
    if (selected)
        return current_node->get_index();
    else
        throw invalid_argument("No jump selected");
}

int Tree :: get_jump_id(){
    if (selected)
        return current_node->get_jump_ID();
    else
        throw invalid_argument("No jump selected");
}

double Tree :: get_kappa(){
    return head->get_kappa();
}

void Tree :: remove(){
    if (selected)
        current_node->remove();
    else
        throw invalid_argument("No jump selected");
}

void Tree :: update_kappa(vector<double> &kappa_in){
    if (!selected)
        throw invalid_argument("No jump selected");
    current_node->update_kappa(kappa_in);
    selected = false;
}

string Tree :: get_structure(){
    return head->get_structure();
}
