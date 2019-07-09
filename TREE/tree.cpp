#include<iostream>
#include<cmath>
#include<cstdio>
#include<string>
#include "tree.h"
#include<assert.h>
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
        delete major, minor;
    }
}

void Node :: remove()
{
    if(parent!=NULL)
        delete this;
    else
    {
        assert(("Leaves are not NULL", (major==NULL && minor==NULL)));
        kappa_tot = 0;
    }
}

void Node :: set_kappa(double *kappa_in, int index_in, bool propagate_kappa)
{
    if (index_in>=0)
        index = index_in;
    kappa = kappa_in;
    if (propagate_kappa)
        add_kappa(kappa_sum());
    else
    	kappa_tot = kappa_sum();
}

double Node :: kappa_sum(int n_max)
{
    assert(isleaf());
    double sum=0;
    for(int i=0; i<n_max; i++)
        sum += kappa[i];
    assert(sum>0);
    return sum;
}

int Node :: get_index()
{
    assert(index>=0);
	return index;
}

void Node :: update_kappa()
{
    assert(isleaf());
    add_kappa(-kappa_tot);
    kappa_tot = 0;
    add_kappa(kappa_sum());
    selected = false;
    assert(kappa_tot==kappa_sum());
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

void Node :: append_kappa(double *kappa_in, int index_in)
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

Node * Node :: create_node(double *kappa_in, int index_in, bool propagate_kappa)
{
    Node *new_node = new Node;
    new_node->parent = this;
    new_node->set_kappa(kappa_in, index_in, propagate_kappa);
    return new_node;
}

void Node :: get_structure() // only relevant for debugging
{
    if(isleaf())
    {
        cout<<"Leaf: "<<index<<" "<<parent->index<<" "<<kappa_tot<<endl;
        return;
    }
    cout<<"Node: "<<index<<" "<<minor->index<<" "<<major->index<<" "<<kappa_tot<<endl;
    minor->get_structure();
    major->get_structure();
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
    assert(xi<kappa_sum());
    for(jump_ID=0; !selected; jump_ID++)
    {
        xi -= kappa[jump_ID];
        if(xi<0)
            selected = true;
    }
    assert(("Jump ID surpassed the maximum number of jump possibilities", --jump_ID<N_MAX));
    return this;
}

int Node :: get_jump_ID()
{
	return jump_ID;
}

Node* Node :: choose_event(double xi)
{
    if(parent==NULL)
        xi *= get_kappa();
    assert(xi<get_kappa());
    if(isleaf())
        return return_chosen_event(xi);
    assert(("Kappa value is not consistent", abs(kappa_tot-minor->get_kappa()-major->get_kappa())<1.0e-4));
    assert(("Leaves probably not properly updated", (!minor->selected && !major->selected)));
        
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

