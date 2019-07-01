#include<iostream>
#include<vector>
#include<cmath>
#include<numeric>

using namespace std;

int main(){
    vector<double> v (1e7);
    for (int i=0; i<v.size(); i++)
        v.at(i) = sin(i);
    double w[v.size()];
    //w[0] = v[0];
    //for(int i=1; i<v.size(); i++)
    //    w[i] = w[i-1]+v[i];
    partial_sum(v.begin(), v.end(), w);
}
