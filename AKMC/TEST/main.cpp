#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

int main(){
    fstream eingabe;
    eingabe.open("init.dat", ios::in);
    string line;
    int pos;
    double vv;
    vector<double> value;
    value.assign(1000000, 0);
    while(getline(eingabe, line))
    {
        pos = line.rfind("#");
        if (pos<line.length())
            line = line.substr(0, pos);
        stringstream ss(line);
        for(int i=0; ss>>vv; i++)
            //value.at(i) = vv;
            value.push_back(vv);
    }
}
