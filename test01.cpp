#include<iostream>
#include<fstream>
#include<string>
#include<vector>
using namespace std ;
struct Cood{
    double x;
    double y;
    double z;
};
class Chain{
public:
    string Chain_id;
    int Res_num;
    int Atom_num;
};
class Res:public Chain{
public:
    string Res_name;
    int Res_id;
};
class Atom:public Res {
    public:
    string Atom_name;
    int Atom_id;
    Cood cood;
};

void readpdb(string filename){
    ifstream file;
    file.open(filename,ios::in);
    if(!file.is_open()){
        cout<<"文件无法打开"<<endl;
        return;
    }
    string buf;
    int rId=-999;
    int count=0;
    Atom atom_tmp ;
    Chain chain_tm;
    Res res_tmp;
    vector<Atom> A;
    string buf;
    while (getline(file,buf)){
        if(buf.substr(0,4)=="ATOM"){
            atom_tmp.Atom_id=atoi(buf.substr(6,5).c_str());
            atom_tmp.Atom_name=buf.substr(12,4);
            atom_tmp.Res_name=buf.substr(17,3);
            atom_tmp.Chain_id=buf.substr(21,1);
            atom_tmp.Res_id=atoi(buf.substr(22,4).c_str());
            atom_tmp.cood.x= atof(buf.substr(30,8).c_str());
            atom_tmp.cood.y= atof(buf.substr(38,8).c_str());
            atom_tmp.cood.z= atof(buf.substr(46,8).c_str());
        }
        A.push_back(atom_tmp);
    }
}