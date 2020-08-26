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
class ATOM{
friend ostream&    operator<<(ostream &out, ATOM a);

public:
    int atom_id;
    string atom_name;
    string res_name;
    string chain_id;
    int res_id;
    Cood cood;
};
class Res{
public:
    vector<ATOM> Atom;
    int res_id;
    int atom_num;
    string res_name;
};
class Chain{
public:
    vector<Res> Residues;
};
ostream&    operator<<(ostream &out, ATOM a){
    out << a.atom_id<<" "<<a.atom_name<<" "
        <<a.res_name<<" "<<a.chain_id<<" "
        <<a.res_id<<" "<<a.cood.x<<" "<<
        a.cood.y<<" "<<a.cood.z;
        return out;
}


void test01(string f){
    ifstream file;
    file.open(f, ios::in);
    if (!file.is_open()){
        cout<<endl;
        return ;
    }
    string buf;
    int rId=-999;
    int count=0;
    ATOM tmp;
    Res res_tmp;
    vector<Res> Residues;

    while (getline(file,buf)){
        if(buf.substr(0,4)=="ATOM"){
            tmp.atom_id=atoi(buf.substr(6,5).c_str());
            tmp.atom_name=buf.substr(12,4);
            tmp.res_name=buf.substr(17,3);
            tmp.chain_id=buf.substr(21,1);
            tmp.res_id=atoi(buf.substr(22,4).c_str());
            tmp.cood.x= atof(buf.substr(30,8).c_str());
            tmp.cood.y= atof(buf.substr(38,8).c_str());
            tmp.cood.z= atof(buf.substr(46,8).c_str());
            if (count==0){
                rId=tmp.res_id;
                res_tmp.res_id=tmp.res_id;
                res_tmp.res_name=tmp.res_name;
            }
            count++;
            if(rId == tmp.res_id ){
                res_tmp.Atom.push_back(tmp);
            }
            else{
                Residues.push_back(res_tmp);
                res_tmp.Atom.clear();
                rId=tmp.res_id;
                res_tmp.Atom.push_back(tmp);
                res_tmp.atom_num=count;
                count=0;
            }
            
            cout<<buf<<endl;
        }
    }
    if(! res_tmp.Atom.empty()){
        res_tmp.atom_num=++count;
        Residues.push_back(res_tmp); 
    }
    /**/
    for (vector<Res>::iterator it= Residues.begin();it!=Residues.end();it++){
        cout<<"NEW RES:"<<endl;
        cout<<"Atom NUM:"<<it->atom_num<<"\tRes_ID"<<it->res_id<<"\tRname"<<it->res_name<<endl;
        for(vector<ATOM>::iterator it2= it->Atom.begin();it2!=it->Atom.end();it2++){
            cout<<*it2<<endl;
        }
    }
    /**/
    file.close();
}

int main(){
    test01("tmp");
    return 0;
}