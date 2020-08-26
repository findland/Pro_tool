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
class Atom{
public:
    string Atom_name;
    int Atom_ID;
    Cood cor;
};
class Residue{
public:
    string Res_name="*";
    int Res_ID=-999;
    int Atom_num=0;
    vector<Atom> Atoms;

    Cood Center();
    void init(){
        this->Atom_num=0;
        this->Atoms.clear();
    }
    void add_Atom(Atom a){
        this->Atom_num++;
        this->Atoms.push_back(a);
    }
};
class Chain{
public:
    string Chain_name="*";
    int Res_num=0;
    int Atom_num=0;
    vector<Residue>  Residues;

    Cood Center();
    void init(){
        this->Res_num=0;
        this->Atom_num=0;
        this->Residues.clear();
    }
    void add_Res(Residue a){
        this->Res_num++;
        this->Atom_num+=a.Atom_num;
        this->Residues.push_back(a);
    }

};
class Complex{
public:
    int Chain_num=0;
    int Res_num=0;
    int Atom_num=0;
    vector<Chain> Chains;

    void add_Chain(Chain a){
        this->Chain_num++;
        this->Res_num+=a.Res_num;
        this->Atom_num+=a.Atom_num;
        Chains.push_back(a);
    }
    Cood Center();

};

Complex readpdb(string filename){

    string buf;
    Atom Atom_tmp;
    Residue Residue_tmp;
    Chain Chain_tmp;
    Complex pro;
    string Res_name;
    string Chain_name;
    int Res_ID;

    ifstream file;
    file.open(filename,ios::in);
    if(!file.is_open()){
        cout<<filename<<"打开失败"<<endl;
        return pro;
    }

    while (getline(file,buf)){
        if(buf.substr(0,4)=="ATOM"){
            Atom_tmp.Atom_ID=atoi(buf.substr(6,5).c_str());
            Atom_tmp.Atom_name=buf.substr(12,4);
            Res_name=buf.substr(17,3);
            Chain_name=buf.substr(21,1);
            Res_ID=atoi(buf.substr(22,4).c_str());
            Atom_tmp.cor.x= atof(buf.substr(30,8).c_str());
            Atom_tmp.cor.y= atof(buf.substr(38,8).c_str());
            Atom_tmp.cor.z= atof(buf.substr(46,8).c_str());
        
            if(Residue_tmp.Res_ID==-999){
                Residue_tmp.Res_name=Res_name;
                Residue_tmp.Res_ID=Res_ID;
                Residue_tmp.add_Atom(Atom_tmp);
            }
            else if(Residue_tmp.Res_ID!=Res_ID){
                //判断是否是初次储存Chain,初次则改名，加残基
                if(Chain_tmp.Chain_name=="*"){
                    Chain_tmp.Chain_name=Chain_name;
                    Chain_tmp.add_Res(Residue_tmp);
                }
                // 判断是否换链,如果换链，
                // 则加入当前残基，
                // 并把当前链加入到complex
                // 初始化，并改链名
                else if (Chain_tmp.Chain_name!=Chain_name){
                    Chain_tmp.add_Res(Residue_tmp);
                    pro.add_Chain(Chain_tmp);
                    Chain_tmp.init();
                    Chain_tmp.Chain_name=Chain_name;
                    
                }
                else{
                    Chain_tmp.add_Res(Residue_tmp);
                }
                //将当前的residue加入到Chain中
                //初始化residue，重新赋值
                Residue_tmp.init();
                Residue_tmp.Res_name=Res_name;
                Residue_tmp.Res_ID=Res_ID;
                Residue_tmp.add_Atom(Atom_tmp);
            }
            else{
                Residue_tmp.add_Atom(Atom_tmp);
            }
        }
    }
    //读取后最后一个残基收尾
    Chain_tmp.add_Res(Residue_tmp);
    pro.add_Chain(Chain_tmp);
    return pro;
}

int main(){
    Complex a = readpdb("6acg.pdb");
    for(vector<Chain>::iterator i =a.Chains.begin();i!=a.Chains.end();i++){
        for(vector<Residue>::iterator j =i->Residues.begin();j!=i->Residues.end();j++){
            for(vector<Atom>::iterator k =j->Atoms.begin();k!=j->Atoms.end();k++){
                cout<<k->Atom_ID<<" "<<k->Atom_name<<" "
                    <<j->Res_name<<" "<<i->Chain_name<<" "
                    <<j->Res_ID<<" "<<k->cor.x<<" "<<k->cor.y<<" "<<k->cor.z<<endl;
            }
        }
    }
    return 0;
}