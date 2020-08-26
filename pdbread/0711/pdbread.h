#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstring>
#include<cmath>
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
    void init();
    bool is_stand_residue();
    void add_Atom(Atom a);
};


class Chain{
public:
    string Chain_name="*";
    int Res_num=0;
    int Atom_num=0;
    vector<Residue>  Residues;

    Cood Center();
    void init();
    void add_Res(Residue a);
};
class Complex{
public:
    int Chain_num=0;
    int Res_num=0;
    int Atom_num=0;
    vector<Chain> Chains;

    void add_Chain(Chain a);
    Complex(){}
    Complex(string filename,string Chain_ID="");
    Complex(const Complex &com , string Cname="");
    Cood Center();
};

double distance(const Cood &a,const Cood &b);
ostream& operator<<(ostream& out ,const Cood &a);

const string residue_name[20]={
    "GLY",    "ALA",    "VAL",    "LEU",    
    "ILE",    "PHE",    "TRP",    "TYR", 
    "ASP",    "ASN",    "GLU",    "LYS",    
    "GLN",    "MET",    "SER",    "THR",    
    "CYS",    "PRO",    "HIS",    "ARG"};
const string residue_name_short[20]={
    "G",    "A",    "V",    "L",    "I",  
    "F",    "W",    "Y",    "D",    "N",   
    "E",    "K",    "Q",    "M",    "S",    
    "T",    "C",    "P",    "H",    "R"};
int find_res(string a);
string res_name_l2s(string a);
string res_name_s2l(string a);

string res_atom[20]={
" N   CA  C   O  ",
" N   CA  C   O   CB ",
" N   CA  C   O   CB  CG1 CG2",
" N   CA  C   O   CB  CG  CD1 CD2",
" N   CA  C   O   CB  CG1 CG2 CD1",
" N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ ",
" N   CA  C   O   CB  CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
" N   CA  C   O   CB  CG  CD1 CD2 CE1 CE2 CZ  OH ",
" N   CA  C   O   CB  CG  OD1 OD2",
" N   CA  C   O   CB  CG  OD1 ND2",
" N   CA  C   O   CB  CG  CD  OE1 OE2",
" N   CA  C   O   CB  CG  CD  CE  NZ ",
" N   CA  C   O   CB  CG  CD  OE1 NE2",
" N   CA  C   O   CB  CG  SD  CE ",
" N   CA  C   O   CB  OG ",
" N   CA  C   O   CB  OG1 CG2",
" N   CA  C   O   CB  SG ",
" N   CA  C   O   CB  CG  CD ",
" N   CA  C   O   CB  CG  ND1 CD2 CE1 NE2",
" N   CA  C   O   CB  CG  CD  NE  CZ  NH1 NH2"
};
int residue_stand_size[20]={4,5,7,8,8,11,14,12,8,8,9,9,9,8,6,7,6,7,10,11};

bool is_stand_atom(string Resname,string Atomname);
void show_chain_head(const Complex &com);
double Res_center_dis ( Residue res1, Residue res2);
double Res_min_dis( Residue res1 , Residue res2);
string get_chain_name(const Complex &com);
string get_chain_name(string filename);
//计算两个复合体之间的contact，contact由自己定义
int  get_contact(const Complex &com1,const Complex &com2,double dis,string outfile,string mod);
