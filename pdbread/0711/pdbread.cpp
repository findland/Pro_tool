#include"pdbread.h"
using namespace std;

void Residue::init(){
    this->Atom_num=0;
    this->Atoms.clear();
}
void Residue::add_Atom(Atom a){
    this->Atom_num++;
    this->Atoms.push_back(a);
}
bool Residue::is_stand_residue(){
    int res_type =find_res(Res_name);
    int count =0;
    for(int i =0 ; i<res_atom[res_type].length();i+=4){
        for(vector<Atom>::iterator j = Atoms.begin();j!=Atoms.end();j++){
            if(res_atom[res_type].substr(i,4)==j->Atom_name){
                count++;
                break;
            }
        }
    }
    if (count == Atom_num && Atom_num == residue_stand_size[res_type])
        return true;
    else 
        return false ;
}

Cood Residue::Center(){
    double Cx=0;
    double Cy=0;
    double Cz=0;
    int count=0;
    for(vector<Atom>::iterator i = Atoms.begin();i!=Atoms.end() ; i++){
        Cx+=i->cor.x;
        Cy+=i->cor.y;
        Cz+=i->cor.z;
        count++;
    }
    Cood center;
    center.x=Cx/count;
    center.y=Cy/count;
    center.z=Cz/count;
    return center;
}

void Chain::init(){
    this->Res_num=0;
    this->Atom_num=0;
    this->Residues.clear();
}

void Chain::add_Res(Residue a){
    this->Res_num++;
    this->Atom_num+=a.Atom_num;
    this->Residues.push_back(a);
}

Cood Chain::Center(){
    double Cx=0;
    double Cy=0;
    double Cz=0;
    int count=0;
    for(vector<Residue>::iterator i = Residues.begin();i!=Residues.end() ; i++){
        for(vector<Atom>::iterator j = i->Atoms.begin();j!=i->Atoms.end() ;j++){
            Cx+=j->cor.x;
            Cy+=j->cor.y;
            Cz+=j->cor.z;
            count++;
        }
    }
    Cood center;
    center.x=Cx/count;
    center.y=Cy/count;
    center.x=Cz/count;
    return center;
}

Complex::Complex(string filename,string Chain_ID){

    if(Chain_ID==""){
        Chain_ID=get_chain_name(filename);
    }
    string buf;
    Atom Atom_tmp;
    Residue Residue_tmp;
    Chain Chain_tmp;
    string Res_name;
    string Chain_name;
    int Res_ID;

    ifstream file;
    file.open(filename,ios::in);
    if(!file.is_open()){
        cout<<filename<<"打开失败"<<endl;
    }

    while (getline(file,buf)){
        if(buf.substr(0,4)=="ATOM"){
            Atom_tmp.Atom_ID=atoi(buf.substr(6,5).c_str());
            Atom_tmp.Atom_name=buf.substr(12,4);
            Res_name=buf.substr(17,3);
            if(!is_stand_atom(Res_name,Atom_tmp.Atom_name)) continue;
            Chain_name=buf.substr(21,1);
            if(Chain_ID.find(Chain_name)==string::npos){
                continue;
            }// 读取特定链
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
                    if(Residue_tmp.is_stand_residue()){
                        Chain_tmp.Chain_name=Chain_name;
                        Chain_tmp.add_Res(Residue_tmp);
                    }
                }
                // 判断是否换链,如果换链，
                // 则加入当前残基，
                // 并把当前链加入到complex
                // 初始化，并改链名
                else if (Chain_tmp.Chain_name!=Chain_name){
                    if(Residue_tmp.is_stand_residue()){
                        Chain_tmp.add_Res(Residue_tmp);
                        add_Chain(Chain_tmp);
                        Chain_tmp.init();
                        Chain_tmp.Chain_name=Chain_name;
                    }
                    
                }
                else{
                    if(Residue_tmp.is_stand_residue()){
                        Chain_tmp.add_Res(Residue_tmp);
                    }
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
    if(Residue_tmp.is_stand_residue()){
        Chain_tmp.add_Res(Residue_tmp);
    }
    add_Chain(Chain_tmp);
}

Complex::Complex(const Complex &com, string Cname){
    if(Cname==""){
        Cname=get_chain_name(com);
    }
    for(int i = 0 ; i<Cname.length();i++){
        int flag=0;
        for(vector<Chain>::const_iterator j =com.Chains.begin();j!=com.Chains.end();j++){
            if (Cname.substr(i,1)== j->Chain_name )
            {
                add_Chain(*j);
                flag=1;
            }
        }
        if (flag ==0){cout<<"在原复合体中没有"<<Cname.substr(i,1)<<"链"<<endl;}
    }
}


void Complex::add_Chain(Chain a){
    this->Chain_num++;
    this->Res_num+=a.Res_num;
    this->Atom_num+=a.Atom_num;
    Chains.push_back(a);
}
Cood Complex::Center(){
    double Cx=0;
    double Cy=0;
    double Cz=0;
    int count=0;
    for(vector<Chain>::iterator i = Chains.begin();i!=Chains.end() ; i++){
        for(vector<Residue>::iterator j = i->Residues.begin();j!=i->Residues.end() ; j++){
            for(vector<Atom>::iterator k = j->Atoms.begin();k!=j->Atoms.end() ;k++){
                Cx+=k->cor.x;
                Cy+=k->cor.y;
                Cz+=k->cor.z;
                count++;
            }
        }
    }
    Cood center;
    center.x=Cx/count;
    center.y=Cy/count;
    center.x=Cz/count;
    return center;
}



double distance(const Cood &a ,const Cood &b){
    double dis = sqrt((a.x-b.x)*(a.x-b.x) +(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
    return dis;
}

// ostream& operator<<(ostream& out,Cood &a){
//     out<<a.x<<" "<<a.y<<" "<<a.z;
//     return out;
// }
ostream& operator<<(ostream& out,const Cood &a){
    out<<a.x<<" "<<a.y<<" "<<a.z;
    return out;
}


string get_chain_name(const Complex &com){
    string Cname="";
    for(vector<Chain>::const_iterator i =com.Chains.begin();i!=com.Chains.end();i++){
        Cname+=i->Chain_name;
    }
    return  Cname; 
}


string get_chain_name(string filename){
    ifstream file;
    file.open(filename,ios::in);
    string chain_name="";
    string tmp="";
    string buf;
    while(getline(file,buf)){
        if(buf.substr(0,4)=="ATOM"){
            string c_name=buf.substr(21,1);
            if(c_name!=tmp){
                tmp=c_name;
                chain_name+=tmp;
            }
        }//if(buf.substr(0,4)=="ATOM"){
    }// while(getline(file,buf)){
    return chain_name;
}

int find_res(string a){
    int i;
    for (i =0;i<20;i++){
        if (residue_name[i]==a || residue_name_short[i]==a)
            break;
    }
    return i;
}
string res_name_l2s(string a){
    int res_num=find_res(a);
    if (res_num==20){
        cout<<"该残基为非标准残基"<<endl;
        return "";
    }
    return residue_name_short[res_num];
}
string res_name_s2l(string a){
    int res_num=find_res(a);
    if (res_num==20){
        cout<<"该残基为非标准残基"<<endl;
        return "";
    }
    return residue_name[res_num];
}
bool is_stand_atom(string Resname ,string Atomname){
    int restype=find_res(Resname);
    if (restype==20) return false ;
    for (int i =0 ;i <res_atom[restype].length();i+=4){
        if(Atomname== res_atom[restype].substr(i,4)) return true;
    }
    return false;
}




void show_chain_head(const Complex &com){
    for(vector<Chain>::const_iterator i =com.Chains.begin();i!=com.Chains.end();i++){
        int count=0;
        for(vector<Residue>::const_iterator j =i->Residues.begin();j!=i->Residues.end();j++){
            for(vector<Atom>::const_iterator k =j->Atoms.begin();k!=j->Atoms.end();k++){
                cout<<k->Atom_ID<<" "<<k->Atom_name<<" "
                    <<j->Res_name<<" "<<i->Chain_name<<" "
                    <<j->Res_ID<<" "<<k->cor<<endl;
            }
        count++;
        if(count==2) break;
        }
    }
}

//计算两个残基之间的中心距离，
double Res_center_dis ( Residue res1, Residue res2){
    return distance(res1.Center(),res2.Center());
}
//计算两个残基之间的最小距离
double Res_min_dis( Residue res1 , Residue res2){
    double min=Res_center_dis(res1,res2);
    for(vector<Atom>::iterator i =res1.Atoms.begin();i!=res1.Atoms.end();i++){
        for(vector<Atom>::iterator j =res2.Atoms.begin();j!=res2.Atoms.end();j++){
            double tmp = distance(i->cor,j->cor);
            if (min >tmp ){
                min=tmp;
            }
        }
    }
    return min;
}

//　计算两个复合体中的残基之间的contact，截断距离由用户指定
int  get_contact(const Complex &com1,const Complex &com2,double dis,string outfile,string mod="center"){
    int count=0;
    ofstream data_out;
    data_out.open(outfile,ios::out);
    double (*pf)(Residue,Residue);//使用函数指针代替两种函数
    if(mod=="center"){
        pf = Res_center_dis;
    }
    else if (mod=="min"){
        pf = &Res_min_dis;
    }
    //在给函数指针赋值时，目标的函数名前
    //加不加取址符&都可以，因为函数名本身就代表了地址
    else{
        cout<<"方法输入错误！"<<endl;
        return count;
    }
    for(vector<Chain>::const_iterator C1=com1.Chains.begin();C1!=com1.Chains.end();C1++)
        for(vector<Residue>::const_iterator R1 =C1->Residues.begin();R1!=C1->Residues.end();R1++){
    for(vector<Chain>::const_iterator C2=com2.Chains.begin();C2!=com2.Chains.end();C2++)
        for(vector<Residue>::const_iterator R2 =C2->Residues.begin();R2!=C2->Residues.end();R2++){
            if(pf(*R1,*R2)<=dis){
                data_out<<C1->Chain_name<<" "<<R1->Res_ID<<" "<<R1->Res_name
                <<"\t"<<C2->Chain_name<<" "<<R2->Res_ID<<" "<<R2->Res_name 
                <<endl;
                count++;
            }
        }
    }
    data_out <<"共存在"<<count<<"个配对"<<endl;
    data_out.close();
    return count ;
}


int main(){
    Complex a("../6acg.pdb");
    show_chain_head(a);
    Complex b(a,"A");
    Complex c(a,"C");
    get_contact(b,c,5,"contact_AC.dat","min");


    // for (vector<Chain>::iterator i =a.Chains.begin();i!=a.Chains.end();i++){
    //     for (vector<Residue>::iterator j =i->Residues.begin();j!=i->Residues.end();j++){
    //         if(j->is_stand_residue())
    //             cout<<j->Res_ID<<"\t"<<j->Res_name<<":yes"<<endl;
    //         else 
    //             cout<<j->Res_ID<<"\t"<<j->Res_name<<":no"<<endl;
    //     }
    // }


    /*
    for(vector<Chain>::iterator i =a.Chains.begin();i!=a.Chains.end();i++){
        for(vector<Residue>::iterator j =i->Residues.begin();j!=i->Residues.end();j++){
            for(vector<Atom>::iterator k =j->Atoms.begin();k!=j->Atoms.end();k++){
                cout<<k->Atom_ID<<" "<<k->Atom_name<<" "
                    <<j->Res_name<<" "<<i->Chain_name<<" "
                    <<j->Res_ID<<" "<<k->cor.x<<" "<<k->cor.y<<" "<<k->cor.z<<endl;
            }
        }
    }
    */
   /*
   for(int l=0 ;l<20 ;l++){
       int flag=0;
        for(vector<Chain>::iterator i =a.Chains.begin();i!=a.Chains.end();i++){
            for(vector<Residue>::iterator j =i->Residues.begin();j!=i->Residues.end();j++){
                if (residue_name[l]==j->Res_name && flag ==0 ){
                    // cout<<j->Res_name<< "\t\"";
                    cout <<"\"";
                    for(vector<Atom>::iterator k =j->Atoms.begin();k!=j->Atoms.end();k++){
                            cout<<k->Atom_name;
                    }
                    cout <<"\","<<endl;
                    flag=1;
                }
            }// loop in residue
        } //loop in chain
    }
    */
    // for (int i=0;i<20;i++){
    //     cout << res_atom[i].length()/4<<",";
    // }
    // cout <<endl;
   
   //cout<<get_chain_name("../6acg.pdb")<<endl;
   // test01();
    return 0;
}