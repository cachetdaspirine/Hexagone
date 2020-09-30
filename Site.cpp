#include "Header.h"

using namespace std;

Site::Site(int i, int j,map<int,int> dimension,Site* Neigh)
{
  dim=dimension;
  I=i;
  J=j;
  Ineigh=ISiteAdjacency(I,J);
  Jneigh=JSiteAdjacency(I,J);
  if(Neigh){Compute_G(this,Neigh);}
  else
    {
      X=I*0.5;
      if((I+J)%2==0){
	Y=J*0.866+0.5774;}
      else{
	Y=J*0.866+0.2887;}	
    }
  // if(Neigh){
  //   cout<<"I,J,X,Y : "<<I<<" "<<J<<" "<<X<<" "<<Y<<" "<<Neigh->g_Xg()<<" "<<Neigh->g_Yg()<<endl;}
  // else{    cout<<"I,J,X,Y : "<<I<<" "<<J<<" "<<X<<" "<<Y<<endl;}
}
void Site::RemakeDim(map<int,int> dimension){dim=dimension;}
vector<int> Site::g_nodes_to_create() const
{
  return g_nodes_from_site(I,J);
}
vector<int> Site::g_nodes_to_output() const
{
  return g_nodes_from_site_output(I,J);
}

int Site::g_I() const{return I;}
int Site::g_J() const{return J;}
void Site::set_G(double Xg,double Yg){X=Xg;Y=Yg;}
double Site::g_Xg() const {return X;}
double Site::g_Yg() const {return Y;}
int Site::g_dim(int const index)const {return dim.at(index);}
bool Site::g_down() const{if((I+J)%2==0){return true;}else{return false;}}
