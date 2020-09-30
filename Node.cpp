#include "Header.h"

using namespace std;

Node::Node(){}

Node::Node(Site* S, int k,double eps)
{
  map<int,vector<int>> all(get_all(S->g_I(),S->g_J(),k));
  for(auto& it : all)
    {
      I[it.first] = it.second[0];
      J[it.first] = it.second[1];
    }
  bool down((S->g_I()+S->g_J())%2==0);
  //cout<<S->g_Xg()<<" "<<S->g_Yg()<<endl;
  SetInitialPosition(X,Y,k,down,eps,S->g_Xg(),S->g_Yg());
  dim=S->g_dim(k);
  IX=-1;
  IY=-1;
}

Node::~Node(){}//cout<<"delete "<<this<<endl;}

map<int,int> Node::g_I() const{return I;}
map<int,int> Node::g_J() const{return J;}
int Node::g_dim() const{return dim;}

int Node::g_IX() const{return IX;}
int Node::g_IY() const{return IY;}

double Node::g_X() const{return X;}
double Node::g_Y() const{return Y;}

void Node::set_X(double x){X=x;}
void Node::set_Y(double y){Y=y;}

void Node::set_IX(int ix){IX=ix;}
void Node::set_IY(int iy){IY=iy;}

void Node::ResetPosition(int type)
{
  bool down(false);  
  if(type==0){
    if((I[0]+J[0])%2==0)
      {down=true;
	SetInitialPosition(X,Y,0,down,0.001,I[0]*0.5,J[0]*0.866+0.5774);
      }
    else{SetInitialPosition(X,Y,0,down,0.001,I[0]*0.5,J[0]*0.866+0.2887);}
  }
  else if(type==6){
    if((I[6]+J[6])%2==0)
      {down=true;
	SetInitialPosition(X,Y,0,down,0.001,I[6]*0.5,J[6]*0.866+0.5774);}
    else{SetInitialPosition(X,Y,0,down,0.001,I[6]*0.5,J[6]*0.866+0.2887);}
  }
}
