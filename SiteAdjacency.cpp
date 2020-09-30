#include "Header.h"
/*
________                            _____                  
___  __ \______________ ____      _____(_)_______ _______ _
__  / / /__  ___/_  __ `/__ | /| / /__  / __  __ \__  __ `/
_  /_/ / _  /    / /_/ / __ |/ |/ / _  /  _  / / /_  /_/ / 
/_____/  /_/     \__,_/  ____/|__/  /_/   /_/ /_/ _\__, /  
                                                  /____/

       ________       _____ ______               
______ ___  __/       __  /____  /_ _____        
_  __ \__  /_         _  __/__  __ \_  _ \       
/ /_/ /_  __/         / /_  _  / / //  __/       
\____/ /_/            \__/  /_/ /_/ \___/        
                                             

        _____                        _____                       
__________  /_____________  ___________  /_____  _______________ 
__  ___/_  __/__  ___/_  / / /_  ___/_  __/_  / / /__  ___/_  _ \
_(__  ) / /_  _  /    / /_/ / / /__  / /_  / /_/ / _  /    /  __/
/____/  \__/  /_/     \__,_/  \___/  \__/  \__,_/  /_/     \___/ 
                                                                                                                                                                              
           /\    ____________________
          /  \   \   0    /\   1    /
         / 2  \   \ i-1  /  \ i+1  /
        /  i   \   \,j  /    \,j  /
       / ,j+1   \   \  / i,j  \  /
      /__________\   \/________\/
     /\          /\   \   2    /
    /  \   i,j  /  \   \  i   /
   /i-1 \      / 1  \   \,j-1/
  / ,j   \    / i+1  \   \  /
 /   0    \  /   ,j   \   \/
/__________\/__________\


                  6            7_________________________11
                 /\             \                        /
                /  \             \    1___________ 5    /
               /    \             \    \          /    /
              /      \             \    \  i,j   /    /
             /        \             \    \      /    /
            /     0    \             \    \    /    /
           /     /\     \             \    \  /    /
          /     /  \     \             \    \/    /
         /     /    \     \             \    3   /
        /     /      \     \             \      /
       /     /   i,j  \     \             \    /
      /   2 /__________\4    \             \  /
    8/________________________\10           \/
                                             9

*/

std::vector<int> ISiteAdjacency(int i,int j){
  bool up;
  if((i+j)%2==0){up=false;}
  else{up=true;}
  std::vector<int> Res;
  Res.resize(3);
  if(up){
    Res[0] = i-1;
    Res[1] = i+1;
    Res[2] = i;
  }
  else{
    Res[0] = i-1;
    Res[1] = i+1;
    Res[2] = i;
  }
  return Res;
}
std::vector<int> JSiteAdjacency(int i,int j){
  bool up;
  if((i+j)%2==0){up=false;}
  else{up=true;}
  std::vector<int> Res;
  Res.resize(3);
  if(up){
    Res[0] = j;
    Res[1] = j;
    Res[2] = j-1;
  }
  else{
    Res[0] = j;
    Res[1] = j;
    Res[2] = j+1;
  }
  return Res;
}

std::vector<int> g_nodes_from_site(const int i, const int j)
{
  std::vector<int> Res;
  Res.resize(6);
  if((i+j)%2==0)
    {
      Res[0] = 1;
      Res[1] = 3;
      Res[2] = 5;
      Res[3] = 7;
      Res[4] = 9;
      Res[5] = 11;
    }
  else
    {
      Res[0] = 0;
      Res[1] = 2;
      Res[2] = 4;
      Res[3] = 6;
      Res[4] = 8;
      Res[5] = 10;
    }
  return Res;
}

std::vector<int> g_nodes_from_site_output(const int i, const int j)
{
  std::vector<int> Res;
  Res.resize(3);
  if((i+j)%2==0)
    {
      Res[0] = 7;
      Res[1] = 9;
      Res[2] = 11;
    }
  else
    {
      Res[0] = 6;
      Res[1] = 8;
      Res[2] = 10;
    }
  return Res;
}
int Nneigh(int i,int j,const vector<int>& State,int Lx)
{
  vector<int> IN(ISiteAdjacency(i, j));
  vector<int> JN(JSiteAdjacency(i, j));
  int Nneigh(0);
  for(int n=0;n<IN.size();n++){
    if(State[IN[n]+JN[n]*Lx]==1){Nneigh++;}
  }
  return Nneigh;
}
set<int> GetNeighNodeSharedWithTwo(int i, int j, const vector<int>& State, int Lx)
{
  vector<int> IN(ISiteAdjacency(i, j));
  vector<int> JN(JSiteAdjacency(i, j));
  set<int> Res;
  vector<int> NeighToNode(3);
  if((i+j)%2==0){
    NeighToNode={3,1,5};
  }
  else{
    NeighToNode={0,2,4};
  }
  for( int n=0;n<IN.size();n++){
    if(State[IN[n]+JN[n]*Lx]==1){
      if(State[IN[(n+1)%3]+JN[(n+1)%3]*Lx]==0){
	Res.insert(Res.begin(),NeighNode(n,NeighToNode[n+(n+1)%3-1]));
      }
      if(State[IN[(n+2)%3]+JN[(n+2)%3]*Lx]==0){
	Res.insert(Res.begin(),NeighNode(n,NeighToNode[n+(n+2)%3-1]));
      }
    }
  }
  return Res;
}
int MagicDim(int i, int j)
{
  int ij(i+3*j);
  if( ij==3 | ij==1){return 7;}
  if(ij==9 | ij==13){return 11;}
  if( ij==16 | ij==0){return 5;}
  if(ij==10 | ij==6){return 15;}
  if(ij==17 | ij==14){return 9;}
  if(ij==5 | ij==8 ){return 13;}
  cout<<"Magic dimension cannot be found i,j="<<i<<" "<<j<<endl;
  exit(0);
}
int NeighNode(int i, int j)
{
  //cout<<"NeighNode"<<endl;
  //cout<<i<<" "<<j<<endl;
  int ij(i+3*j);
  if(ij==1|ij==8){return 1;}
  if(ij==5|ij==10){return 2;}
  if(ij==13|ij==6){return 3;}
  if(ij==17|ij==9){return 4;}
  if(ij==14|ij==0){return 5;}
  if(ij==16|ij==3){return 0;}
  cout<<"cannot find the neighbors correspondance of this node IJ="<<i<<" "<<j<<endl;
  exit(0);
}
void Compute_G(Site* This,Site* Neigh)
{
  if(Neigh->g_J()==This->g_J()){
    if(Neigh->g_I()==This->g_I()-1){
      if(This->g_down()){This->set_G(Neigh->g_Xg()+0.5,Neigh->g_Yg()+0.2887);return;}
      else{This->set_G(Neigh->g_Xg()+0.5,Neigh->g_Yg()-0.2887);return;}
    }
    else if(Neigh->g_I()==This->g_I()+1){
      if(This->g_down()){This->set_G(Neigh->g_Xg()-0.5,Neigh->g_Yg()+0.2887);return;}
      else{This->set_G(Neigh->g_Xg()-0.5,Neigh->g_Yg()-0.2887);return;}      
    }
  }
  else if(Neigh->g_J()==This->g_J()-1){
    This->set_G(Neigh->g_Xg(), Neigh->g_Yg()+0.5774);return;
  }
  else if(Neigh->g_J()==This->g_J()+1){
    This->set_G(Neigh->g_Xg(), Neigh->g_Yg()- 0.5774);return;
  }
  cout<<"the site given is not a neighbor : "<<endl;
  cout<<"neighbor : "<<Neigh->g_I()<<" "<<Neigh->g_J()<<endl;
  cout<<"this pointer : "<<This->g_I()<<" "<<This->g_J()<<endl;
  exit(0);
}
set<int> GetNodeSharedWithTwo(int i, int j, const vector<int>& State, int Lx,map<int,int>& Dim)
{
  vector<int> IN(ISiteAdjacency(i, j));
  vector<int> JN(JSiteAdjacency(i, j));
  set<int> Res;
  vector<int> NeighToNode(3);
  if((i+j)%2==0){
    NeighToNode={3,1,5};
  }
  else{
    NeighToNode={0,2,4};
  }
  for( int n=0;n<IN.size();n++){
    if(State[IN[n]+JN[n]*Lx]==1){
      if(State[IN[(n+1)%3]+JN[(n+1)%3]*Lx]==0){
	Res.insert(Res.begin(),NeighToNode[n+(n+1)%3-1]);
	Dim[NeighToNode[n+(n+1)%3-1]]=MagicDim(n,NeighToNode[n+(n+1)%3-1]);
      }
      if(State[IN[(n+2)%3]+JN[(n+2)%3]*Lx]==0){
	Res.insert(Res.begin(),NeighToNode[n+(n+2)%3-1]);
	Dim[NeighToNode[n+(n+2)%3-1]]=MagicDim(n,NeighToNode[n+(n+2)%3-1]);
      }
    }
  }
   return Res;
}
map<int,int> set_dim(int SiteNum,const vector<int>& State,int Lx,int Ly)
{
  map<int,int> dim;
  vector<int> Ind(3);
  int i(SiteNum%Lx),j(SiteNum/Lx);
  int Neigh(Nneigh(i,j,State,Lx));
  vector<int> IN(ISiteAdjacency(i, j));
  vector<int> JN(JSiteAdjacency(i, j));
  
  if((i+j)%2==0){Ind[0] = 5;Ind[1] = 1;Ind[2] = 3;}
  else{Ind[0] = 4;Ind[1] = 2;Ind[2] = 0;}
  
  for(auto& it : Ind){dim[it]=0;dim[it+6]=0;}  
  if(Neigh==0){for(auto& it : Ind){dim[it]=-SiteNum;dim[it+6]=-SiteNum;}return dim;}
  
  if(Neigh==1)
    {
      for(int n=0;n<IN.size();n++){
	if(State[IN[n]+Lx*JN[n]]==1)
	  {dim[Ind[n]]=n+1;dim[Ind[n]+6]=n+1;}//dim 1,2,3
      }
    }
  map<int,int> Dim;
  set<int> sharedIJ(GetNodeSharedWithTwo(i,j,State,Lx,Dim));
  set<int> sharedNeigh,finalshared;
  for(int n=0;n<IN.size();n++){
    if(State[IN[n]+JN[n]*Lx]==1){
      set<int> sharedN(GetNeighNodeSharedWithTwo(IN[n],JN[n],State,Lx));
      sharedNeigh.insert(sharedN.begin(),sharedN.end());
    }
  }
  set_intersection(sharedIJ.begin(),sharedIJ.end(),sharedNeigh.begin(),sharedNeigh.end(),inserter(finalshared,finalshared.begin()));
  for(auto& it : finalshared){
    dim[it]=Dim.at(it);
    dim[it+6]=Dim.at(it);
  }
  return dim;
}

int g_Nnodes(){return 6;}// number of node per site
