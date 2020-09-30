#include "Header.h"

using namespace std;
//This function return a list of pair of index that should make a spring
//becarefull each spring has to be made only once per site!
vector< pair<int,int> > GetSpringAdjacency(int i, int j){
  vector<pair<int,int>> Res;
  if((i+j)%2==0){
    Res.push_back({1,3});
    Res.push_back({1,5});
    Res.push_back({5,3});

    Res.push_back({7,9});
    Res.push_back({7,11});
    Res.push_back({11,9});
    /*
    Res.push_back({3,9});
    Res.push_back({1,7});
    Res.push_back({5,11});
    */
    Res.push_back({1,9});
    Res.push_back({1,11});

    Res.push_back({5,7});
    Res.push_back({5,9});

    Res.push_back({3,7});
    Res.push_back({3,11});
  }
  else{
    Res.push_back({0,2});
    Res.push_back({0,4});
    Res.push_back({2,4});

    Res.push_back({6,8});
    Res.push_back({6,10});
    Res.push_back({8,10});
    /*
    Res.push_back({0,6});
    Res.push_back({2,8});
    Res.push_back({4,10});
    */
    Res.push_back({0,8});
    Res.push_back({0,10});

    Res.push_back({4,6});
    Res.push_back({4,8});

    Res.push_back({2,6});
    Res.push_back({2,10});
  }
  return Res;
}
vector<vector<int>> GetSpring3Adjacency(int i, int j){
  vector<vector<int>> Res;
  if((i+j)%2==0){
    Res.push_back({1,3,5});
    Res.push_back({7,9,11});
  }
  else{
    Res.push_back({0,2,4});
    Res.push_back({6,8,10});
  }
  return Res;
}
double getK(int index1, int index2,System* system){
  if(index1==1 | index1==3 | index1==5 | index1==0 | index1==2 | index1==4){
    if(index2==1 | index2==0 | index2==2 | index2==3 | index2==4 |index2==5){
      return system->K1;
    }
  }
  else if(index1==6 | index1==7 | index1==8 | index1==9 | index1==10 | index1==11){
    if(index2==6 | index2==7 | index2==8 | index2==9 | index2==10 |index2==11){
      return system->K1;
    }
  }
  return system->K2; // in the C function that create system we'll find a way to adjust K
}
double getKvol(int index1, int index2, int index3,System* system){return system->Kvol;}
double getA0(int index1, int index2, int index3,System* system){
  if(index1==1 & index2==3 & index3==5){
    return sqrt(3)/4.*pow(1-system->eps,2);
  }
  else if(index1==0 & index2==2 & index3==4){
    return sqrt(3)/4*pow(1-system->eps,2);
  }
  else if(index1==7 & index2==9 & index3==11){
    return sqrt(3)/4*pow(1+system->eps,2);
  }
  else if(index1==6 & index2==8 & index3==10){
    return sqrt(3)/4*pow(1+system->eps,2);
  }
  else{
    cout<<"Error in the building of Spring3 cannot give a Rest area"<<endl;
    cout<<index1<<" "<<index2<<" "<<index3<<endl;
    exit(0);
  }
}
double getL0(int index1, int index2,System* system){
  if(index1==1 | index1==3 | index1==5 | index1==0 | index1==2 | index1==4){
    if(index2==1 | index2==0 | index2==2 | index2==3 | index2==4 |index2==5){
      return 1-system->eps;
    }
  }
  else if(index1==6 | index1==7 | index1==8 | index1==9 | index1==10 | index1==11){
    if(index2==6 | index2==7 | index2==8 | index2==9 | index2==10 |index2==11){
      return 1+system->eps;
    }
  }
  return sqrt(pow(system->eps * 2 / sqrt(3)+ sqrt(3)/2*(1-system->eps),2)+pow((1-system->eps)/2,2));
}
