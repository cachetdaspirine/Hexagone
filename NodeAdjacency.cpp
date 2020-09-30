#include "Header.h"
/*
This file contain the structure of how the nodes are arranged.
The different function need to be adapted if we want to change
the structure of the particles.
 */
//------------------------------------------------------------
// First : we give a certain node and the function returns all
// The different index this node can have
//------------------------------------------------------------
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
                                                                                                                                                                              
      ____________                      0
     /\          /\                    /\
    /  \ i,j+1,3/  \                  /  \
   /i-1 \      /    \                /    \
  /j+1,4 \    / i+1  \              / i,j  \    (i+j)%2==1
 /        \  /,j+1,2  \            /        \
/__________\/__________\         2/__________\4
\          /\   i+1    / 
 \ i-1    /  \  ,j,1  /          1____________5
  \,j,5  /    \      /            \          /
   \    /      \    /              \  i,j   /
    \  /  i,j,0 \  /                \      /    (i+j)%2==0
     \/__________\/                  \    /
                                      \  /
                                       \/
                                        3


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
__________                 ______                        _____ _____                             
___  ____/____  __________ ___  /______ ________ ______ ___  /____(_)______ _______        _____ 
__  __/   __  |/_/___  __ \__  / _  __ `/__  __ \_  __ `/_  __/__  / _  __ \__  __ \       ___(_)
_  /___   __>  <  __  /_/ /_  /  / /_/ / _  / / // /_/ / / /_  _  /  / /_/ /_  / / /       ___   
/_____/   /_/|_|  _  .___/ /_/   \__,_/  /_/ /_/ \__,_/  \__/  /_/   \____/ /_/ /_/        _(_)  
                  /_/                                                                            
a site i,j has 3 nodes => 6 labeling depending  if the site is up or down.
Eache node belong to 6 sites, the next function  gives all the indices for
each node position. i,j referrs to the site position and k is the position
of the node with respect to a given site i,j.
*/

std::map<int,std::vector<int>> get_all(int i, int j, int k){
  // If you have a Site i,j and ask for all the adress of its node k this is what
  // this function gives you.
  std::map<int,std::vector<int>> Res;
  int N(0);
  if(k==0 | k==6){
    if(k==6){N=6;}
    Res[0+N] = {i,j};
    Res[1+N] = {i+1,j};
    Res[2+N] = {i+1,j+1};
    Res[3+N] = {i,j+1};
    Res[4+N] = {i-1,j+1};
    Res[5+N] = {i-1,j};
  }
  else if(k==5 | k==11){
    if(k==11){N=6;}
    Res[0+N] = {i+1,j};
    Res[1+N] = {i+2,j};
    Res[2+N] = {i+2,j+1};
    Res[3+N] = {i+1,j+1};
    Res[4+N] = {i,j+1};
    Res[5+N] = {i,j};
  }
  else if(k==4 | k==10){
    if(k==10){N=6;}
    Res[0+N] = {i+1,j-1};
    Res[1+N] = {i+2,j-1};
    Res[2+N] = {i+2,j};
    Res[3+N] = {i+1,j};
    Res[4+N] = {i,j};
    Res[5+N] = {i,j-1};
  }
  else if(k==3 | k==9){
    if(k==9){N=6;}
    Res[0+N] = {i,j-1};
    Res[1+N] = {i+1,j-1};
    Res[2+N] = {i+1,j};
    Res[3+N] = {i,j};
    Res[4+N] = {i-1,j};
    Res[5+N] = {i-1,j-1};
  }
  else if(k==1 | k==7){
    if(k==7){N=6;}
    Res[0+N] = {i-1,j,0};
    Res[1+N] = {i,j,0};
    Res[2+N] = {i,j+1,2};
    Res[3+N] = {i-1,j+1,3};
    Res[4+N] = {i-2,j+1,4};
    Res[5+N] = {i-2,j,5};
  }
  else if(k==2 | k==8){
    if(k==8){N=6;}
    Res[0+N] = {i-1,j-1};
    Res[1+N] = {i,j-1};
    Res[2+N] = {i,j};
    Res[3+N] = {i-1,j};
    Res[4+N] = {i-2,j};
    Res[5+N] = {i-2,j-1};
  }
  return Res;
};

void SetInitialPosition(double& X, double& Y, int NodeIndex, bool down,double eps, double Xg, double Yg){
  X = Xg;
  Y = Yg;

  if(down){// down
    if (NodeIndex==3) {
      Y -= 0.5774 * (1-eps);
    } else if (NodeIndex == 1) {
      X -= 0.5 * (1-eps);
      Y += 0.2887 * (1-eps);
    } else if(NodeIndex == 5){
      X += 0.5 * (1-eps);
      Y += 0.2887 * (1-eps);
    }
    if (NodeIndex==9) {
      Y -= 0.5774 * (1+eps);
    } else if (NodeIndex == 7) {
      X -= 0.5 * (1+eps);
      Y += 0.2887 * (1+eps);
    } else if(NodeIndex == 11){
      X += 0.5 * (1+eps);
      Y += 0.2887 * (1+eps);
    }
  }
  else{//up
    if (NodeIndex == 0) {
      Y += 0.5774 *(1-eps);
    } else if (NodeIndex == 4) {
      X += 0.5 * (1-eps);
      Y -= 0.2887 * (1-eps);
    } else if(NodeIndex ==2){
      X -= 0.5 * (1-eps);
      Y -= 0.2887 * (1-eps);
    }
    if (NodeIndex == 6) {
      Y += 0.5774 *(1+eps);
    } else if (NodeIndex == 10) {
      X += 0.5 * (1+eps);
      Y -= 0.2887 * (1+eps);
    } else if(NodeIndex ==8){
      X -= 0.5 * (1+eps);
      Y -= 0.2887 * (1+eps);
    }

  }
};
