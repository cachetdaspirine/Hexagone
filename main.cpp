#include "Header.h"
int main(int argc, char* argv[])
{
  int array[5*5]={0,0,0,0,0,
		  0,0,0,0,0,
		  0,1,0,0,0,
		  0,0,1,0,0,
                  0,0,0,0,0};
  System* system=new System(array,5,5,0.1,1.,1.,1.);
  cout<<system->get_Energy()<<endl;
  //system->OutputSpring("aight.txt");
  int array2[5*5]={0,0,0,0,0,
		  0,0,0,1,0,
		  0,0,1,1,0,
		  0,1,1,0,0,
		  0,0,0,0,0};
  System* sys2=new System(*system);
  system->UpdateEnergy(array2,5,5);
  sys2->UpdateEnergy(array2,5,5);
  cout<<system->get_Energy()<<endl;
  system->OutputSpring("aight2.txt");
  delete system;
  delete sys2;
  return 0;
}
