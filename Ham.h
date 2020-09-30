#ifndef Ham_h
#define Ham_h


struct Ham{
  std::map<std::pair<Node*, Node*>, Spring*> springs;
  std::vector<Spring3*> springs3;
  double Eflip=0;
  void CheckSteadiness(VecDoub_I &x, double EmaxSpring, double EmaxSpring3)
  {
    Doub Nrj(0.);
    for(auto& it: springs){
      Nrj=it.second->ComputeNRJ(x);
      if(Nrj>EmaxSpring){EmaxSpring=Nrj;}
    }
    for(auto& it: springs3){
      Nrj=it->ComputeNRJ(x,Eflip);
      if(Nrj>EmaxSpring3){EmaxSpring3=Nrj;}
    }    
  }
  Doub operator() (VecDoub_I &x)
  {
    Doub Nrj(0);
    for(auto& it: springs){
      Nrj+=it.second->ComputeNRJ(x);
    }
    for(auto& it: springs3){
      Nrj+=it->ComputeNRJ(x,Eflip);
    }
    return Nrj;
  }
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    for (int i = 0; i < deriv.size(); i++) {deriv[i] = 0.;}
    for(auto& it : springs){
      it.second->ComputeDerivative(x,deriv);
    }
    for(auto& it: springs3){
      it->ComputeDerivative(x,deriv);
    }
  }
};
#endif
