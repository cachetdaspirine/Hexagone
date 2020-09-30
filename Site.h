#ifndef Site_h
#define Site_h

class Site
{
 public:
  Site(int i, int j,std::map<int,int> dimension,Site* Neigh);
  int g_I() const;
  int g_J() const;
  void set_G(double Xg, double Yg);
  double g_Xg() const;
  double g_Yg() const;
  int g_dim(int const index) const;
  std::vector<int> g_nodes_to_create() const;
  std::vector<int> g_nodes_to_output() const;
  bool g_down() const;
  void RemakeDim(std::map<int,int> dimension);
 private:
  int I,J;
  std::vector<int> Ineigh;
  std::vector<int> Jneigh;
  std::map<int,int> dim;
  double X,Y;
  
  
  
};
#endif
