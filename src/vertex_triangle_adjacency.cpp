#include "vertex_triangle_adjacency.h"

using namespace std;
void vertex_triangle_adjacency(
  const Eigen::MatrixXi & F,
  const int num_vertices,
  std::vector<std::vector<int> > & VF)
{
  VF.resize((unsigned long) num_vertices);
  for (int i=0;i<F.rows();i++){
    VF[F(i,0)].push_back(i);

    VF[F(i,1)].push_back(i);

    VF[F(i,2)].push_back(i);
  }
}

