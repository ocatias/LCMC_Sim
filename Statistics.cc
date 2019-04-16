#include "Statistics.h"

Eigen::Matrix3d Statistics::orderTensor(box *boxes, int N, int relevantBaseIndex)
{
  Eigen::Matrix3d tensor;
  tensor << 0,0,0,0,0,0,0,0,0;
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      for(int n = 0; n < N; n++)
      {
        tensor(i,j) +=  3./2*boxes[n].base[relevantBaseIndex][i]*boxes[n].base[relevantBaseIndex][j];
        if(i == j)
          tensor(i,j) -= 1./2;
      }
      tensor(i,j) /= N;
    }
  }

  return tensor;
}

double Statistics::orderParameter(box *boxes, int N, int relevantBaseIndex)
{
  Eigen::Matrix3d tensor = orderTensor(boxes, N, relevantBaseIndex);
  Eigen::EigenSolver<Eigen::Matrix3d> es(tensor);

  double S = double(es.eigenvalues().real().coeff(0));
  for(int i = 0; i < 3; i++)
  {
    if(es.eigenvalues().real().coeff(i) > S)
      S = double(es.eigenvalues().real().coeff(i));
  }
  return S;
}
