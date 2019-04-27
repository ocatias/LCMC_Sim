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

tuple<double, Eigen::Vector3d> Statistics::orderPair(box *boxes, int N, Eigen::Matrix3d tensor)
{
  Eigen::EigenSolver<Eigen::Matrix3d> es(tensor);


  tuple<double, Eigen::Vector3d> S = make_tuple(double(es.eigenvalues().real().coeff(0)), es.eigenvectors().real().col(0));
  for(int i = 0; i < 3; i++)
  {
    if(es.eigenvalues().real().coeff(i) > get<0>(S))
      S = make_tuple(double(es.eigenvalues().real().coeff(i)), es.eigenvectors().real().col(i));
  }
  return make_tuple(get<0>(S), get<1>(S));
}

tuple<valarray<double>, valarray<double>> Statistics::orderParameter(box *boxes, int N)
{
  Eigen::Matrix3d Q0, Q1, Q2;
  Q0 = orderTensor(boxes, N, 0);
  Q1 = orderTensor(boxes, N, 1);
  Q2 = orderTensor(boxes, N, 2);


  tuple<double, Eigen::Vector3d> S20, S21, S22;
  S20 = orderPair(boxes, N, Q0);
  S21 = orderPair(boxes, N, Q1);
  S22 = orderPair(boxes, N, Q2);

  double B20, B21, B22;
  B20 = 1./3 * (double((get<1>(S21)).transpose()*(Q1*get<1>(S21))) + double((get<1>(S22)).transpose()*(Q2*get<1>(S22))) - double((get<1>(S22)).transpose()*(Q1*get<1>(S22))) - double((get<1>(S21)).transpose()*(Q2*get<1>(S21))));
  B21 = 1./3 * (double((get<1>(S20)).transpose()*(Q0*get<1>(S20))) + double((get<1>(S22)).transpose()*(Q2*get<1>(S22))) - double((get<1>(S22)).transpose()*(Q0*get<1>(S22))) - double((get<1>(S20)).transpose()*(Q2*get<1>(S20))));
  B22 = 1./3 * (double((get<1>(S21)).transpose()*(Q1*get<1>(S21))) + double((get<1>(S20)).transpose()*(Q0*get<1>(S20))) - double((get<1>(S20)).transpose()*(Q1*get<1>(S20))) - double((get<1>(S21)).transpose()*(Q0*get<1>(S21))));


  valarray<double> returnValue1 = {get<0>(S20), get<0>(S21), get<0>(S22)};
  valarray<double> returnValue2 = {B20, B21, B22};

  return make_tuple(returnValue1, returnValue2);
}
