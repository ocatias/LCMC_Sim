#include <iostream>
#include "vector3d.h"
//#include "BoxCollision2.h"
#include <fstream>
#include <cmath>
#include <tuple>
#include <assert.h>
#include "RandomMove.h"
#include <chrono>
#include <iomanip>
#include <string>
#include <ctime>

using namespace std;

const int N = 5000;


int main()
{
  ofstream file;
  file.open("Output/EVRGauss.txt");

  double L = 1;

  for(int n = 0; n < N; n++)
  {
    /*
    double x,y,z;
    do {
      x = (RandomMove::randf() - 0.5)*2*L;
      y = (RandomMove::randf() - 0.5)*2*L;
      z = (RandomMove::randf() - 0.5)*2*L;

    } while(sqrt(x*x + y*y + z*z) > L);


    double length = sqrt(x*x + y*y + z*z);

    x /= length;
    y /= length;
    z /= length;


    file << x  << "\t" << y << "\t" << z << endl;
    */


    vector3d e0 (1,0,0);
    vector3d output;

    tuple<double, double, double, double> quaternion;
    //do {
      quaternion = RandomMove::randomQuaternion(1);

    //} while(acos(get<0>(quaternion)) > 0.01);
    output = RandomMove::rotateByQuaternion(e0, quaternion);

    //tuple<double, double, double> eulerAngles = RandomMove::randomEulerAngles(1);
    //output = RandomMove::rotateByEulerAngles(e0, eulerAngles);


    file << output.x << "\t" << output.y << "\t" << output.z << endl;
  }
}
