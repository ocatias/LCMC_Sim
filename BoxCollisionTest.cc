#include <iostream>
#include "vector3d.h"
#include "box.h"
//#include "BoxCollision2.h"
#include "BoxCollision.h"
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

const int N = 1000000000;

//WIDTH/LENGTH/HEIGTH
double w = 2;
double l = 2;
double h = 2;

int hits = 0;
int tries = 0;

int mistakes = 0;

int main()
{
srand(time(NULL));
  double L =  7;
  cout << M_PI;
  vector3d e0(1,0,0);
  vector3d e1(0,1,0);
  vector3d e2(0,0,1);
  box centerBox(vector3d(0,0,0), w/2, l/2, h/2, e0, e1, e2);


  double reducedVolume = 2*(w*h*l + 2*(w*h + w*l + h*l)*1./4*(w + h + l)) / (L*L*L*4./3*M_PI);
  cout << "Simulating " << N << " tries to calculate excluded volume" << endl;
  cout << "Analytical reduced volume: " << reducedVolume << endl;
  cout << std::left << std::setw(15) << "\nTotal Diff."<< std::left << std::setw(15) << " % Diff." << "Progress" << endl;
  auto start = chrono::high_resolution_clock::now();

  for(int n = 0; n < N; n++)
  {
    if(n % 100000 == 0)
      cout <<  "\r                                             \r" << std::left << std::setw(15) <<  reducedVolume-double(hits)/tries << std::setw(15) <<  to_string(abs(reducedVolume-double(hits)/tries)/reducedVolume*100) + " %" << double(n)/N*100 << "% " << mistakes << flush;


    double x,y,z;
    do {
      x = (RandomMove::randf() - 0.5)*2*L;
      y = (RandomMove::randf() - 0.5)*2*L;
      z = (RandomMove::randf() - 0.5)*2*L;

    } while(sqrt(x*x + y*y + z*z) > L);

    vector3d centerVec(x,y,z);
    box collisionBox(centerVec, w/2, l/2, h/2, e0, e1, e2);


    tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
    //cout << "a"  <<collisionBox.base[1] << e1 << endl;
    //cout << "q" << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion) << ", " << endl;
    collisionBox.base[0] = RandomMove::rotateByQuaternion(e0, quaternion).normalize();
		collisionBox.base[1] = RandomMove::rotateByQuaternion(e1, quaternion).normalize();
		collisionBox.base[2] = RandomMove::rotateByQuaternion(e2, quaternion).normalize();
    collisionBox.checkBase();


    /*
    tuple<double, double, double> eulerAngles = RandomMove::randomEulerAngles(1);
    collisionBox.base[0] = RandomMove::rotateByEulerAngles(collisionBox.base[0].normalize(), eulerAngles).normalize();
    collisionBox.base[1] = RandomMove::rotateByEulerAngles(collisionBox.base[1].normalize(), eulerAngles).normalize();
    collisionBox.base[2] = RandomMove::rotateByEulerAngles(collisionBox.base[2].normalize(), eulerAngles).normalize();
    collisionBox.checkBase();
    */

    //vector3d base1 [3]  = {e0, e1, e2};
    //vector3d base2 [3]= {RandomMove::rotateByEulerAngles(e0, eulerAngles), RandomMove::rotateByEulerAngles(e1, eulerAngles), RandomMove::rotateByEulerAngles(e2, eulerAngles)};
    //vector3d base2 [3]= {RandomMove::rotateByQuaternion(e0, quaternion).normalize(), RandomMove::rotateByQuaternion(e1, quaternion).normalize(), RandomMove::rotateByQuaternion(e2, quaternion).normalize()};
    //cout << "b"  <<collisionBox.base[1] << base2[1] << endl;
    //cout << "q" << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion) << ", " << endl;
    //cout << RandomMove::rotateByQuaternion(e1, quaternion).normalize()<< endl;
    bool algo1Result {BoxCollision::isColliding(centerBox, collisionBox)};

    /*
    bool algo2Result {BoxCollision2::isColliding(base1, base2,vector3d(0,0,0), centerVec, vector3d(w/2, l/2, h/2), vector3d(w/2, l/2, h/2))};
    if(algo1Result != algo2Result)
    {
      mistakes++;
      cout  << "\nerror" << algo1Result <<  algo2Result << endl;
      cout << collisionBox.center << centerVec << endl;
      cout << base1[0] << centerBox.base[0] << base1[1] << centerBox.base[1] << base1[2] << centerBox.base[2] << endl;
      cout << collisionBox.base[0] << base2[0] << collisionBox.base[1] << base2[1] << collisionBox.base[2] << base2[2] << endl;
      cout << endl << endl << endl << endl;

    }
    */

    if(algo1Result)
      hits++;

    tries++;
  }
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(stop - start);

  cout << "\r                                                                         \x1b[A\r                                                                        \x1b[A\r                                                                         \r";
  cout << "\nSimulated reduced excluded volume " << double(hits)/tries << " " << endl;
  cout << "Difference " <<  reducedVolume-double(hits)/tries << endl;
  cout << "Difference " <<  abs(reducedVolume-double(hits)/tries)/reducedVolume*100 << " %" << endl;

  cout << "\nTime for MC simulation " << duration.count() << " seconds" << endl;

  //cout << "TEST COMPLETED SUCCESSFULLY" << endl;
}
