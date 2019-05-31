#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <fstream>
#include <math.h>
#include <tuple>
#include "RandomMove.h"



using namespace std;


double r = 3;
double R = 7;
int N =  1000000000;
int hits = 0;
int tries = 0;

bool isColliding(vector3d c1, vector3d c2)
{
  if((c2-c1).length() < 2*r)
    return true;
  return false;
}

int main()
{
  vector3d center(0,0,0);

  for(int n = 0; n < N; n++)
  {
    if(n % 1000000 == 0)
      cout << double(n)/N*100 << "% " << flush << "\r              \r";

    double x,y,z;
    do {
      x = (RandomMove::randf() - 0.5)*2*R;
      y = (RandomMove::randf() - 0.5)*2*R;
      z = (RandomMove::randf() - 0.5)*2*R;

    } while(sqrt(x*x + y*y + z*z) > R);

    if(isColliding(center, vector3d(x,y,z)))
      hits++;
    tries++;
  }

  cout << "Reduced excluded volume " << double(hits)/tries*4.*M_PI*R*R*R/3 << " " << endl;
  cout << "Analytical reduced volume " << 2*(4.*M_PI*r*r*r/3. + 4*M_PI*r*r*r) << endl;
  cout << "TEST COMPLETED SUCCESSFULLY" << endl;}
