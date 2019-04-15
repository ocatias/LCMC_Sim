#include <iostream>
#include "vector3d.h"
#include <fstream>
#include <math.h>
#include <tuple>
#include <assert.h>

using namespace std;

int main()
{
  cout << "VECTOR3D TEST" << endl;
  vector3d a(1,2,3);
  vector3d b(1,0,-3);
  vector3d a2(1,2,3);
  vector3d c(-1.1,0,19);
  vector3d right(1,0,0);
  vector3d forward(0,1,0);
  vector3d up(0,0,1);
  vector3d zero(0,0,0);

  //Comparison
  cout << "Comparison operator: ";
    assert(a==a);
    assert(!(a==b));
    assert(a!=b);
    assert(b!=a);
  cout << "Check" << endl;

  //Scalar Product
  cout << "Scalar Product: ";
    assert(right * up == 0);
    assert(up * right == 0);
    assert(right * right == 1);
    assert(a * c == 55.9);
    assert(c * a == 55.9);
  cout << "Check" << endl;

  //Vector Product
  cout << "Vector Product: ";
    assert(a*0 == zero);
    assert(a*4.1 == vector3d(4.100000, 8.200000, 12.300000));
    assert(-4.1*a == vector3d(-4.1,-8.2,-12.3));
    assert(b*-0.1 == vector3d(-0.1, 0, 0.3));
    assert(b*-0.1 == -0.1*b);
  cout << "Check" << endl;

  //Cross Product
  cout << "Vector Product: ";
    assert(right%right == zero);
    assert(right%forward == up);
    assert(forward%right == up*-1);
    assert(a%c == vector3d(38, -22.3, 2.2));
  cout << "Check" << endl;

  //Addition
  cout << "Addition: ";
    assert(a+a == 2*a);
    assert(b+a == vector3d(2,2,0));
    assert(b+a == a+b);
    assert(a-a == zero);
    assert(b-a == vector3d(0,-2,-6));
  cout << "Check" << endl;

  //Length
  cout << "Length: ";
    assert(right.length() == 1);
    assert((up*-1).length() == 1);
    assert(a2.length() == sqrt(14));
  cout << "Check" << endl;

  //Normalization
  cout << "Normalization: ";
    assert(forward.normalize() == forward);
    assert(a2.normalize() == a2*(1/sqrt(14)));
    assert(a2.normalize(3.7) == a2*(3.7/sqrt(14)));
  cout << "Check" << endl;

  //Misc
  cout << "Misc: ";
    assert(vector3d() == zero);
  cout << "Check" << endl;

  cout << "Subscript [] operator: ";
    assert(a2[0] == 1);
    assert(a2[1] == 2);
    assert(a2[2] == 3);
    assert(c[0] == -1.1);
  cout << "Check" << endl;

  cout << "TEST COMPLETED SUCCESSFULLY" << endl;
}
