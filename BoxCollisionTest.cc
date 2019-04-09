#include <iostream>
#include "vector3d.h"
#include "box.h"
#include "BoxCollision.h"
#include <fstream>
#include <math.h>
#include <tuple>
#include <assert.h>

using namespace std;

int main()
{
  vector3d zero;
  vector3d right(1.1,0,0);
  vector3d forward(0,1.1, 0);
  vector3d up(0,0,1.1);

  box cubeCenter;

  box cubeUp(up, right, forward, up);
  box cubeRight(right, 10, 10, 10, right, forward, up);
  box cubeDown(-1*up, right, forward, up);
  cout << cubeCenter << cubeRight << cubeDown;

  cout << "BOX COLLISION TEST" << endl;
    assert(BoxCollision::isColliding(cubeCenter, cubeCenter) == true);
    assert(BoxCollision::isColliding(cubeCenter, cubeUp) == false);
    assert(BoxCollision::isColliding(cubeDown, cubeUp) == false);
    assert(BoxCollision::isColliding(cubeCenter, cubeRight) == true);

  cout << "TEST COMPLETED SUCCESSFULLY" << endl;
}
