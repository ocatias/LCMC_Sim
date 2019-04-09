#ifndef BOX_COLLISION_H
#define BOX_COLLISION_H

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <fstream>
#include <math.h>
#include <tuple>

using namespace std;

class BoxCollision
{
	public:
		//scaledBase = BaseVectors * BoxDimensions
		//e.g.: (l*e_x, w*e_y, h*e_z) with l,w,z ... half-dimensions ('radii') of the box
		static bool checkProjection(vector3d T, vector3d L, vector3d scaledBase1[3], vector3d scaledBase2[3]);

		//Returns true if the boxes are colliding
		static bool isColliding(box box1, box box2);
};

#endif
