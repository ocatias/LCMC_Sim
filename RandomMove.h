#ifndef RANDOM_MOVE_H
#define RANDOM_MOVE_H

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <fstream>
#include <math.h>
#include <tuple>

using namespace std;

class RandomMove
{
	public:
		static double randf();
		static tuple<double, double, double, double> randomQuaternion();
		static vector3d getRandomStep();
		static tuple<double, double, double, double> hamiltonProduct(tuple<double, double, double, double> a, tuple<double, double, double, double> b);
		static vector3d rotateByQuaternion(vector3d inputVector, tuple<double, double, double, double> leftQuaternion);


};

#endif
