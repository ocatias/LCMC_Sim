#ifndef RANDOM_MOVE_H
#define RANDOM_MOVE_H

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <fstream>
#include <cmath>
#include <tuple>
#include <assert.h>
#include <random>

using namespace std;

class RandomMove
{
	public:
		static double randf();
		static tuple<double, double, double, double> randomQuaternion(double delta = 1);
		static vector3d getRandomStep();
		static tuple<double, double, double, double> hamiltonProduct(tuple<double, double, double, double> a, tuple<double, double, double, double> b);
		static vector3d rotateByQuaternion(vector3d inputVector, tuple<double, double, double, double> leftQuaternion);
		static vector3d rotateByEulerAngles(vector3d vec, tuple<double, double, double> eulerAngles);
		static vector3d rotateByEulerAngles(vector3d vec, double alpha, double beta, double gamma);
		static tuple<double, double, double> randomEulerAngles(double delta);
};

#endif
