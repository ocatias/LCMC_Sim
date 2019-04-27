#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <tuple>
#include <stdexcept>
#include <cmath>
#include "Eigen/Dense"
#include <valarray>

using namespace std;

class Statistics {
	private:
	public:
		static Eigen::Matrix3d orderTensor(box *boxes, int N, int relevantBaseIndex);
		static tuple<double, Eigen::Vector3d> orderPair(box *boxes, int N, Eigen::Matrix3d tensor);
		static tuple<valarray<double>, valarray<double>> orderParameter(box *boxes, int N);
};

#endif
