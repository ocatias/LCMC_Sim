#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include <tuple>
#include <stdexcept>
#include <cmath>
#include "Eigen/Dense"

using namespace std;

class Statistics {
	private:
	public:
		static Eigen::Matrix3d orderTensor(box *boxes, int N, int relevantBaseIndex);
		static double orderParameter(box *boxes, int N, int relevantBaseIndex);
};

#endif
