#ifndef BOX_H
#define BOX_H

#include <iostream>
#include "vector3d.h"
#include <tuple>
#include <stdexcept>
#include <cmath>

using namespace std;

class box {
	private:
	public:
		void checkBase();

		vector3d center;
		vector3d base[3];
		vector3d edges[8];
		double halfRatio[3];
 		double standardHalfRatio {0.5};

		void updateEdges();
		box(vector3d c, double r1, double r2, double r3, vector3d b1, vector3d b2, vector3d b3);
		box(vector3d c, vector3d b1, vector3d b2, vector3d b3);
		box();
		friend ostream& operator<<(ostream& os, const box& b);
		void print() { cout << center << ": " << base[0] << ", " << base[1] << ", " << base[2]; }
		void printEdges();

		void translate(vector3d translation);
		tuple<vector3d, vector3d, vector3d> getScaledBase() {checkBase(); return make_tuple(base[0]*halfRatio[0], base[1]*halfRatio[1], base[2]*halfRatio[2]);}
};

#endif
