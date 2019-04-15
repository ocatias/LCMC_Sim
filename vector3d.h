#ifndef VECTOR_3D
#define VECTOR_3D
#include <iostream>
#include <cmath>
#include <assert.h>

using namespace std;

class vector3d {
	private:

	public:
		double x, y, z;
		vector3d(double x , double y, double z) : x(x), y(y), z(z)   {   }
		vector3d() : x(0), y(0), z(0) {}
		//vector3d(vector3d vec) : x(vec.x), y(vec.y), z(vec.z) {}
		friend ostream& operator<<(ostream& os, const vector3d& vec);
		friend vector3d operator+(const vector3d& a, const vector3d& b);
		friend vector3d operator-(const vector3d& a, const vector3d& b);
		friend vector3d operator*(const vector3d& a, const double& b);
		friend vector3d operator*(const double& a, const vector3d& b);
		friend double operator*(const vector3d& a, const vector3d& b);
		friend vector3d operator%(const vector3d& a, const vector3d& b);
		friend bool operator==(const vector3d& lhs, const vector3d& rhs);
		friend bool operator!=(const vector3d& lhs, const vector3d& rhs);
		double operator[](int i);
		vector3d normalize(double r = 1);
		double length();
};

#endif
