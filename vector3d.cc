#include "vector3d.h"

vector3d vector3d::normalize(double r)
{
	return vector3d(x / length() * r, y / length() * r, z / length() * r);
}

double vector3d::length()
{
	return sqrt(x*x + y*y + z*z);
}

//CROSS PRODUCT
vector3d operator%(const vector3d& a, const vector3d& b)
{
	vector3d vec (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
	return vec;
}

//Scalar Product
double operator*(const vector3d& a, const vector3d& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

ostream& operator<<(ostream& os, const vector3d& vec)
{
    os << "(" << to_string(vec.x) << ", " << to_string(vec.y) << ", " << to_string(vec.z) << ")";
    return os;
}

vector3d operator*(const double& a, const vector3d& b)
{
	return b*a;
}

vector3d operator*(const vector3d& a, const double& b)
{
	vector3d vec (double(a.x * b), double(a.y * b), double(a.z * b));
	return vec;
}

vector3d operator-(const vector3d& a, const vector3d& b)
{
	return a + (b * (-1));
}

vector3d operator+(const vector3d& a, const vector3d& b)
{
	vector3d c(a.x + b.x, a.y + b.y, a.z + b.z);
	return c;
}

double vector3d::operator[](int i)
{
	assert(!((i > 2) || (i < 0)));
	if(i == 0)
		return x;
	else if(i == 1)
		return y;
	else if(i == 2)
		return z;
}


bool operator==(const vector3d& lhs, const vector3d& rhs)
{
	double epsilon = 0.00001;

	return (abs(lhs.x - rhs.x) <= epsilon) && (abs(lhs.y - rhs.y) <= epsilon) && (abs(lhs.z - rhs.z) <= epsilon);
}

bool operator!=(const vector3d& lhs, const vector3d& rhs)
{
	return !(lhs==rhs);
}
