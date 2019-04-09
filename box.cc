/*
#include "box.h"
#include <iostream>
#include "vector3d.h"
#include <tuple>
*/
#include "box.h"

ostream& operator<<(ostream& os, const box& b)
{
	os << b.center << ", " << vector3d(b.halfRatio[0], b.halfRatio[1], b.halfRatio[2]) << ", " << b.base[0] << ", " << b.base[1] << ", " << b.base[2] << endl;
	//os << acos(b.base[0]*vector3d(1,0,0)) * 180/M_PI << "°, " << acos(b.base[1]*vector3d(0,1,0)) * 180/M_PI << "°, " << acos(b.base[2]*vector3d(0,0,1)) * 180/M_PI << "°" << endl;
	return os;
}

void box::updateEdges()
{
	edges[0] = center + base[0] + base[1] + base[2];
	edges[1] = center + base[0] + base[1] - base[2];
	edges[2] = center + base[0] - base[1] + base[2];
	edges[3] = center + base[0] - base[1] - base[2];
	edges[4] = center - base[0] + base[1] + base[2];
	edges[5] = center - base[0] + base[1] - base[2];
	edges[6] = center - base[0] - base[1] + base[2];
	edges[7] = center - base[0] - base[1] - base[2];
}

void box::checkBase(vector3d baseVectors[3])
{
	/*
	if ((baseVectors[0] * baseVectors[1] != 0) || (baseVectors[1] * baseVectors[2] != 0) || (baseVectors[2] * baseVectors[0] != 0))
	{
		cout << "BASE IS NOT ORTHOGONAL: " << (baseVectors[0] * baseVectors[1]) << ", " << (baseVectors[1] * baseVectors[2]) << ", " << (baseVectors[2] * baseVectors[0])<<endl;
		cout << baseVectors[0] << ", " << baseVectors[1] << ", " << baseVectors[2] << endl;
	}
	*/
		//throw invalid_argument( "Base is not orthogonal" );
	if ((baseVectors[0].length() != 1) || (baseVectors[1].length() != 1) || (baseVectors[2].length() != 1))
	{
		base[0] = base[0].normalize();
		base[1] = base[1].normalize();
		base[2] = base[2].normalize();
	}
	//if (baseVectors[0]%baseVectors[1] != baseVectors[2])
	//	throw invalid_argument( "Base is not right handed" );
}

box::box(vector3d c, double r1, double r2, double r3, vector3d b1, vector3d b2, vector3d b3)
{
	//cout << "Full constructor";
	//cout << c << ", " << to_string(r1) << ", " << to_string(r2) << ", " << to_string(r3) << ", " << b1 << ", " << b2 << ", " << b3 << endl;
	center = c;
	halfRatio[0] = r1;
	halfRatio[1] = r2;
	halfRatio[2] = r3;
	base[0] = b1.normalize();
	base[1] = b2.normalize();
	base[2] = b3.normalize();
	checkBase(base);
	updateEdges();
	//cout << endl << b1.normalize() << base[1] << endl;
}

box::box()
{
	center = vector3d();
	halfRatio[0] = standardHalfRatio;
	halfRatio[1] = standardHalfRatio;
	halfRatio[2] = standardHalfRatio;
	base[0] = vector3d(1,0,0);
	base[1] = vector3d(0,1,0);
	base[2] = vector3d(0,0,1);
	checkBase(base);
	updateEdges();
}

box::box(vector3d c, vector3d b1, vector3d b2, vector3d b3)
{
	//box(c, standardHalfRatio, standardHalfRatio, standardHalfRatio, b1, b2, b3);
	center = c;
	halfRatio[0] = standardHalfRatio;
	halfRatio[1] = standardHalfRatio;
	halfRatio[2] = standardHalfRatio;
	base[0] = b1.normalize();
	base[1] = b2.normalize();
	base[2] = b3.normalize();
	checkBase(base);
	updateEdges();
}

void box::translate(vector3d translation)
{
	center = center + translation;
	updateEdges();
}


void box::printEdges()
{
	for(int i = 0; i <= 7; i++) {
		cout << edges[i];

		if (i != 7) cout << ", ";
	}
}
