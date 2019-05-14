#include "BoxCollision.h"

//scaledBase = BaseVectors * BoxDimensions
//e.g.: (l*e_x, w*e_y, h*e_z) with l,w,z ... half-dimensions ('radii') of the box
bool BoxCollision::checkProjection(vector3d T, vector3d L, vector3d scaledBase1[3], vector3d scaledBase2[3])
{
	double projection1 {0};
	double projection2 {0};

	for(int i = 0; i <= 2; i++)
	{
		//cout << i << ": " << scaledBase1[i] << ",  " << scaledBase2[i] << endl;
		projection1 += abs(scaledBase1[i]*L);
		projection2 += abs(scaledBase2[i]*L);
	}
	//<< to_string(scaledBase1) << to_string(scaledBase2)
	//cout << "In " << T << L << endl;
	//cout << to_string(abs(T*L)) << ", " << to_string(projection1) << ", " << to_string(projection2) << endl;
	return abs(T*L) > ((projection1 + projection2) + 0.0000000000001);
}

//Returns true if the boxes are colliding
bool BoxCollision::isColliding(box box1, box box2)
{
	//cout << "Checking collision: " << box1 << endl << box2 << endl;
	//cout << "Scaled Base1: " << get<0>(box1.getScaledBase()) << ", " << get<1>(box1.getScaledBase()) << ", " << get<2>(box1.getScaledBase()) << endl;
	//cout << "Scaled Base2: " << get<0>(box2.getScaledBase()) << ", " << get<1>(box2.getScaledBase()) << ", " << get<2>(box2.getScaledBase()) << endl;

	vector3d base1[3] = box1.base;
	vector3d base2[3] = box2.base;
	vector3d T (box1.center - box2.center);

	if(T.length() > sqrt(4*box1.halfRatio[0]*box1.halfRatio[0] + 4*box1.halfRatio[1]*box1.halfRatio[1] + 4*box1.halfRatio[2]*box1.halfRatio[2]))
		return false;

	vector3d scaledBase1[3] {get<0>(box1.getScaledBase()), get<1>(box1.getScaledBase()), get<2>(box1.getScaledBase())};
	vector3d scaledBase2[3] {get<0>(box2.getScaledBase()), get<1>(box2.getScaledBase()), get<2>(box2.getScaledBase())};
	return !(
		checkProjection(T, base1[0], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[1], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[2], scaledBase1, scaledBase2) ||
		checkProjection(T, base2[0], scaledBase1, scaledBase2) ||
		checkProjection(T, base2[1], scaledBase1, scaledBase2) ||
		checkProjection(T, base2[2], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[0]%base2[0], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[0]%base2[1], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[0]%base2[2], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[1]%base2[0], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[1]%base2[1], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[1]%base2[2], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[2]%base2[0], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[2]%base2[1], scaledBase1, scaledBase2) ||
		checkProjection(T, base1[2]%base2[2], scaledBase1, scaledBase2));
}
