#include "RandomMove.h"

double RandomMove::randf()
{
	return ((double) rand() / (RAND_MAX));
}


tuple<double, double, double, double> RandomMove::randomQuaternion()
{
	double s = randf();
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2*M_PI*randf();
	double theta2 = 2*M_PI*randf();
	double w = cos(theta2)*sigma2;
	double x = sin(theta1)*sigma1;
	double y = cos(theta1)*sigma1;
	double z = sin(theta2)*sigma2;
	//cout << w << ", " << x << ", " << y << ", " << z << endl;
	//cout << "Quaternion size = " << to_string(sqrt(w*w + x*x + y*y + z*z)) << endl;
	return make_tuple(w, x, y, z);
}

vector3d RandomMove::getRandomStep()
{
	double theta =  M_PI*randf();
	double phi =  2*M_PI*randf();
	vector3d step(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	return step;
}

tuple<double, double, double, double> RandomMove::hamiltonProduct(tuple<double, double, double, double> a, tuple<double, double, double, double> b)
{
	//double w = get<0>(a)*get<0>(b) - get<1>(a)*get<1>(b) - get<2>(a)*get<2>(b) - get<3>(a)*get<3>(b);
	double w = get<0>(a)*get<0>(b) - get<1>(a)*get<1>(b) - get<2>(a)*get<2>(b) - get<3>(a)*get<3>(b);
	double x = get<0>(a)*get<1>(b) + get<1>(a)*get<0>(b) + get<2>(a)*get<3>(b) - get<3>(a)*get<2>(b);
	double y = get<0>(a)*get<2>(b) - get<1>(a)*get<3>(b) + get<2>(a)*get<0>(b) + get<3>(a)*get<1>(b);
	double z = get<0>(a)*get<3>(b) + get<1>(a)*get<2>(b) - get<2>(a)*get<1>(b) + get<3>(a)*get<0>(b);
	//cout << w << ", " << x << ", " << y << ", " << z << endl;
	return make_tuple(w,x,y,z);
}

vector3d RandomMove::rotateByQuaternion(vector3d inputVector, tuple<double, double, double, double> leftQuaternion)
{
	tuple<double, double, double, double> vec = make_tuple(0, inputVector.x, inputVector.y, inputVector.z);
	tuple<double, double, double, double> rightQuaternion = make_tuple(get<0>(leftQuaternion), -get<1>(leftQuaternion), -get<2>(leftQuaternion), -get<3>(leftQuaternion));
	//cout << "L " <<  get<0>(leftQuaternion) << ", " <<  get<1>(leftQuaternion) << ", " <<  get<2>(leftQuaternion) << ", "  <<  get<3>(leftQuaternion) << endl;
	//cout << "R " <<  get<0>(rightQuaternion) << ", " <<  get<1>(rightQuaternion) << ", " <<  get<2>(rightQuaternion) << ", "  <<  get<3>(rightQuaternion) << endl;
	tuple<double, double, double, double> postTrafo =  hamiltonProduct(hamiltonProduct(leftQuaternion, vec), rightQuaternion);
	vector3d outputVector (get<1>(postTrafo), get<2>(postTrafo), get<3>(postTrafo));
	return outputVector;
}
