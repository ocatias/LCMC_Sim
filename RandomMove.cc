#include "RandomMove.h"

double RandomMove::randf()
{
	return ((double) rand() / (RAND_MAX));
}


tuple<double, double, double, double> RandomMove::randomQuaternion(double delta)
{
	//if(delta > 1)
	//	delta = 1;

	/*
	std::random_device rd{};
  std::mt19937 gen{rd()};

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> d{0,1};
 	double x = d(gen);
	double y = d(gen);
	double z = d(gen);
	double w = d(gen);
	double l = sqrt(w*w + x*x + y*y + z*z);
	return make_tuple(w/l, x/l, y/l, z/l);
	*/

	double s = randf()*delta;
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2*M_PI*randf();
	double theta2 = 2*M_PI*randf()*delta;
	double w = cos(theta2)*sigma2;
	double x = sin(theta1)*sigma1;
	double y = cos(theta1)*sigma1;
	double z = sin(theta2)*sigma2;
	//cout << w << ", " << x << ", " << y << ", " << z << endl;
	//cout << "Quaternion size = " << to_string(sqrt(w*w + x*x + y*y + z*z)) << endl;
	return make_tuple(w, x, y, z);
}

//delta is the maximum of each individual rotation in percent
//Propably slightly wrong
tuple<double, double, double> RandomMove::randomEulerAngles(double delta)
{
	assert(!((delta > 1) ||(delta <= 0)));
	double theta = 2*M_PI*randf() - M_PI;
	double phi = acos(1 - 2*randf()) + M_PI/2.;
	if(randf() < 0.5)
	{
		if(phi < M_PI)
			phi = phi + M_PI;
		else
			phi = phi - M_PI;
	}
	double eta = 2*M_PI*randf() - M_PI;
	return make_tuple(theta, phi, eta);
	//return make_tuple(randf()*delta*2*M_PI, randf()*delta*2*M_PI, randf()*delta*2*M_PI);
}

//Propably slightly wrong
vector3d RandomMove::rotateByEulerAngles(vector3d vec, double gamma, double beta, double alpha)
{
	//return vector3d(
	//	(cos(alpha)*cos(gamma) - sin(alpha)*cos(beta)*sin(gamma))*vec.x + (-cos(alpha)*sin(gamma) - sin(alpha)*cos(beta)*cos(gamma))*vec.y + sin(alpha)*sin(beta)*vec.z,
	//	(sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*sin(gamma))*vec.x + (-sin(alpha)*sin(gamma) + cos(alpha)*cos(beta)*cos(gamma))*vec.y - cos(alpha)*sin(beta)*vec.z,
	//	 sin(beta)*sin(gamma)*vec.x + sin(beta)*cos(gamma)*vec.y + cos(beta)*vec.z);
	 return vector3d(
			 cos(alpha)*cos(beta)*vec.x + (cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma))*vec.y + (cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma))*vec.z,
			 sin(alpha)*cos(beta)*vec.x + (sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma))*vec.y + (sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma))*vec.z,
				-sin(beta)*vec.x + cos(beta)*sin(gamma)*vec.y + cos(beta)*cos(gamma)*vec.z);

}

//Propably slightly wrong
vector3d RandomMove::rotateByEulerAngles(vector3d vec, tuple<double, double, double> eulerAngles)
{
	double alpha = get<0>(eulerAngles);
	double beta = get<1>(eulerAngles);
	double gamma = get<2>(eulerAngles);
	return rotateByEulerAngles(vec, alpha, beta, gamma);
}


vector3d RandomMove::getRandomStep()
{
	double x, y, z, length;
	do {
		x = (RandomMove::randf() - 0.5)*2;
		y = (RandomMove::randf() - 0.5)*2;
		z = (RandomMove::randf() - 0.5)*2;
		length = sqrt(x*x + y*y + z*z);
	} while(length > 1);

	x /= length;
	y /= length;
	z /= length;

	return vector3d(x, y, z);
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
	tuple<double, double, double, double> postTrafo =  hamiltonProduct(leftQuaternion, hamiltonProduct(vec, rightQuaternion));
	vector3d outputVector (get<1>(postTrafo), get<2>(postTrafo), get<3>(postTrafo));
	return outputVector;
}
