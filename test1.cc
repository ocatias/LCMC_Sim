#ifndef TEST
#define TEST

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include "BoxCollision.h"
#include <fstream>
#include <math.h>
#include <tuple>

using namespace std;

const int N = 1000;
const double L = 10;
const int TIMESTEPS = 0;
const double DELTA_X = 1;

box boxes[N];

double randf()
{
	return ((double) rand() / (RAND_MAX));
}


tuple<double, double, double, double> randomQuaternion()
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

vector3d getRandomStep()
{
	double theta =  M_PI*randf();
	double phi =  2*M_PI*randf();
	vector3d step(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	return step;
}

tuple<double, double, double, double> hamiltonProduct(tuple<double, double, double, double> a, tuple<double, double, double, double> b)
{
	//double w = get<0>(a)*get<0>(b) - get<1>(a)*get<1>(b) - get<2>(a)*get<2>(b) - get<3>(a)*get<3>(b);
	double w = get<0>(a)*get<0>(b) - get<1>(a)*get<1>(b) - get<2>(a)*get<2>(b) - get<3>(a)*get<3>(b);
	double x = get<0>(a)*get<1>(b) + get<1>(a)*get<0>(b) + get<2>(a)*get<3>(b) - get<3>(a)*get<2>(b);
	double y = get<0>(a)*get<2>(b) - get<1>(a)*get<3>(b) + get<2>(a)*get<0>(b) + get<3>(a)*get<1>(b);
	double z = get<0>(a)*get<3>(b) + get<1>(a)*get<2>(b) - get<2>(a)*get<1>(b) + get<3>(a)*get<0>(b);
	//cout << w << ", " << x << ", " << y << ", " << z << endl;
	return make_tuple(w,x,y,z);
}

vector3d rotateByQuaternion(vector3d inputVector, tuple<double, double, double, double> leftQuaternion)
{
	tuple<double, double, double, double> vec = make_tuple(0, inputVector.x, inputVector.y, inputVector.z);
	tuple<double, double, double, double> rightQuaternion = make_tuple(get<0>(leftQuaternion), -get<1>(leftQuaternion), -get<2>(leftQuaternion), -get<3>(leftQuaternion));
	//cout << "L " <<  get<0>(leftQuaternion) << ", " <<  get<1>(leftQuaternion) << ", " <<  get<2>(leftQuaternion) << ", "  <<  get<3>(leftQuaternion) << endl;
	//cout << "R " <<  get<0>(rightQuaternion) << ", " <<  get<1>(rightQuaternion) << ", " <<  get<2>(rightQuaternion) << ", "  <<  get<3>(rightQuaternion) << endl;
	tuple<double, double, double, double> postTrafo =  hamiltonProduct(hamiltonProduct(leftQuaternion, vec), rightQuaternion);
	vector3d outputVector (get<1>(postTrafo), get<2>(postTrafo), get<3>(postTrafo));
	return outputVector;
}

vector3d rotateByRandomQuaternion(vector3d inputVector)
{
	return rotateByQuaternion(inputVector, randomQuaternion());
}

void writeStateToFile(ostream &file)
{

	for( int i = 0; i <= N-1; i++)
	{
		file << boxes[i];
	}
	file << endl << endl;
}

bool tryMove(int particleNr, vector3d translationVector)
{
	box trialBox(boxes[particleNr].center+translationVector,
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[3],
		boxes[particleNr].base[0], boxes[particleNr].base[1], boxes[particleNr].base[2]);
	for(int i = 0; i < N; i++)
	{
		if (i == particleNr)
			continue;

		if (BoxCollision::isColliding(trialBox, boxes[i]))
			return false;
	}
	return true;
}

bool tryRotate(int particleNr, tuple<double, double, double, double> quaternion)
{
	box trialBox(boxes[particleNr].center,
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[3],
		rotateByQuaternion(boxes[particleNr].base[0], quaternion),
		rotateByQuaternion(boxes[particleNr].base[1], quaternion),
		rotateByQuaternion(boxes[particleNr].base[2], quaternion));
	for(int i = 0; i < N; i++)
	{
		if (i == particleNr)
			continue;

		if (BoxCollision::isColliding(trialBox, boxes[i]))
			return false;
	}
	return true;
}

int main()
{
	srand(2);
	ofstream file;
	file.open ("output.txt");

	vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);

	int counter = 0;

	for( int i = 0; i < 10; i++)
	{
		for( int j = 0; j < 10; j++)
		{
			for( int k = 0; k < 10; k++)
			{
				boxes[counter] = box(zeroVec + right*i*4.5 + up*j*4.5 + forward*k*4.5, 0.3, 2, 0.3, right, forward, up);
				//cout << counter << " : " <<  to_string(i) << " " << to_string(j) << " " << to_string(k) << endl;
				counter++;
			}
		}
	}
	//ANFANGSWERTE SCHREIBEN
	//writeStateToFile(file);


	//TEST QUATERNION ROTATION, SHOULD RETURN (0,0,-1)
	/*
	vector3d P (1,0,0);
	cout << rotateByQuaternion(P, make_tuple(0.70710678, 0.0,  0.70710678, 0.0));
	*/


	cout << "TEST RANDOM STEP, SHOULD RETURN (0,0,0) : ";
	vector3d sumVec(0,0,0);
	for(int i = 0; i <= 100000; i++)
	{
		sumVec = sumVec + getRandomStep()*DELTA_X;
	}
	cout << sumVec << endl;

	cout << "TEST RANDOM ROTATION, SHOULD RETURN (0,0,0) : ";
	vector3d sumVec2(0,0,0);
	vector3d unit(1,0,0);
	for(int i = 0; i <= 100000; i++)
	{
		sumVec2 = sumVec2 + rotateByQuaternion(unit, randomQuaternion());
	}
	cout << sumVec2 << endl;


	//vector3d vec(1,0,0);
	//vector3d vec2(0,1,0);

	//tuple<double, double, double, double> quaternion {randomQuaternion()};
	//cout << " Quaternion: " << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion)  << " ";
	//cout << rotateByQuaternion(vec, quaternion)<< endl;
	//cout << rotateByQuaternion(vec2, quaternion) << endl;
	//cout << rotateByQuaternion(vec, quaternion) * rotateByQuaternion(vec2, quaternion);

	for(int t = 0; t < TIMESTEPS; t++)
	{
		if(t%100 == 0)
			cout << t << endl;

		int n = rand() % (N);
		//cout << "T: Selected " << to_string(n);
		vector3d translationVector (getRandomStep()*DELTA_X);
		if (tryMove(n, translationVector))
		{
			boxes[n].translate(getRandomStep()*DELTA_X);
			//cout << " Allowed" << endl;
		}
		else
		{
			//cout << " Denied" << endl;
		}

		int n2 = rand() % (N);
		//cout << "R: Selected " << to_string(n2);
		tuple<double, double, double, double> quaternion {randomQuaternion()};
		//cout << " Quaternion: " << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion)  << " ";
		//cout << " Base:" <<  boxes[n2].base[0] << ", "  <<  boxes[n2].base[1] << ", "  <<  boxes[n2].base[2] << " ";
		if (tryRotate(n2, quaternion))
		{
			boxes[n2].base[0] = rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			boxes[n2].base[1] = rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			boxes[n2].base[2] = rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();
			//cout << " Allowed" << endl;
		}
		else
		{
			//cout << " Denied" << endl;
		}

	}
	writeStateToFile(file);
	file.close();
	return 0;
}



#endif
