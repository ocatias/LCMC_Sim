#ifndef TEST
#define TEST

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include "BoxCollision.h"
#include <fstream>
#include <math.h>
#include <tuple>
#include "RandomMove.h"

using namespace std;

const int N = 1000;
const double L = 10;
const int TIMESTEPS = 1000;
const double DELTA_X = 1;

box boxes[N];

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
		RandomMove::rotateByQuaternion(boxes[particleNr].base[0], quaternion),
		RandomMove::rotateByQuaternion(boxes[particleNr].base[1], quaternion),
		RandomMove::rotateByQuaternion(boxes[particleNr].base[2], quaternion));
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

	//vector3d vec(1,0,0);
	//vector3d vec2(0,1,0);

	//tuple<double, double, double, double> quaternion {randomQuaternion()};
	//cout << " Quaternion: " << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion)  << " ";
	//cout << rotateByQuaternion(vec, quaternion)<< endl;
	//cout << rotateByQuaternion(vec2, quaternion) << endl;
	//cout << rotateByQuaternion(vec, quaternion) * rotateByQuaternion(vec2, quaternion);

	for(int t = 1; t <= TIMESTEPS; t++)
	{
		if(t%100 == 0)
			cout << t << endl;

		int n = rand() % (N);
		//cout << "T: Selected " << to_string(n);
		vector3d translationVector (RandomMove::getRandomStep()*DELTA_X);
		if (tryMove(n, translationVector))
		{
			boxes[n].translate(RandomMove::getRandomStep()*DELTA_X);
			//cout << " Allowed" << endl;
		}


		int n2 = rand() % (N);
		//cout << "R: Selected " << to_string(n2);
		tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
		//cout << " Quaternion: " << get<0>(quaternion) << ", " << get<1>(quaternion) << ", " << get<2>(quaternion) << ", " << get<3>(quaternion)  << " ";
		//cout << " Base:" <<  boxes[n2].base[0] << ", "  <<  boxes[n2].base[1] << ", "  <<  boxes[n2].base[2] << " ";
		if (tryRotate(n2, quaternion))
		{
			boxes[n2].base[0] = RandomMove::rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			boxes[n2].base[1] = RandomMove::rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			boxes[n2].base[2] = RandomMove::rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();
			//cout << " Allowed" << endl;
		}

	}
	writeStateToFile(file);
	file.close();
	return 0;
}

#endif
