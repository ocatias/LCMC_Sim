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
#include <chrono>

using namespace std;

const int N = 300;
//const double L = 10;
const int TIMESTEPS = 10000000;
const double DELTA_X = 0.3;

const double R = 10;
const double H = 10;

double w = 0.5;
double l = 2;
double h = 2;

int acceptedMoves = 0;
int deniedMoves = 0;
int acceptedRotations = 0;
int deniedRotations = 0;
box boxes[N];

void writeStateToFile(ostream &file)
{

	for( int i = 0; i <= N-1; i++)
	{
		file << boxes[i];
	}
	file << endl << endl;
}

bool isOutsideCylinder(box particle)
{
	for(int i = 0; i < 8; i++)
	{
		if(sqrt(particle.edges[i].x*particle.edges[i].x + particle.edges[i].z*particle.edges[i].z) > R
		|| abs(particle.edges[i].y) > H/2)
			return true;
	}
	return false;
}

bool tryMove(int particleNr, vector3d translationVector)
{
	box trialBox(boxes[particleNr].center+translationVector,
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[2],
		boxes[particleNr].base[0], boxes[particleNr].base[1], boxes[particleNr].base[2]);

	if(isOutsideCylinder(trialBox))
		return false;

	for(int i = 0; i < 8; i++)
	{
		if(trialBox.edges[i].y > 5 || trialBox.edges[i].y < -5)
			cout << "ERROR";
	}

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
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[2],
		RandomMove::rotateByQuaternion(boxes[particleNr].base[0], quaternion),
		RandomMove::rotateByQuaternion(boxes[particleNr].base[1], quaternion),
		RandomMove::rotateByQuaternion(boxes[particleNr].base[2], quaternion));
	trialBox.updateEdges();

	if(isOutsideCylinder(trialBox))
			return false;

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
	ofstream fileOut;
	fileOut.open ("Output/output.txt");
	ofstream fileIn;
	fileIn.open ("Output/input.txt");

	vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);

	int counter = 0;

	cout << "Placing " << N << " particles:" << endl;
	int i {0};
	while(i < N)
	{
		double phi {RandomMove::randf()*2*M_PI};
		double radius {RandomMove::randf()*R};
		double x {cos(phi)*radius};
		double z {sin(phi)*radius};
		double y {RandomMove::randf()*H - H/2};
		box newBox (zeroVec + right*x + up*z + forward*y, w/2, h/2, l/2, right, forward, up);

		tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
		newBox.base[0] = RandomMove::rotateByQuaternion(newBox.base[0], quaternion).normalize();
		newBox.base[1] = RandomMove::rotateByQuaternion(newBox.base[1], quaternion).normalize();
		newBox.base[2] = RandomMove::rotateByQuaternion(newBox.base[2], quaternion).normalize();

		bool isAllowed = true;
		for(int j = 0; j < i; j++)
		{
			if(BoxCollision::isColliding(newBox, boxes[j]))
				isAllowed = false;
		}

		if(!isOutsideCylinder(newBox) && isAllowed)
		{
			boxes[i] = newBox;
			i++;
			cout << i << "\r" << flush;
		}
	}

	writeStateToFile(fileIn);

	auto start = chrono::high_resolution_clock::now();

	cout << flush << endl << "Starting MC with " << TIMESTEPS << " steps:" << endl;
	for(int t = 1; t <= TIMESTEPS; t++)
	{
		if(t%1000 == 0)
		{
			cout << t << "\r" << flush;
		}

		int n = rand() % (N);
		vector3d translationVector (RandomMove::getRandomStep()*DELTA_X);
		if (tryMove(n, translationVector))
		{
			boxes[n].translate(translationVector);
			acceptedMoves++;
		}
		else
			deniedMoves++;

		int n2 = rand() % (N);
		tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
		if (tryRotate(n2, quaternion))
		{
			boxes[n2].base[0] = RandomMove::rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			boxes[n2].base[1] = RandomMove::rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			boxes[n2].base[2] = RandomMove::rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();
			boxes[n2].updateEdges();
			acceptedRotations++;
		}
		else
			deniedRotations++;
	}
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(stop - start);

	writeStateToFile(fileOut);
	cout << endl;
	cout << "Time for MC simulation " << duration.count() << " seconds" << endl;
	cout << "Accepted Moves: " << float(acceptedMoves)/(acceptedMoves+deniedMoves)*100 << "%" << endl;
	cout << "Accepted Rotations: " << float(acceptedRotations)/(acceptedRotations+deniedRotations)*100 << "%" << endl;

	//cout << acceptedMoves << "/" << acceptedMoves+deniedMoves << endl;
	//cout << acceptedRotations << "/" << acceptedRotations+deniedRotations << endl;

	fileIn.close();
	fileOut.close();
	return 0;
}

#endif
