#ifndef CYLINDER_MC
#define CYLINDER_MC

#include <iostream>
#include "vector3d.h"
#include "box.h"
#include "BoxCollision.h"
#include <fstream>
#include <math.h>
#include <tuple>
#include "RandomMove.h"
#include <chrono>
#include "Statistics.h"
#include <iomanip>

using namespace std;

const int N = 42;
const int TIMESTEPS = 50000000;
const int LOGINTERVALL = TIMESTEPS/10;
const double DELTA_X = 0.05;
const double DELTA_ANGLE = 0.001;

const double R = 5;
const double H = 5;

double w = 0.5;
double h = 1;
double l = 1;

const int relevantBaseIndex = 0;

int acceptedMoves = 0;
int deniedMoves = 0;
int acceptedRotations = 0;
int deniedRotations = 0;
box boxes[N];

void writeStateToFile(ostream &file, int timestep = 0, double S = 0, double deltaS = 0, double acceptedT = 0, double acceptedR = 0)
{
	file << "cylinder, " << R << ", " << H << ", " << timestep << ", " << S << ", " << deltaS << ", " << acceptedT << ", " << acceptedR << endl;
	for( int i = 0; i <= N-1; i++)
	{
		file << boxes[i];
	}
}

bool isOutsideCylinder(box particle)
{
	particle.updateEdges();
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

bool tryRotate(int particleNr, tuple<double, double, double> eulerAngles)
{
	box trialBox(boxes[particleNr].center,
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[2],
		RandomMove::rotateByEulerAngles(boxes[particleNr].base[0], eulerAngles),
		RandomMove::rotateByEulerAngles(boxes[particleNr].base[1], eulerAngles),
		RandomMove::rotateByEulerAngles(boxes[particleNr].base[2], eulerAngles));
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
	srand(3);
	ofstream fileOut;
	fileOut.open ("Output/output.txt");
	fileOut << N << ", " << relevantBaseIndex << endl;

	vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);

	int counter = 0;
	cout << "Cylinder:  R = " << R << ", H = " << H << ", V = " << R*R*M_PI << endl;
	cout << "Density " << N*w*l*h/(R*R*H*M_PI)*100.0 << "%" << endl;

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
			cout << i << flush << "\r";
		}
	}

	double prevS = Statistics::orderParameter(boxes, N, relevantBaseIndex);
	writeStateToFile(fileOut, 0, prevS, 0);

	int currAcceptedTranslations = 0;
	int currAcceptedRotations = 0;
	int currDeniedTranslations = 0;
	int currDeniedRotations = 0;

	auto start = chrono::high_resolution_clock::now();
	cout << flush << endl << "Starting MC with " << TIMESTEPS << " steps:" << endl;
	for(int t = 1; t <= TIMESTEPS; t++)
	{
		if(t%10000 == 0)
		{
			cout << "% "<< float(t)/TIMESTEPS*100 << flush <<"\r" << "                  " << "\r";
		}

		int n = rand() % (N);
		vector3d translationVector (RandomMove::getRandomStep()*DELTA_X);
		if (tryMove(n, translationVector))
		{
			boxes[n].translate(translationVector);
			acceptedMoves++;
			currAcceptedTranslations++;
		}
		else
		{
			deniedMoves++;
			currDeniedTranslations++;
		}

		int n2 = rand() % (N);
		//tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
		tuple<double, double, double> eulerAngles {RandomMove::randomEulerAngles(DELTA_ANGLE)};
		if (tryRotate(n2, eulerAngles))
		{
			//boxes[n2].base[0] = RandomMove::rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			//boxes[n2].base[1] = RandomMove::rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			//boxes[n2].base[2] = RandomMove::rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();
			boxes[n2].base[0] = RandomMove::rotateByEulerAngles(boxes[n2].base[0], eulerAngles).normalize();
			boxes[n2].base[1] = RandomMove::rotateByEulerAngles(boxes[n2].base[1], eulerAngles).normalize();
			boxes[n2].base[2] = RandomMove::rotateByEulerAngles(boxes[n2].base[2], eulerAngles).normalize();

			boxes[n2].updateEdges();
			acceptedRotations++;
			currAcceptedRotations++;
		}
		else
		{
			deniedRotations++;
			currDeniedRotations++;
		}

		if((t%LOGINTERVALL == 0) && (t != TIMESTEPS))
		{
			double currS = Statistics::orderParameter(boxes, N, relevantBaseIndex);
			cout << std::left <<  "S = " << std::setw(20) << currS;
			cout  << "ΔS = " << currS-prevS  << endl ;
			writeStateToFile(fileOut, t, currS, currS-prevS,
				double(currAcceptedTranslations)/(currAcceptedTranslations + currDeniedTranslations)*100,
				double(currAcceptedRotations)/(currAcceptedRotations+currDeniedRotations)*100);
			prevS = currS;
		}
	}
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(stop - start);

	double currS = Statistics::orderParameter(boxes, N, relevantBaseIndex);
	writeStateToFile(fileOut, TIMESTEPS, currS, currS-prevS,
		double(currAcceptedTranslations)/(currAcceptedTranslations + currDeniedTranslations)*100,
		double(currAcceptedRotations)/(currAcceptedRotations+currDeniedRotations)*100);
	cout << endl;
	cout << "Time for MC simulation " << duration.count() << " seconds" << endl;
	cout << "Accepted Moves: " << float(acceptedMoves)/(acceptedMoves+deniedMoves)*100 << "%" << endl;
	cout << "Accepted Rotations: " << float(acceptedRotations)/(acceptedRotations+deniedRotations)*100 << "%" << endl;

	//cout << acceptedMoves << "/" << acceptedMoves+deniedMoves << endl;
	//cout << acceptedRotations << "/" << acceptedRotations+deniedRotations << endl;

	//fileIn.close();
	fileOut.close();
	return 0;
}

#endif
