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
#include <fstream>
#include <valarray>

using namespace std;

int N = 1000;
unsigned long long TIMESTEPS = 2000000000;
//unsigned long long TIMESTEPS = 1000000;

unsigned long long LOGINTERVALL = TIMESTEPS/50;
const double DELTA_X = 0.1;
const double DELTA_ANGLE = 0.5;
//const int UPDATEINTERVALL = 1000000;
const int UPDATEINTERVALL = 10000;

const double R = 18;

double w = 1;
double l = 20;
double h = 1;

/*
const double R = 45;
const double H = 15;

double w = 5;
double l = 12;
double h = 1;
*/

const int relevantBaseIndex = 1;

int acceptedMoves = 0;
int deniedMoves = 0;
int acceptedRotations = 0;
int deniedRotations = 0;
box *boxes;

void writeStateToFile(ostream &file, unsigned long long timestep = 0, valarray<double> S = {0, 0, 0}, valarray<double> B = {0, 0, 0}, valarray<double> deltaS = {0, 0, 0}, double acceptedT = 0, double acceptedR = 0)
{
	file << "sphere, " << R << ", " << R << ", " << timestep << ", " << S[0] << ", " << S[1] << ", " << S[2] << ", " << B[0] << ", " << B[1] << ", " << B[2] <<  ", " << deltaS[0] << ", " << deltaS[1] << ", " << deltaS[2] << ", ";
	file << acceptedT << ", " << acceptedR << endl;
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
		if(sqrt(particle.edges[i].x*particle.edges[i].x + particle.edges[i].y*particle.edges[i].y + particle.edges[i].z*particle.edges[i].z) > R
			||  particle.edges[i].y < 0)
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

/*
	for(int i = 0; i < 8; i++)
	{
		if(trialBox.edges[i].y > H/2 || trialBox.edges[i].y < -H/2)
			cout << "ERROR";
	}
*/

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

int main(int argc, char** argv)
{
	srand(3);
	ofstream fileOut;
	//fileOut.open ("Output/output.txt");

	string fileName = argv[1];
	int j = 2;

	while(FILE *file = fopen(("Output/"+fileName+".txt").c_str(), "r"))
	{
		fclose(file);
		fileName = argv[1] + to_string(j);
		j++;
	}
	fileOut.open("Output/"+fileName+".txt");



	if(argc > 1)
	{
		l = stof(argv[2]);
		h = stof(argv[3]);
		w = stof(argv[4]);

		float density = stof(argv[5]);
		N = density*(R*R*R*4/3*M_PI/2)/(w*l*h);
		cout << "Using input parameters: w = " << w << ", l = " << l << ", h = " << h << ", N = " << N << endl;
	}


	boxes = new box[N];

	double density =  N*w*l*h/(R*R*R*4/3*M_PI/2)*100.0;

	fileOut << N << ", " << relevantBaseIndex << ", " << density << endl;

	vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);

	int counter = 0;
	cout << "Cylinder:  R = " << R << ", H = " << R << ", V = " << R*R*R*4/3*M_PI << endl;
	cout << "Density " << density << "%" << endl;

	cout << "Placing " << N << " particles:" << endl;
	int i {0};
	while(i < N)
	{

/*
		//Completely random positions
		double phi {RandomMove::randf()*2*M_PI};
		double radius {RandomMove::randf()*R};
		double x {cos(phi)*radius};
		double z {sin(phi)*radius};
		double y {RandomMove::randf()*H - H/2};
		box newBox (zeroVec + right*x + up*z + forward*y, w/2, l/2, h/2, right, forward, up);

		tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion()};
		newBox.base[0] = RandomMove::rotateByQuaternion(newBox.base[0], quaternion).normalize();
		newBox.base[1] = RandomMove::rotateByQuaternion(newBox.base[1], quaternion).normalize();
		newBox.base[2] = RandomMove::rotateByQuaternion(newBox.base[2], quaternion).normalize();
*/

/*
		//Upright particles
		double phi {RandomMove::randf()*2*M_PI};
		double radius {RandomMove::randf()*R};
		double x {cos(phi)*radius};
		double z {sin(phi)*radius};
		box newBox (zeroVec + right*x + up*z, w/2, l/2, h/2, right, forward, up);
*/

		//Lying particles
		int maxParticlesY = R/h;

		double yDis = (2*R - maxParticlesY*h)/(maxParticlesY+1);
		double y = 0 +yDis + h/2 + (h+yDis)*(int(RandomMove::randf()*maxParticlesY)  % (maxParticlesY+1));

		double phi = acos(y/R);
		double length = R*sin(phi);
		//cout << zLength << endl;
		int maxParticlesX = 2*(length)/(w + 0.01);
		int maxParticlesZ = 2*(length)/(l + 0.01);

		if ((maxParticlesX <= 0) || (maxParticlesZ <= 0))
			continue;

		double xDis = (2*(length) - 0.01 - maxParticlesX*w)/(maxParticlesX+1);
		double zDis = (2*(length) - 0.01 - maxParticlesZ*l)/(maxParticlesZ+1);
		double x = length - w/2 - zDis -(w+0.01+xDis)*(int(RandomMove::randf()*maxParticlesX)%(maxParticlesX+1)) ;
		double z = length - l/2 - zDis -(l+0.01+zDis)*(int(RandomMove::randf()*maxParticlesZ)%(maxParticlesZ+1)) ;

		box newBox (zeroVec + y*forward + x*right + up*z, w/2, l/2, h/2, right, forward, up);


		//cout << x << " " << y << " " << z << endl;


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

	auto[prevS, _currB] = Statistics::orderParameter(boxes, N);

	writeStateToFile(fileOut, 0, prevS, _currB, {0,0,0});

	unsigned long long currAcceptedTranslations = 0;
	unsigned long long currAcceptedRotations = 0;
	unsigned long long currDeniedTranslations = 0;
	unsigned long long currDeniedRotations = 0;

	auto start = chrono::high_resolution_clock::now();
	cout << flush << endl << "Starting MC with " << TIMESTEPS << " steps:" << endl;
	cout << "Delta X = " << DELTA_X << ", Delta Angle = " << DELTA_ANGLE << endl;
	cout << "S = " << prevS[0] << ", " << prevS[1] << ", " << prevS[2] << endl;
	for(unsigned long long t = 1; t <= TIMESTEPS; t++)
	{
		if(t%UPDATEINTERVALL == 0)
		{
			cout << "% "<< static_cast<long double>(t)/TIMESTEPS*100 << flush <<"\r" << "                  " << "\r";
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
		tuple<double, double, double, double> quaternion {RandomMove::randomQuaternion(DELTA_ANGLE)};
		//tuple<double, double, double> eulerAngles {RandomMove::randomEulerAngles(DELTA_ANGLE)};
		if (tryRotate(n2, quaternion))
		{
			boxes[n2].base[0] = RandomMove::rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			boxes[n2].base[1] = RandomMove::rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			boxes[n2].base[2] = RandomMove::rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();
			//boxes[n2].base[0] = RandomMove::rotateByEulerAngles(boxes[n2].base[0], eulerAngles).normalize();
			//boxes[n2].base[1] = RandomMove::rotateByEulerAngles(boxes[n2].base[1], eulerAngles).normalize();
			//boxes[n2].base[2] = RandomMove::rotateByEulerAngles(boxes[n2].base[2], eulerAngles).normalize();

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
    	auto [ currS, currB ] = Statistics::orderParameter(boxes, N);
			cout << std::left <<  "S = " << std::setw(12) << currS[0]
			 	<< "ΔS = " << std::setw(13) << (currS-prevS)[0]
				<< "T% = " << std::setw(9) << static_cast<long double>(currAcceptedTranslations)/(currAcceptedTranslations + currDeniedTranslations)*100
				<< "R% = " << std::setw(9) << static_cast<long double>(currAcceptedRotations)/(currAcceptedRotations+currDeniedRotations)*100 << endl;
			writeStateToFile(fileOut, t, currS, currB, currS-prevS,
				static_cast<long double>(currAcceptedTranslations)/(currAcceptedTranslations + currDeniedTranslations)*100,
				static_cast<long double>(currAcceptedRotations)/(currAcceptedRotations+currDeniedRotations)*100);
			prevS = currS;
		}
	}
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(stop - start);

	auto [ currS, currB ] = Statistics::orderParameter(boxes, N);
	writeStateToFile(fileOut, TIMESTEPS, currS, currB, currS-prevS,
		static_cast<long double>(currAcceptedTranslations)/(currAcceptedTranslations + currDeniedTranslations)*100,
		static_cast<long double>(currAcceptedRotations)/(currAcceptedRotations+currDeniedRotations)*100);
	cout << endl;
	cout << "Time for MC simulation " << duration.count() << " seconds" << endl;
	cout << "Accepted Moves: " << static_cast<long double>(acceptedMoves)/(acceptedMoves+deniedMoves)*100 << "%" << endl;
	cout << "Accepted Rotations: " << static_cast<long double>(acceptedRotations)/(acceptedRotations+deniedRotations)*100 << "%" << endl;

	//cout << acceptedMoves << "/" << acceptedMoves+deniedMoves << endl;
	//cout << acceptedRotations << "/" << acceptedRotations+deniedRotations << endl;

	//fileIn.close();
	fileOut.close();
	return 0;
}

#endif
