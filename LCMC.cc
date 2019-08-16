#ifndef LCMC_C
#define LCMC_C

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

//Number of timesteps
unsigned long long TIMESTEPS = 200000;


const double DELTA_X = 0.6;
const double DELTA_ANGLE = 0.06;

//How many different frames of the simulation will be saved
unsigned long long LOGINTERVALL = TIMESTEPS/10;
//After how many timesteps the progress-counter will be updated
const int UPDATEINTERVALL = 100000;

const int relevantBaseIndex = 1;

double R = 19;
double H = 9;

double w = 1;
double l = 20;
double h = 3;

double volume;

int acceptedMoves = 0;
int deniedMoves = 0;
int acceptedRotations = 0;
int deniedRotations = 0;
int N = 1000;

string confinementName = "cylinder";

box *boxes;

void writeStateToFile(ostream &file, unsigned long long timestep = 0, valarray<double> S = {0, 0, 0}, valarray<double> B = {0, 0, 0}, valarray<double> deltaS = {0, 0, 0}, double acceptedT = 0, double acceptedR = 0)
{
	file << confinementName << ", " << R << ", " << H << ", " << timestep << ", " << S[0] << ", " << S[1] << ", " << S[2] << ", " << B[0] << ", " << B[1] << ", " << B[2] <<  ", " << deltaS[0] << ", " << deltaS[1] << ", " << deltaS[2] << ", ";
	file << acceptedT << ", " << acceptedR << endl;
	for( int i = 0; i <= N-1; i++)
		file << boxes[i];
}

bool (*isOutside)(box particle);

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

bool isOutsideCube(box particle)
{
	particle.updateEdges();
	for(int i = 0; i < 8; i++)
	{
		if(abs(particle.edges[i].x) > H/2	|| abs(particle.edges[i].y) > H/2 || abs(particle.edges[i].z) > H/2)
			return true;
	}
	return false;
}

bool isOutSideSphere(box particle)
{
	particle.updateEdges();
	for(int i = 0; i < 8; i++)
	{
		if(sqrt(particle.edges[i].x*particle.edges[i].x + particle.edges[i].y*particle.edges[i].y + particle.edges[i].z*particle.edges[i].z) > R)
			return true;
	}
	return false;
}

bool isOutsideHalfSphere(box particle)
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

void (*initializeParticles)(int N);


void initializeParticlesCylinder(int N)
{
  vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);
  int i {0};
  while(i < N)
  {
  /*
    //Upright particles
    double phi {RandomMove::randf()*2*M_PI};
    double radius {RandomMove::randf()*R};
    double x {cos(phi)*radius};
    double z {sin(phi)*radius};
    box newBox (zeroVec + right*x + up*z, w/2, l/2, h/2, right, forward, up);
  */
    //Lying particles
    //double phi {RandomMove::randf()*2*M_PI};
    //double radius {RandomMove::randf()*R};
    //double x {cos(phi)*radius};
    //double z {sin(phi)*radius};
    //int maxParticlesOnTop = H/(h+0.1);
    int maxParticlesOnTop = H/h;

    int xParticles = 2*R/(w+0.3);
    //int yParticles = 2*R/(h+0.0001);

    //cout << -H/2 + h/2 + 0.01 + (h+0.01)*(int(RandomMove::randf()*maxParticlesOnTop)  % (maxParticlesOnTop+1)) << endl;
    //double yDis = (H - 0.1 - maxParticlesOnTop*h)/(maxParticlesOnTop+1);
    //double y = -H/2 +yDis + h/2 + 0.1 + (h+0.1+yDis)*(int(RandomMove::randf()*maxParticlesOnTop)  % (maxParticlesOnTop+1));
    double yDis = (H - maxParticlesOnTop*h)/(maxParticlesOnTop+1);
    double y = -H/2 +yDis + h/2 + (h+yDis)*(int(RandomMove::randf()*maxParticlesOnTop)  % (maxParticlesOnTop+1));

    double x = -R + (w+0.3)*(int(RandomMove::randf()*xParticles)%(xParticles+1)) ;
    double phi = acos(x/R);
    double zLength = R*sin(phi);
    //cout << zLength << endl;
    int zParticles = 2*(zLength - 1)/(l + 0.1);
    if (zParticles <= 0)
      continue;

    double zDis = (2*(zLength) -1 - 0.1 - zParticles*l)/(zParticles+1);
    double z = zLength - l/2 - 0.5 - zDis -(l+0.1+zDis)*(int(RandomMove::randf()*zParticles)%(zParticles+1)) ;

    //cout <<(h+0.00001)*(int((RandomMove::randf()-0.5)*maxParticlesOnTop)  % maxParticlesOnTop)<< endl << endl;
    box newBox (zeroVec + y*forward + x*right + up*z, w/2, l/2, h/2, right, up, forward);

    bool isAllowed = true;
    for(int j = 0; j < i; j++)
    {
      if(BoxCollision::isColliding(newBox, boxes[j]))
        isAllowed = false;
    }

    if(!isOutside(newBox) && isAllowed)
    {
      boxes[i] = newBox;
      i++;
      cout << i << flush << "\r";
    }
  }
}

void initializeParticlesCube(int N)
{
  vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);
	int i {0};
	while(i < N)
	{
		int maxParticlesL = H/w;
		int maxParticlesW = H/l;
		int maxParticlesH = H/h;

		double xDis = (H - maxParticlesL*w)/(maxParticlesL+1) + 0.01;
		double yDis = (H - maxParticlesW*l)/(maxParticlesW+1) + 0.01;
		double zDis = (H - maxParticlesH*h)/(maxParticlesH+1) + 0.01;

		double x = -H/2 + xDis  + w/2 + (w+xDis)*(int(RandomMove::randf()*maxParticlesL)  % (maxParticlesL+1));
		double y = -H/2 + yDis + l/2 + (l+yDis)*(int(RandomMove::randf()*maxParticlesW)  % (maxParticlesW+1));
		double z = -H/2 + zDis + h/2 + (h+zDis)*(int(RandomMove::randf()*maxParticlesH)  % (maxParticlesH+1));

		box newBox (zeroVec + y*forward + x*right + up*z, w/2, l/2, h/2, right, forward, up);

		bool isAllowed = true;
		for(int j = 0; j < i; j++)
		{
			if(BoxCollision::isColliding(newBox, boxes[j]))
				isAllowed = false;
		}

		if(!isOutside(newBox) && isAllowed)
		{
			boxes[i] = newBox;
			i++;
			cout << i << flush << "\r";
		}
	}
}

void initializeParticlesSphere(int N)
{
  vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);
	int i {0};
	while(i < N)
	{
		//Lying particles
		int maxParticlesY = 2*R/h;

		double yDis = (2*R - maxParticlesY*h)/(maxParticlesY+1);
		double y = -R +yDis + h/2 + (h+yDis)*(int(RandomMove::randf()*maxParticlesY)  % (maxParticlesY+1));

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

}

void initializeParticlesHalfSphere(int N)
{
  vector3d zeroVec(0, 0, 0);
	vector3d up(0, 0, 1);
	vector3d forward(0, 1, 0);
	vector3d right(1, 0, 0);
	int i {0};
	while(i < N)
	{
		//Lying particles
		int maxParticlesY = R/h;

		double yDis = (R - maxParticlesY*h)/(maxParticlesY+1);
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
}

bool tryMove(int particleNr, vector3d translationVector)
{
	box trialBox(boxes[particleNr].center+translationVector,
		boxes[particleNr].halfRatio[0], boxes[particleNr].halfRatio[1], boxes[particleNr].halfRatio[2],
		boxes[particleNr].base[0], boxes[particleNr].base[1], boxes[particleNr].base[2]);

	if(isOutside(trialBox))
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

	if(isOutside(trialBox))
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

	//Creating the output file
	ofstream fileOut;
	string fileName = argv[1];
	int j = 2;
	while(FILE *file = fopen(("Output/"+fileName+".txt").c_str(), "r"))
	{
		fclose(file);
		fileName = argv[1] + to_string(j);
		j++;
	}
	fileOut.open("Output/"+fileName+".txt");

	l = stof(argv[2]);
	h = stof(argv[3]);
	float density = stof(argv[4]);

	//Confinement selection
  switch (stoi(argv[5]))
  {
    //Cylinder
    case 0:
      isOutside = isOutsideCylinder;
      R = stof(argv[6]);
      H = stof(argv[7]);
      volume = R*R*M_PI*H;
      initializeParticles = initializeParticlesCylinder;
      confinementName = "cylinder";
			cout << "Cylinder:  R = " << R << ", H = " << H << ", V = " << volume << endl;
      break;

    //Cube
    case 1:
      isOutside = isOutsideCube;
      H = stof(argv[6]);
      R = stof(argv[6]);
      volume = H*H*H;
      initializeParticles = initializeParticlesCube;
      confinementName = "cube";
			cout << "Cube:  L = " << R << ", V = " << volume << endl;
      break;

    //Sphere
    case 2:
      isOutside = isOutSideSphere;
      R = stof(argv[6]);
      volume = R*R*R*4./3*M_PI;
      initializeParticles = initializeParticlesSphere;
      confinementName = "sphere";
			cout << "Sphere:  R = " << R << ", V = " << volume << endl;
      break;

    //Half Sphere
    case 3:
      isOutside = isOutsideHalfSphere;
      R = stof(argv[6]);
      volume = R*R*R*2./3*M_PI;
      initializeParticles = initializeParticlesSphere;
      confinementName = "sphere";
			cout << "Half Sphere:  R = " << R << ", V = " << volume << endl;
      break;
  }

	N = density*volume/(w*l*h);
	cout << "Input parameters: w = " << w << ", l = " << l << ", h = " << h << ", N = " << N << endl;

	boxes = new box[N];

	double actual_density =  N*w*l*h/volume*100.0;

	fileOut << N << ", " << relevantBaseIndex << ", " << actual_density << endl;

	int counter = 0;
	cout << "Density " << actual_density << "%" << endl;
	cout << "Placing " << N << " particles:" << endl;

  initializeParticlesCylinder(N);

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
		if (tryRotate(n2, quaternion))
		{
			boxes[n2].base[0] = RandomMove::rotateByQuaternion(boxes[n2].base[0], quaternion).normalize();
			boxes[n2].base[1] = RandomMove::rotateByQuaternion(boxes[n2].base[1], quaternion).normalize();
			boxes[n2].base[2] = RandomMove::rotateByQuaternion(boxes[n2].base[2], quaternion).normalize();

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

	fileOut.close();
	return 0;
}

#endif
