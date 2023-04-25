//#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <Eigen\Dense>
#include "dprint.h"
#define PI 3.14159265358979323
double conservativeArcCos(double dot);

class randomNumberGen
{
public:
	randomNumberGen() {
		srand((unsigned)time(NULL));
	};
	~randomNumberGen() {};
public:
	double get() { return rand() * 1.0 / RAND_MAX; }
};

double ICP_Energy(Eigen::Matrix3Xd &X, Eigen::Matrix3Xd &Y);
