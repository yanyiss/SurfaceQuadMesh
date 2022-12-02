//#include <stdlib.h>
#include <cmath>
#include <ctime>
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