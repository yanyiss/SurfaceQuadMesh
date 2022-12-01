#include "mathFunctions.h"
double conservativeArcCos(double dot)
{
	if (dot > 1.0) return 0.0;
	else if (dot < -1.0) return PI;
	else return acos(dot);
}