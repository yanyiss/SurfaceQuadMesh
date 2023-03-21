#include "mathFunctions.h"
double conservativeArcCos(double dot)
{
	if (dot > 1.0) return 0.0;
	else if (dot < -1.0) return PI;
	else return acos(dot);
}

double ICP_Energy(Eigen::Matrix3Xd &X, Eigen::Matrix3Xd &Y)
{
	Eigen::Matrix3d H; H.setZero();
	int n = X.cols();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				H(j, k) += X(j, i)*Y(k, i);
	H /= n;
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U, V;
	U = svd.matrixU(); V = svd.matrixV();
	if (U.determinant()*V.determinant() < 0)
	{
		for (int i = 0; i < 3; ++i)
			U(i, 2) = -U(i, 2);
	}
	H = V * U.transpose();
	double sum = 0;
	for (int i = 0; i < n; ++i)
	{
		sum += (H*X.col(i) - Y.col(i)).squaredNorm();
	}
	return sum / n;
}