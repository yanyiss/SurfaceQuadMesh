#pragma once
#include <Eigen/Core>
#include "..\src\Dependency\triangle\triangle.h"
#include <sstream>
#include <iomanip>
#include "GeneralMathMethod.h"

void triangulate(const Eigen::Matrix2Xd &all_pts, const Eigen::Matrix2Xi &E, const double area_threshold, TriMesh &mesh);

void triangulate(const Eigen::Matrix2Xd &all_pts, const std::vector<Eigen::Matrix2Xi> &E, const double area_threshold, TriMesh &mesh);