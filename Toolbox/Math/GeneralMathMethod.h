#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "..\src\MeshViewer\MeshDefinition.h"

namespace GeneralMathMethod
{
#define PI 3.1415926535897932

	using namespace Eigen;
	using namespace std;
	typedef Matrix2Xd Polygon;
	typedef vector<Polygon> Polygons;

#pragma region polygon method
	//polygon顶点数为polygon.rows(),要求输入点逆时针排序
	Eigen::Vector2d ComputePolygonInteriorPoint(Polygon &polygon);
	bool IsInPolygon(Polygon &polygon, Eigen::Vector2d &p);
	double ComputePolygonArea(Polygon &polygon);
	double ComputePolygonPerimeter(Polygon &polygon);
	Eigen::Vector4d ComputePolygonBound(Polygon &polygon);//返回(x_min,x_max,y_min,y_max)

	int ComputePolygonsOuterBoundIndex(Polygons &polygons);//计算多个多边形中的外圈序号
	//以下若干多边形，要求第一个是外圈，各多边形圈不相交
	//唯一外圈各点逆时针排序，内圈皆为顺时针排序
	bool IsInPolygons(Polygons &polygons, Eigen::Vector2d &p);
	double ComputePolygonsArea(Polygons &polygons);//外圈所占面积减去内圈所占面积
#pragma endregion

#ifdef OPENMESH_TRIMESH_ARRAY_KERNEL_HH
	double ComputeVoronoiArea(TriMesh* mesh, OpenMesh::VertexHandle v);
	double ComputeGaussCurvature(TriMesh* mesh, OpenMesh::VertexHandle v);//计算高斯曲率
	double ComputeMeanCurvature(TriMesh* mesh, OpenMesh::VertexHandle v);
#endif OPENMESH_TRIMESH_ARRAY_KERNEL_HH
}