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
	//polygon������Ϊpolygon.rows(),Ҫ���������ʱ������
	Eigen::Vector2d ComputePolygonInteriorPoint(Polygon &polygon);
	bool IsInPolygon(Polygon &polygon, Eigen::Vector2d &p);
	double ComputePolygonArea(Polygon &polygon);
	double ComputePolygonPerimeter(Polygon &polygon);
	Eigen::Vector4d ComputePolygonBound(Polygon &polygon);//����(x_min,x_max,y_min,y_max)

	int ComputePolygonsOuterBoundIndex(Polygons &polygons);//������������е���Ȧ���
	//�������ɶ���Σ�Ҫ���һ������Ȧ���������Ȧ���ཻ
	//Ψһ��Ȧ������ʱ��������Ȧ��Ϊ˳ʱ������
	bool IsInPolygons(Polygons &polygons, Eigen::Vector2d &p);
	double ComputePolygonsArea(Polygons &polygons);//��Ȧ��ռ�����ȥ��Ȧ��ռ���
#pragma endregion

#ifdef OPENMESH_TRIMESH_ARRAY_KERNEL_HH
	double ComputeVoronoiArea(TriMesh* mesh, OpenMesh::VertexHandle v);
	double ComputeGaussCurvature(TriMesh* mesh, OpenMesh::VertexHandle v);//�����˹����
	double ComputeMeanCurvature(TriMesh* mesh, OpenMesh::VertexHandle v);
#endif OPENMESH_TRIMESH_ARRAY_KERNEL_HH
}