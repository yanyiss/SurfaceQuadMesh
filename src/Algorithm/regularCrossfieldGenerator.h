#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include <complex>
#include <Eigen\Sparse>

#define PI 3.14159265358979323
#include "..\Toolbox\dprint.h"
#include "..\Toolbox\filesOperator.h"
class regularCrossfieldGenerator
{
public:
	regularCrossfieldGenerator(TriMesh* m) : mesh(m) { calculateMeshFaceBase(); };
	~regularCrossfieldGenerator(std::vector<OpenMesh::Vec3d>& crossfield);
public:
	void run();
	void calculateMeshFaceBase();

private:
	TriMesh* mesh;
	std::vector<OpenMesh::Vec3d> faceBase;
};

