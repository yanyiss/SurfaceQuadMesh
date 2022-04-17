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
	~regularCrossfieldGenerator() {};
public:
	void run(std::vector<OpenMesh::Vec3d>& crossfield);
	void calculateMeshFaceBase();

private:
	TriMesh* mesh;
	std::vector<OpenMesh::Vec3d> faceBase;
};

