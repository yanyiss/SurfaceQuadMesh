#include "..\MeshViewer\MeshDefinition.h"
#include <complex>
#include <Eigen\Sparse>

#define PI 3.14159265358979323
#include "..\Toolbox\dprint.h"
#include "..\Toolbox\filesOperator.h"
#include "StatisticsMostValues.h"
/*
generate crossfield for quad layout
should satisfy the following properties
1. align with principal direction at least in most regions
2. singularities should not appear in where smaller principal curvature is almost zero
3. smooth in most regions but at singularities
4. may decided later in umblic regions

*/
class crossField
{
public:
	crossField(Mesh *m) : mesh(m) {
		position.resize(3, m->n_vertices());
		for (auto &tv : m->vertices())
		{
			auto &p = mesh->point(tv);
			position.col(tv.idx()) << p[0], p[1], p[2];
		}
		calculateMeshFaceBase(); 
	};
	crossField(const crossField& cf) = delete;
	~crossField() {};

private:
	void calculateMeshFaceBase();
public:
	void runPolynomial();
	//void runIteration(std::vector<OpenMesh::Vec3d>& crossfield);

	void setCurvatureConstraint();
	Eigen::VectorXi& getConstraintFace() {
		return constraintId;
	}
	Eigen::Matrix3Xd& getCrossField()
	{
		return crossfield;
	}
	Eigen::VectorXi& getMatching();
	Eigen::VectorXi& getSingularity();

private:
	typedef std::complex<double> COMPLEX;
	Mesh* mesh;
	//std::vector<OpenMesh::Vec3d> faceBase;
	//std::vector<int> constraintId;
	//std::vector<OpenMesh::Vec3d> constraintVector;
	//std::vector<OpenMesh::Vec3d> crossfield;

	std::vector<COMPLEX> ec4;
	Eigen::Matrix3Xd faceBase;
	Eigen::VectorXi constraintId;
	Eigen::Matrix3Xd constraintVector;
	Eigen::Matrix3Xd crossfield;
	Eigen::Matrix3Xd position;

	Eigen::VectorXi matching;
	Eigen::VectorXi singularity;
};


