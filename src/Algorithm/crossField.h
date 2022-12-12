#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include <complex>
#include <Eigen\Sparse>

#include "..\Toolbox\mathFunctions.h"
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
	crossField(Mesh *m) {
		init(m);
		setCurvatureConstraint();
		runPolynomial();
	};
	crossField(std::string& field_file);
	crossField(const crossField& cf) = delete;
	~crossField() {};

private:
	void calculateMeshFaceBase();
public:
	void init(Mesh *m) 
	{
		mesh = m;
		position.resize(3, m->n_vertices());
		for (auto& tv : m->vertices())
		{
			auto& p = mesh->point(tv);
			position.col(tv.idx()) << p[0], p[1], p[2];
		}
		calculateMeshFaceBase();
	}
	void runPolynomial();
	void runLocalOpt(std::vector<FaceHandle>& opt_face, std::deque<bool> &grow_dir, std::deque<bool>& opt_flag, std::deque<bool>& constraint_flag);
	//void runIteration(std::vector<OpenMesh::Vec3d>& crossfield);

	void setCurvatureConstraint();
	void setOuterConstraint(std::deque<bool> &cons_flag, Eigen::Matrix3Xd &cons_direction);
	void setMatching();
	void setSingularity();
	void setNormal();

	void write_field(std::string &field_file);

	const std::vector<int>& getConstraintFace() { return constraintId; }
	Eigen::Matrix3Xd& getCrossField() { return crossfield; }
	const Eigen::Matrix3Xd& getPosition() { return position; }
	const Eigen::Matrix3Xd& getNormal() { return normal; }
	const Eigen::Matrix3Xd& getFaceBase() { return faceBase; }
	std::vector<int>& getMatching() { return matching; }
	const std::vector<int>& getSingularity() { return singularity; }
	

//private:
	typedef std::complex<double> COMPLEX;
	Mesh* mesh = nullptr;
	//std::vector<OpenMesh::Vec3d> faceBase;
	//std::vector<int> constraintId;
	//std::vector<OpenMesh::Vec3d> constraintVector;
	//std::vector<OpenMesh::Vec3d> crossfield;

	std::vector<COMPLEX> ec4;
	Eigen::Matrix3Xd faceBase;
	std::vector<int> constraintId;
	Eigen::Matrix3Xd constraintVector;
	Eigen::Matrix3Xd crossfield;
	Eigen::Matrix3Xd position;
	Eigen::Matrix3Xd normal;

	std::vector<int> matching;
	std::vector<int> singularity;
};


