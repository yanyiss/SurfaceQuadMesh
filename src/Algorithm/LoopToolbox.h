#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "crossField.h"
//#include "../Dependency/CPS/CPS_AABBTree.h"
#include <queue>
namespace LoopGen
{
	class LocalParametrization
	{
	public:
		LocalParametrization(){};
		~LocalParametrization(){};
	public:

	private:

	};


	class LoopGen
	{
	public:
		LoopGen(Mesh& mesh_) :mesh(&mesh_) { /*InitializeAABBTREE();*/ };
		~LoopGen(){};
	public:


	//private:
		Mesh* mesh;
		//ClosestPointSearch::AABBTree* aabbtree;
		Eigen::Matrix4Xd weight;
		crossField* cf = nullptr;

		struct PointOnHalfedge{
			HalfedgeHandle h;
			double c;
			PointOnHalfedge(){}
			PointOnHalfedge(HalfedgeHandle h_, double c_) :h(h_), c(c_) {}
		};
		typedef std::vector<PointOnHalfedge> PlaneLoop;
		struct InfoOnVertex {
			VertexHandle v;
			std::map<InfoOnVertex*, int> mark;
			std::vector<VertexHandle> loop;
			PlaneLoop pl;
			double plane[4];
			double energy;
			bool operator>(const InfoOnVertex& right) const
			{
				return energy > right.energy;
			}
		};
		std::vector<double> eov;
		std::vector<double> simi_e;
		std::vector<InfoOnVertex> InfoOnMesh;
		std::priority_queue<InfoOnVertex, std::vector<InfoOnVertex>, std::greater<InfoOnVertex>> pq;
		//Eigen::Matrix3Xd loop0, loop1;
		//Eigen::Matrix3Xd fragment0, fragment1;

		//void InitializeAABBTREE();
		void InitializeField();
		void InitializeGraphWeight();
		void InitializePQ();
		void ConstructSubMesh(Mesh &submesh);

		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle> &loop, int shift = 0);
		double RefineLoop(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift);
		double ComputeAdjVertexSimilarity(InfoOnVertex& iov0, InfoOnVertex& iov1);

		void GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3]);
		double EvaluateSimilarity(Eigen::Matrix3Xd &loop0, Eigen::Matrix3Xd &loop1, double u, int begin_seg);
	};

	void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
	double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);
}