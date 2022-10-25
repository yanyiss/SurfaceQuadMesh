#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "crossField.h"
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
		LoopGen(Mesh &mesh_):mesh(&mesh_) { };
		~LoopGen(){};
	public:


	//private:
		Mesh* mesh;
		Eigen::Matrix4Xd weight;
		crossField* cf = nullptr;

		struct InfoOnVertex {
			VertexHandle v;
			std::vector<VertexHandle> loop[2];
			double plane[2][4];
			double energy;
			bool operator>(const InfoOnVertex& right) const
			{
				return energy > right.energy;
			}
		};
		std::vector<double> eov;
		std::priority_queue<InfoOnVertex, std::vector<InfoOnVertex>, std::greater<InfoOnVertex>> pq;

		void InitializeField();
		void InitializeGraphWeight();
		void InitializePQ();
		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle> &loop, int shift = 0);

		void GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3]);

	};

	void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
	double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);
}