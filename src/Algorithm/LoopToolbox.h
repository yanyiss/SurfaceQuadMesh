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

		struct EnergyOnVertex {
			VertexHandle v;
			double energy;
			bool operator>(const EnergyOnVertex& right) const
			{
				return energy > right.energy;
			}
		};
		std::priority_queue<EnergyOnVertex, std::vector<EnergyOnVertex>, std::greater<EnergyOnVertex>> pq;

		void InitializeField();
		void InitializeGraphWeight();
		void InitializePQ();
		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<int> &loop, int shift = 0);
	};
}