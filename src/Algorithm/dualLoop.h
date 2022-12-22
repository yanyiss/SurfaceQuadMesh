#pragma once
#include <queue>
#include "crossField.h"
namespace QuadLayout
{
#if 0
	using namespace Eigen;
	using std::vector;


	class dualLoop
	{
	public:
		dualLoop(Mesh *m) :mesh(m) {
			//cf = new crossField(mesh);
			cf->runPolynomial();
			/*initDijkstraMark();
			computeGraphWeight(900);*/
		};
		dualLoop(const dualLoop& dl) = delete;
		~dualLoop()
		{
			/*if (cf) { delete cf; cf = nullptr; };
			if (visited) { delete visited; visited = nullptr; }
			if (distance) { delete distance; distance = nullptr; }*/
		}

	private:
#pragma region Dijkstra anisotropic propagation
		//void initDijkstraMark()
		//{
		//	if (visited) delete visited;
		//	if (distance) delete distance;
		//	if (prev) delete prev;
		//	visited = new bool[mesh->n_vertices()];
		//	distance = new double[mesh->n_vertices()];
		//	prev = new OpenMesh::SmartHalfedgeHandle[mesh->n_vertices()];
		//}
		//void recoverDijkstraMark()
		//{
		//	for (int i = 0; i < mesh->n_vertices(); ++i)
		//	{
		//		visited[i] = false;
		//		distance[i] = DBL_MAX;
		//		//prev[i] = -1;
		//	}
		//}

		Matrix4Xd weight;
		void computeFieldAlignedWeight();
		//void computeGraphWeight(double alpha);
	public:
		bool DijkstraLoop(int vid, int cfid, std::vector<int> &loop, double length);//选定顶点和初始方向 cfid = 1, 2，输出路径和路径长度

		//bool* visited = nullptr;
		//double* distance = nullptr;
		//OpenMesh::SmartHalfedgeHandle* prev = nullptr;
		//MatrixXd weight;
		std::vector<int> halfedgeAlignId;

		struct VertexPQ
		{
			int id;
			//bool visited;
			double dist;
			int count;
			VertexPQ(){}
			VertexPQ(int id_, double dist_, int count_) :id(id_), dist(dist_), count(count_) {}
			bool operator>(const VertexPQ &x) const { return dist > x.dist; }
		};
#pragma endregion

	private:
		crossField* cf = nullptr;
		Mesh *mesh;
	};
#endif
}