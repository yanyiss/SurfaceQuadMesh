#pragma once
#include <vector>
#include <deque>
#include <map>
#include <queue>

#include "..\MeshViewer\MeshDefinition.h"

#define YYSS_INFINITE 1.0e12
#define YYSS_FAIRLY_SMALL 1.0e-3

namespace LoopGen
{
	struct PointOnHalfedge {
		HalfedgeHandle h;
		double c = 0.0;
		PointOnHalfedge() {}
		PointOnHalfedge(HalfedgeHandle h_, double c_) :h(h_), c(c_) {}
	};
	typedef std::vector<PointOnHalfedge> PlaneLoop;

	struct InfoOnVertex {
		int id;
		std::map<int, int> mark;//这里的int记录了平行loop的朝向关系，0表示同向，1表示反向
		//std::map<VertexHandle, InfoOnVertex*> adj;
		PlaneLoop pl;
		//double plane[4] = { 0, 0, 0, 0 };
		double energy = 1.0e12;
		/*bool operator>(const InfoOnVertex& right) const
		{
			return energy > right.energy;
		}*/
	};

	struct info_pair
	{
		InfoOnVertex* iov = nullptr;
		double energy = 0;
		info_pair() {}
		info_pair(InfoOnVertex* iov_, double e)
		{
			iov = iov_; energy = e;
		}
		bool operator>(const info_pair& right) const
		{
			return energy > right.energy;
		}
	};
	typedef std::priority_queue<info_pair, std::vector<info_pair>, std::greater<info_pair>> info_pair_pq;

	struct cylinder
	{
		int id;
		std::vector<VertexHandle> vertices;
		std::vector<HalfedgeHandle> bounds;
		std::vector<FaceHandle> faces;
		std::deque<bool> vertice_flag;
		std::deque<bool> bound_flag;
		std::deque<bool> face_flag;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];
	};
	struct cylinder_set
	{
		std::vector<cylinder> cylinders;
		std::deque<bool> set_vertice_flag;
		std::vector<std::vector<int>> vertex_cylinder_map;
		void push_back(cylinder& cy);

	};

}