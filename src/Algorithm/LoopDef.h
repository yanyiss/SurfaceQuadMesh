#pragma once
#include <vector>
#include <deque>
#include <map>
#include <queue>

#include "..\MeshViewer\MeshDefinition.h"
#include "..\Toolbox\dprint.h"

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
		std::map<int, int> mark;//�����int��¼��ƽ��loop�ĳ����ϵ��0��ʾͬ��1��ʾ����
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
		int id = -1;
		Mesh* mesh = nullptr;
		std::vector<VertexHandle> vertices;
		std::vector<std::vector<HalfedgeHandle>> bounds;
		std::vector<FaceHandle> faces;
		std::deque<bool> vertice_flag;
		std::deque<bool> info_on_region;
		std::deque<bool> bound_flag;
		std::deque<bool> face_flag;
		std::deque<bool> cut_v_flag;
		std::deque<bool> cut_f_flag;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

		cylinder() {}
		cylinder(cylinder&& rhs);
		void SetBound();
		OpenMesh::Vec3d GetUGrad(FaceHandle fh);
		OpenMesh::Vec3d GetVGrad(FaceHandle fh);
	};
	struct cylinder_set
	{
		std::vector<cylinder> cylinders;
		std::deque<bool> set_vertice_flag;
		std::vector<std::vector<int>> iov_cylinder_map;
		//void push_back(cylinder& cy);

	};

}