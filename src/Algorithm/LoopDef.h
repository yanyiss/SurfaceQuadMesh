#pragma once
#include <vector>
#include <map>
#include <queue>

#include "..\MeshViewer\MeshDefinition.h"
#include "..\Toolbox\dprint.h"
#include "covering_space.h"


namespace LoopGen
{
	struct PointOnHalfedgeLayer {
		//HalfedgeHandle h;
		HalfedgeLayer* hl;
		double c = 0.0;
		PointOnHalfedgeLayer() {}
		PointOnHalfedgeLayer(HalfedgeLayer* hl_, double c_) :hl(hl_), c(c_) {}
		template <typename MX>
		OpenMesh::Vec3d point(MX &mx)
		{
			return c * mx.mesh->point(mx.verticelayers[hl->from].v)
				+ (1 - c) * mx.mesh->point(mx.verticelayers[hl->to].v);
		}
	};
	typedef std::vector<PointOnHalfedgeLayer> PlaneLoop;
	typedef std::vector<PlaneLoop> PLS;

	//struct InfoOnVertex {
	//	int id;
	//	std::map<int, int> mark;//这里的int记录了平行loop的朝向关系，0表示同向，1表示反向
	//	//std::map<VertexHandle, InfoOnVertex*> adj;
	//	//PlaneLoop pl;

	//	double energy = YYSS_INFINITE;
	//};
	enum DIRECTION {
		Forward,
		Reverse
	};
	struct InfoOnVertex
	{
		int id;
		//std::map<int, int> mark;
		//PlaneLoop pl;
		int plid;
		DIRECTION dir;
		double energy = YYSS_INFINITE;
	};

	struct info_pair
	{
		VertexLayer* vl = nullptr;
		double energy = 0;
		info_pair() {}
		info_pair(VertexLayer* vl_, double e) : vl(vl_), energy(e) { }
		bool operator>(const info_pair& right) const
		{
			return energy > right.energy;
		}
	};
	typedef std::priority_queue<info_pair, std::vector<info_pair>, std::greater<info_pair>> info_pair_pq;

	struct cylinder
	{
		int id = -1;
		std::vector<VertexLayer*> vertices;
		std::vector<FaceLayer*> faces;
		std::vector<VertexLayer*> cut;
		std::vector<std::vector<HalfedgeLayer*>> bounds;
		BoolVector vertice_flag;
		BoolVector face_flag; 

		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

		cylinder::cylinder() {}
		cylinder::~cylinder() {}
		void set_bound();
		void parametrize(M4 &m4, const Eigen::Matrix3Xd& normal);
	};

	/*struct cylinder
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
	};*/
	struct cylinder_set
	{
		std::vector<cylinder> cylinders;
		//void push_back(cylinder& cy);

	};

}