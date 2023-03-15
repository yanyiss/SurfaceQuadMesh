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
		inline OpenMesh::Vec3d point(MX &mx)
		{
			if (hl)
				return c * mx.mesh->point(mx.verticelayers[hl->from].v) +
				(1 - c) * mx.mesh->point(mx.verticelayers[hl->to].v);
			else
				return mx.mesh->point(mx.mesh->vertex_handle(c));
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

	struct vl_pair
	{
		VertexLayer* vl = nullptr;
		double data = 0;
		vl_pair() {}
		vl_pair(VertexLayer* vl_, double e) : vl(vl_), data(e) { }
		bool operator>(const vl_pair& right) const
		{
			return data > right.data;
		}
	};
	typedef std::priority_queue<vl_pair, std::vector<vl_pair>, std::greater<vl_pair>> vl_pair_pq;

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
		std::vector<std::pair<int, int>> handle_to_layer;

		cylinder::cylinder() {}
		cylinder::~cylinder() {}
		void set_bound();
		void parametrize(M4 &m4, const Eigen::Matrix3Xd& normal);
	};

	struct cylinder_set
	{
		std::vector<cylinder> cylinders;
		std::vector<PlaneLoop> all_path;
		BoolVector bound_edge_flag;
		std::vector<std::pair<int, int>> vertex_bound_index;
		BoolVector has_vertex;
		BoolVector has_face;
		std::vector<BoolVector> tangential_intersection;
	};

	struct disk
	{
		std::vector<VertexHandle> new_vertex;
		std::vector<FaceHandle> new_face;
		std::vector<VertexHandle> region_vertex;
		std::vector<FaceHandle> region_face;

		BoolVector new_f_flag;
		BoolVector new_v_flag;
		BoolVector region_f_flag;
		BoolVector region_v_flag;
		BoolVector constraint_f_flag;
		std::vector<int> grow_dir;

		std::vector<HalfedgeHandle> bounds;
		PlaneLoop initial_path;
		//std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;

		std::pair<int, int> from_bound;
		std::pair<int, int> to_bound;

		Eigen::VectorXd normal_similarity_angle;
		double length[3] = { YYSS_INFINITE, -YYSS_INFINITE, -YYSS_INFINITE };
		bool has_nsa = false;
	};

	struct disk_set
	{
		std::vector<disk> disks;
	};

	/*struct disk
	{
		std::vector<VertexHandle> new_vertex;
		std::vector<FaceHandle> new_face;
		std::vector<VertexHandle> region_vertex;
		std::vector<FaceHandle> region_face;

		BoolVector new_f_flag;
		BoolVector new_v_flag;
		BoolVector region_f_flag;
		BoolVector region_v_flag;
		BoolVector constraint_f_flag;
		std::vector<int> grow_dir;

		std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

		std::pair<int, int> from_bound;
		std::pair<int, int> to_bound;

		Eigen::VectorXd normal_similarity_angle;
		double umx[2] = { 0, 0 };
		double length[3] = { YYSS_INFINITE, -YYSS_INFINITE, -YYSS_INFINITE };
		bool has_nsa = false;

		inline double GetU(int vid) { return uv[0](vidmap[vid]); }
		inline double GetV(int vid) { return uv[1](vidmap[vid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }
	};*/
}