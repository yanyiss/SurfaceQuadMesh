#pragma once
#include <vector>
#include <map>
#include <queue>

#include "..\MeshViewer\MeshDefinition.h"
#include "..\Toolbox\dprint.h"
#include "covering_space.h"

#define USE_NEW_SIMILARITY_ENERGY 1

namespace LoopGen
{
	struct PointOnHalfedgeLayer {
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

	struct LayerNode
	{
		int id;//顶点编号
		int count;//加入队列的次数
		double dist;//与源点的距离
		LayerNode() {}
		LayerNode(int id_, double dist_, int count_) :id(id_), dist(dist_), count(count_) {}
		bool operator>(const LayerNode& x) const { return dist > x.dist; }
	};
	typedef std::priority_queue<LayerNode, std::vector<LayerNode>, std::greater<LayerNode>> layernode_pq;

	enum DIRECTION {
		Forward,
		Reverse
	};

	struct InfoOnVertex
	{
		int id;
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

	struct spread_info
	{
		M4* m4;

		std::vector<VertexLayer*> new_vertex;
		std::vector<FaceLayer*> new_face;

		BoolVector new_f_flag;
		BoolVector new_v_flag;
		std::vector<int> grow_dir;

		std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;

		Eigen::Matrix3Xd normal_sampling;
		bool has_ns = false;
	};

	struct region
	{
		//最基础信息
		int id = -1;
		std::vector<VertexLayer*> vertices;
		BoolVector vertice_flag;
		//次基础信息
		std::vector<FaceLayer*> faces;
		BoolVector face_flag;
		std::vector<std::vector<HalfedgeLayer*>> bounds;
		std::vector<std::pair<int, int>> handle_to_layer;

		region() {}
		~region() {}
		void set_face(M4 &m4);
		virtual void set_bound() {}
	};

	struct cylinder : region
	{
		std::vector<VertexLayer*> cut;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

		cylinder() {}
		cylinder(cylinder &&cy);
		cylinder& operator=(cylinder&& cy);
		
		inline double GetRegularU(int vlid) { return uv[0](vidmap[vlid]) - std::floor(uv[0](vidmap[vlid])); }
		inline double GetU(int vlid) { return uv[0](vidmap[vlid]); }
		inline double GetV(int vlid) { return uv[1](vidmap[vlid]); }
		void set_bound() override;
		void parametrize(M4 &m4, const Eigen::Matrix3Xd& normal);
	};

	/*struct cylinder_set
	{
		std::vector<cylinder> cylinders;
		std::vector<PlaneLoop> all_path;
		BoolVector bound_edge_flag;
		std::vector<std::pair<int, int>> vertex_bound_index;
		BoolVector has_vertex;
		BoolVector has_face;
		std::vector<BoolVector> tangential_intersection;

		void ProcessOverlap(M4 &m4);
	};
*/
	struct disk : region
	{
		std::pair<int, int> from_bound;
		std::pair<int, int> to_bound;
		BoolVector constraint_f_flag;

		disk(){}
		disk(disk &&dk);
		disk& operator=(disk&& dk);
		void set_bound() override;
	};

	/*struct disk_set
	{
		std::vector<disk> disks;

		void ProcessOverlap(M4 &m4);
	};*/

	struct region_set
	{
		std::vector<region*> regions;
		std::vector<PlaneLoop> bridge;
		BoolVector bound_halfedgelayer_flag;
		BoolVector has_verticelayer;
		BoolVector has_facelayer;
		std::vector<std::pair<int, int>> verticelayer_bound_index;

		int disk_mark = -1;
		region_set() {}
		~region_set() { for (auto rg : regions) delete rg; }
		void ProcessOverlap(M4 &m4, int type);
		void clearRegion() { for (auto rg : regions) if (rg) { delete rg; rg = nullptr; } }
	};
	extern region_set rset;

#define cyPtr(rg) static_cast<cylinder*>(rg)
#define dkPtr(rg) static_cast<disk*>(rg)
}