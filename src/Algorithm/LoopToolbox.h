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
		LoopGen(Mesh& mesh_) :mesh(&mesh_) { };
		~LoopGen() {};
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

		// 定义了一个能量的优先队列,每个点上一个能量
		void InitializeField();
		void InitializeGraphWeight();
		void InitializePQ();
		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<int>& loop, int shift = 0);

#pragma region 

		Eigen::Matrix3Xd edge_fd;  //四个方向的方向场，插值到边上，索引是edge.idx()*4 + current_direction_index
		std::vector<int> face_M4layer_index_begin;// 存face crossfield 第一个方向 所在的M4层
		int current_direction_index = 3;
		struct LNode {
			VertexHandle v;
			int v_id;
			double distance;
			OpenMesh::SmartHalfedgeHandle from_halfedge;
			bool have_fixed = false; // 用来判断k-ring的dijikstra是否确定最短路径
			int deepth = 0;			//	用来进行k-ring的带深度的广度遍历
			bool operator>(const LNode& right) const
			{
				return distance > right.distance;
			}
		};
		int k_default = 5;
		std::vector<std::vector<LNode>> Local_node;
		void cal_M4_edge_Fd();// 计算每个边的方向场插值，Fd、Fc、Fl
		void cal_Weighted_LocalGraph(int current_direction_index); // 计算在指定方向场上的权重图
		std::vector<LNode> cal_K_ring_node(VertexHandle v);
		void calM4Layers(VertexHandle root);
		std::vector<int> get_face_M4layer_index_begin();
		int find_node_in_vector(VertexHandle v, std::vector<LNode>& n_vector);// 返回node在vector里的位置，没有返回-1
		void cal_local_node();
		void cal_K_ring_dijkstra(VertexHandle root, std::vector<LNode>& n_vector);// 更新k-ring 的点的状态，存距离、路径
		double cal_node_distance(VertexHandle root, VertexHandle target,std::vector<LNode>& n_vector );
		std::vector<OpenMesh::SmartHalfedgeHandle> cal_k_ring_path(VertexHandle root, VertexHandle target, std::vector<LoopGen::LNode>& node_vector);// 获得k-ring里两点间路径，只能在cal_K_ring_dijkstra后使用

		// 判断当前globalnode是否固定
		struct GNode {
			int v_id;				// GNode 的v_id
			int last_v_id;			// GNode 的上一个的v_id
			double distance;		// root 到 GNode 的距离
			bool operator>(const GNode& right) const
			{
				return distance > right.distance;
			}
		};// 存全局的node，只包含v_id、last_v_id、distance
		std::vector<GNode> Global_node; // 存全局dijkstra后的结果
		void cal_global_dijkstra(VertexHandle root, int current_direction_index);
		std::vector<GNode> get_global_node();
		std::vector<OpenMesh::SmartHalfedgeHandle> cal_global_path(VertexHandle root, VertexHandle target);

#pragma endregion

		

	};
}