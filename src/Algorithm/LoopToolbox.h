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

		// ������һ�����������ȶ���,ÿ������һ������
		void InitializeField();
		void InitializeGraphWeight();
		void InitializePQ();
		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<int>& loop, int shift = 0);

#pragma region 

		Eigen::Matrix3Xd edge_fd;  //�ĸ�����ķ��򳡣���ֵ�����ϣ�������edge.idx()*4 + current_direction_index
		std::vector<int> face_M4layer_index_begin;// ��face crossfield ��һ������ ���ڵ�M4��
		int current_direction_index = 3;
		struct LNode {
			VertexHandle v;
			int v_id;
			double distance;
			OpenMesh::SmartHalfedgeHandle from_halfedge;
			bool have_fixed = false; // �����ж�k-ring��dijikstra�Ƿ�ȷ�����·��
			int deepth = 0;			//	��������k-ring�Ĵ���ȵĹ�ȱ���
			bool operator>(const LNode& right) const
			{
				return distance > right.distance;
			}
		};
		int k_default = 5;
		std::vector<std::vector<LNode>> Local_node;
		void cal_M4_edge_Fd();// ����ÿ���ߵķ��򳡲�ֵ��Fd��Fc��Fl
		void cal_Weighted_LocalGraph(int current_direction_index); // ������ָ�������ϵ�Ȩ��ͼ
		std::vector<LNode> cal_K_ring_node(VertexHandle v);
		void calM4Layers(VertexHandle root);
		std::vector<int> get_face_M4layer_index_begin();
		int find_node_in_vector(VertexHandle v, std::vector<LNode>& n_vector);// ����node��vector���λ�ã�û�з���-1
		void cal_local_node();
		void cal_K_ring_dijkstra(VertexHandle root, std::vector<LNode>& n_vector);// ����k-ring �ĵ��״̬������롢·��
		double cal_node_distance(VertexHandle root, VertexHandle target,std::vector<LNode>& n_vector );
		std::vector<OpenMesh::SmartHalfedgeHandle> cal_k_ring_path(VertexHandle root, VertexHandle target, std::vector<LoopGen::LNode>& node_vector);// ���k-ring�������·����ֻ����cal_K_ring_dijkstra��ʹ��

		// �жϵ�ǰglobalnode�Ƿ�̶�
		struct GNode {
			int v_id;				// GNode ��v_id
			int last_v_id;			// GNode ����һ����v_id
			double distance;		// root �� GNode �ľ���
			bool operator>(const GNode& right) const
			{
				return distance > right.distance;
			}
		};// ��ȫ�ֵ�node��ֻ����v_id��last_v_id��distance
		std::vector<GNode> Global_node; // ��ȫ��dijkstra��Ľ��
		void cal_global_dijkstra(VertexHandle root, int current_direction_index);
		std::vector<GNode> get_global_node();
		std::vector<OpenMesh::SmartHalfedgeHandle> cal_global_path(VertexHandle root, VertexHandle target);

#pragma endregion

		

	};
}