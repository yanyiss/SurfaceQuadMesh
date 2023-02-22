#pragma once

#include "crossField.h"
#include "LoopDef.h"
#include "covering_space.h"
#include "LocalParametrization.h"
//#include "../Dependency/CPS/CPS_AABBTree.h"
#include <queue>
namespace LoopGen
{
	class LoopGen
	{
	public:
		LoopGen(Mesh& mesh_) :mesh(&mesh_) { /*InitializeAABBTREE();*/initMeshStatusAndNormal(*mesh); };
		~LoopGen() {
			if (cf) { delete cf; cf = nullptr; }
		};
	public:


	//private:
		//模型基础信息
		std::string model_name;
		Mesh* mesh;
		//模型加工信息
		crossField* cf = nullptr;
		//M2 m2;
		M4 m4;
		timeRecorder tr;

		
		//初始化阶段变量
		std::vector<InfoOnVertex> InfoOnMesh;
		PLS pls;

		//基础函数
		void SetModelName(std::string& name) { model_name = name; truncateFileName(model_name); }

		//初始化阶段函数
		//bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle>& loop, int shift = 0);
		bool FieldAligned_PlanarLoop(/*M4 &m4, */VertexLayer* vl, std::vector<VertexLayer*> &loop);
		double RefineLoopByPlanarity(/*M4 &m4, */std::vector<VertexLayer*>& loop, PlaneLoop& planar_loop);
		void GetPositionFromLoop(const std::vector<VertexLayer*>& loop, Eigen::VectorXd xyz[3]);
		double EvaluateSimilarity(Eigen::Matrix3Xd& loop0, Eigen::Matrix3Xd& loop1, double u, int begin_seg);
		void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
		double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);

		void InitializeField();
		//void InitializeGraphWeight(double alpha = 900);
		void InitializePQ();

		//用于数据交换与绘制的变量
		//std::vector<VertexHandle> region_vertex;
		//std::vector<VertexHandle> sub_vertex; std::vector<FaceHandle> sub_face;
		BoolVector old_face_flag;
		BoolVector old_vert_flag;
		BoolVector new_face_flag;
		BoolVector new_vert_flag;
		Eigen::VectorXd uv_para[2];
		std::vector<int> vertexidmap;
		Eigen::Matrix3Xd xaxis;
		BoolVector cut_vertex_flag;
		std::vector<int> growDIR;
		BoolVector v_cache_flag;
		std::vector<VertexLayer*> v_cache;

		std::deque<bool> bound_edge_flag;
		std::vector<double> similarity_energy;
		std::vector<PlaneLoop> all_plane_loop;
		//std::vector<Vec3d> u0point5;
		std::vector<VertexLayer*> seed_vertex;
		std::vector<VertexLayer*> cut_vertex;

		//算法重要参数
		double energy_threshold = 0.2;
		int extend_layer = 3;

		//优化阶段变量
		cylinder_set cset;
		//优化阶段函数
		void AssembleSimilarityAngle(VertexLayer* vl, Eigen::VectorXd& sa, LocalParametrization& lp, int loop_fragment_num);
		bool RefineLoopByParametrization(VertexLayer* vl, LocalParametrization& lp, BoolVector& visited_v, BoolVector& visited_f);
		void ResetLocalField(LocalParametrization &lp, std::vector<FaceLayer*>& opt_face, BoolVector& opt_flag, BoolVector& constraint_flag);
		double LoopLenGrad(std::vector<VertexLayer*> &vertex_set, LocalParametrization &lp, BoolVector &vertex_flag, int growDir);
		void AssembleIOVLoopEnergy(info_pair_pq &pq);
		bool CheckTopology(std::vector<VertexLayer*>& vertex_set, BoolVector &set_flag, std::vector<int> &grow_dir);

		void ConstructInitialRegion(InfoOnVertex* iov, LocalParametrization &lp);
		bool SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2]);
		bool ConstructRegionCut(VertexLayer* vl, BoolVector& visited, std::vector<VertexLayer*>& cut);
		void ProcessOverlap(std::vector<std::vector<int>> &region_index);
		void OptimizeLoop();

		//Loop迭代阶段
		void ReLoop();
		void LocalOpt();

		//void SetUParaLine(InfoOnVertex& iov, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f);
	};

}