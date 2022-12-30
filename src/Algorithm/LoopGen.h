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
		M2 m2;
		M4 m4;
		//模型加工信息
		crossField* cf = nullptr;
		
		//初始化阶段变量
		std::vector<InfoOnVertex> InfoOnMesh;

		//基础函数
		void SetModelName(std::string& name) { model_name = name; truncateFileName(model_name); }

		//初始化阶段函数
		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle>& loop, int shift = 0);
		double RefineLoopByPlanarity(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift);
		void GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3]);
		double EvaluateSimilarity(Eigen::Matrix3Xd& loop0, Eigen::Matrix3Xd& loop1, double u, int begin_seg);
		void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
		double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);

		void InitializeField();
		//void InitializeGraphWeight(double alpha = 900);
		void InitializePQ();

		//用于数据交换与绘制的变量
		//std::vector<VertexHandle> region_vertex;
		//std::vector<VertexHandle> sub_vertex; std::vector<FaceHandle> sub_face;
		std::deque<bool> optimized_face_flag;
		std::deque<bool> optimized_vert_flag;
		std::deque<bool> bound_edge_flag;
		Eigen::VectorXd uv_para[2];
		std::vector<double> similarity_energy;
		std::vector<int> idmap;
		std::vector<PlaneLoop> all_plane_loop;
		//std::vector<Vec3d> u0point5;
		std::vector<VertexHandle> seed_vertex;
		std::vector<VertexHandle> cut_vertex;
		std::vector<int> growDIR;

		//算法重要参数
		double energy_threshold = 0.2;
		int extend_layer = 3;

		//优化阶段变量
		cylinder_set cset;

		//优化阶段函数
		void AssembleSimilarityAngle(VertexHandle v, Eigen::VectorXd& sa, LocalParametrization& lp, int loop_fragment_num);
		bool RefineLoopByParametrization(VertexHandle v, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f);
		void ResetLocalField(LocalParametrization &lp, std::vector<FaceHandle>& opt_face, std::deque<bool>& opt_flag, std::deque<bool>& constraint_flag);
		double LoopLenGrad(std::vector<VertexHandle> &vertex_set, LocalParametrization &lp, std::deque<bool> &vertex_flag, int growDir);
		void AssembleIOVLoopEnergy(info_pair_pq &pq);
		bool CheckTopology(std::vector<VertexHandle>& vertex_set, std::deque<bool> &set_flag, std::vector<int> &grow_dir);

		void ConstructInitialRegion(InfoOnVertex* iov, LocalParametrization &lp);
		bool SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2]);
		bool ConstructRegionCut(InfoOnVertex* iov, std::deque<bool>& visited, std::vector<VertexHandle>& cut);
		void UpdateIOM();
		void ProcessOverlap();
		void LookForCrossingLoop();
		void OptimizeLoop();


		void SetUParaLine(InfoOnVertex& iov, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f);
	};

}