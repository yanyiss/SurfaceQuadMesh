#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "crossField.h"
//#include "../Dependency/CPS/CPS_AABBTree.h"
#include <queue>
namespace LoopGen
{
	struct PointOnHalfedge {
		HalfedgeHandle h;
		double c;
		PointOnHalfedge() {}
		PointOnHalfedge(HalfedgeHandle h_, double c_) :h(h_), c(c_) {}
	};
	typedef std::vector<PointOnHalfedge> PlaneLoop;
	struct InfoOnVertex {
		VertexHandle v;
		std::map<InfoOnVertex*, int> mark;//这里的int记录了平行loop的朝向关系，0表示同向，1表示反向
		std::vector<VertexHandle> loop;
		PlaneLoop pl;
		double plane[4] = { 0, 0, 0, 0 };
		double energy = 1.0e12;
		/*bool operator>(const InfoOnVertex& right) const
		{
			return energy > right.energy;
		}*/
	};
	struct info_pair
	{
		InfoOnVertex* iov;
		double energy;
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
    
    class LocalParametrization
	{
	public:
		LocalParametrization(Mesh &mesh_, crossField &cf_, VertexHandle v/*, int shift*/) : mesh(&mesh_), cf(&cf_) 
		{ 
			int nf = mesh->n_faces();
			int nv = mesh->n_vertices();
			region_f_flag.resize(nf, false);
			region_v_flag.resize(nv, false);
			region_vertex.push_back(v);
			region_v_flag[v.idx()] = true;
			uv[0].resize(1); uv[0].setZero();
			uv[1].resize(1); uv[1].setZero();
			all_pl.resize(nv);
			vidmap.resize(nv); vidmap[v.idx()] = 0;
			x_axis.resize(3, nf);
			y_axis.resize(3, nf);
			grow_dir.resize(nv, -1);
		};
		~LocalParametrization(){};
	public:
		inline std::vector<VertexHandle>& GetNewVertex() { return new_vertex; }
		inline std::vector<FaceHandle>& GetNewFace() { return new_face; }
		inline std::vector<VertexHandle>& GetRegionVertex() { return region_vertex; }
		inline std::vector<FaceHandle>& GetRegionFace() { return region_face; }
		inline std::vector<VertexHandle>& GetCut() { return cut; }

		inline std::deque<bool>& GetNewFFlag() { return new_f_flag; }
		inline std::deque<bool>& GetNewVFlag() { return new_v_flag; }
		inline std::deque<bool>& GetRegionFFlag() { return region_f_flag; }
		inline std::deque<bool>& GetRegionVFlag() { return region_v_flag; }
		inline std::deque<bool>& GetCutV_Flag() { return cutv_flag; }
		inline std::deque<bool>& GetCutF_Flag() { return cutf_flag; }
		inline std::vector<int>& GetGrowDir() { return grow_dir; }

		inline std::vector<PlaneLoop>& GetAllPL() { return all_pl; }
		inline Eigen::Matrix3Xd& GetXAxis() { return x_axis; }
		inline Eigen::Matrix3Xd& GetYAxis() { return y_axis; }
		inline std::vector<int>& GetVidMap() { return vidmap; }
		inline double GetRegularU(int vid) { return uv[0](vidmap[vid]) - std::floor(uv[0](vidmap[vid])); }
		inline double GetU(int vid) { return uv[0](vidmap[vid]); }
		inline double GetV(int vid) { return uv[1](vidmap[vid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }
		inline Eigen::VectorXd& GetNormalSimilarityAngle() { return normal_similarity_angle; }


	//private:
		Mesh* mesh;
		crossField* cf;

		std::vector<VertexHandle> new_vertex;
		std::vector<FaceHandle> new_face;
		std::vector<VertexHandle> region_vertex;
		std::vector<FaceHandle> region_face;
		std::vector<VertexHandle> cut;

		std::deque<bool> new_f_flag;
		std::deque<bool> new_v_flag;
		std::deque<bool> region_f_flag;
		std::deque<bool> region_v_flag;
		std::deque<bool> cutv_flag;
		std::deque<bool> cutf_flag;
		std::vector<int> grow_dir;

		std::vector<PlaneLoop> all_pl;
		Eigen::Matrix3Xd x_axis;
		Eigen::Matrix3Xd y_axis;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];
		Eigen::VectorXd normal_similarity_angle;
		bool has_nsa = false;
		InfoOnVertex* iov;

		void run();
	};


	class LoopGen
	{
	public:
		LoopGen(Mesh& mesh_) :mesh(&mesh_) { /*InitializeAABBTREE();*/ };
		~LoopGen() {
			if (cf) { delete cf; cf = nullptr; }
		};
	public:


	//private:
		//模型基础信息
		std::string model_name;
		Mesh* mesh;
		//模型加工信息
		Eigen::Matrix4Xd weight;
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
		void InitializeGraphWeight(double alpha = 900);
		void InitializePQ();

		//用于数据交换与绘制的变量
		std::vector<VertexHandle> region_vertex;
		std::vector<VertexHandle> sub_vertex; std::vector<FaceHandle> sub_face;
		Eigen::VectorXd uv_para[2];
		std::vector<double> similarity_energy;
		std::vector<int> idmap;
		std::vector<PlaneLoop> all_plane_loop;
		//std::vector<Vec3d> u0point5;
		std::vector<VertexHandle> seed_vertex;
		std::vector<VertexHandle> cut_vertex;
		std::vector<int> growDIR;

		//算法重要参数
		double energy_threshold = 0.18;
		int extend_layer = 3;

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
		void OptimizeLoop();

		void SetUParaLine(InfoOnVertex& iov, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f);
	};

}