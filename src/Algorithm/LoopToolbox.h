#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "crossField.h"
//#include "../Dependency/CPS/CPS_AABBTree.h"
#include <queue>
namespace LoopGen
{
	class LocalParametrization
	{
	public:
		LocalParametrization(Mesh &mesh_, crossField &cf_, VertexHandle v, int shift) : mesh(&mesh_), cf(&cf_) 
		{ 
			region_f_flag.resize(mesh->n_faces(), false);
			region_v_flag.resize(mesh->n_vertices(), false);
			region_vertex.push_back(v);
			region_v_flag[v.idx()] = true;
			uv[0].resize(1); uv[0].setZero();
			uv[1].resize(1); uv[1].setZero();
			vidmap.resize(mesh->n_vertices()); vidmap[v.idx()] = 0;
			grow_dir.resize(mesh->n_vertices(), false);
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
		inline std::deque<bool>& GetGrowDir() { return grow_dir; }
		inline std::vector<int>& GetVidMap() { return vidmap; }
		inline double GetRegularU(int vid) { return uv[0](vidmap[vid]) - std::floor(uv[0](vidmap[vid])); }
		inline double GetU(int vid) { return uv[0](vidmap[vid]); }
		inline double GetV(int vid) { return uv[1](vidmap[vid]); }
		inline Eigen::VectorXd& GetU() { return uv[0]; }
		inline Eigen::VectorXd& GetV() { return uv[1]; }

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
		std::deque<bool> grow_dir;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

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
		Mesh* mesh;
		//ClosestPointSearch::AABBTree* aabbtree;
		Eigen::Matrix4Xd weight;
		crossField* cf = nullptr;

		struct PointOnHalfedge{
			HalfedgeHandle h;
			double c;
			PointOnHalfedge(){}
			PointOnHalfedge(HalfedgeHandle h_, double c_) :h(h_), c(c_) {}
		};
		typedef std::vector<PointOnHalfedge> PlaneLoop;
		struct InfoOnVertex {
			VertexHandle v;
			std::map<InfoOnVertex*, int> mark;//这里的int记录了平行loop的朝向关系，0表示同向，1表示反向
			std::vector<VertexHandle> loop;
			PlaneLoop pl;
			double plane[4];
			double energy;
			bool operator>(const InfoOnVertex& right) const
			{
				return energy > right.energy;
			}
		};
		std::vector<double> eov;
		//std::vector<std::vector<InfoOnVertex*>> advancing_front[2];
		std::vector<VertexHandle> region_vertex;
		std::vector<int> idmap;
		Eigen::VectorXd uv_para[2];
		std::vector<double> similarity_energy;
		std::vector<InfoOnVertex> InfoOnMesh;
		std::priority_queue<InfoOnVertex, std::vector<InfoOnVertex>, std::greater<InfoOnVertex>> pq;
		std::vector<VertexHandle> sub_vertex; std::vector<FaceHandle> sub_face; std::vector<VertexHandle> sub_cut;
		//Eigen::Matrix3Xd loop0, loop1;
		//Eigen::Matrix3Xd fragment0, fragment1;
		std::string model_name;
		void SetModelName(std::string& name) { model_name = name; truncateFileName(model_name); }

		//void InitializeAABBTREE();
		void InitializeField();
		void InitializeGraphWeight(double alpha = 899);
		void InitializePQ();
		void ConstructInitialRegion(InfoOnVertex* iov, LocalParametrization &lp);
		bool SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2]);
		void ResetField(LocalParametrization& lp);
		void ConstructRegionCut(VertexHandle v, int shift, std::deque<bool>& visited, std::vector<VertexHandle> &cut);
		void OptimizeLoop();
		bool IsGood(InfoOnVertex* iov0, InfoOnVertex* iov1, LocalParametrization& lp, double threshold = 2.0);

		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle> &loop, int shift = 0);
		double RefineLoopByPlanarity(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift);
		bool RefineLoopByParametrization(InfoOnVertex &iov, LocalParametrization& lp, std::deque<bool> &visited_v, std::deque<bool> &visited_f);
		void GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3]);
		double EvaluateSimilarity(Eigen::Matrix3Xd &loop0, Eigen::Matrix3Xd &loop1, double u, int begin_seg);
	};

	void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
	double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);
}