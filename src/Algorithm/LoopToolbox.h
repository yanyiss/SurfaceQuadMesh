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
		LocalParametrization(Mesh &mesh_, crossField &cf_) : mesh(&mesh_), cf(&cf_) { };
		~LocalParametrization(){};
	public:
		inline std::vector<VertexHandle>& GetVertex() { return vertex; }
		inline std::vector<FaceHandle>& GetFace() { return face; }
		inline std::vector<VertexHandle>& GetCut() { return cut; }
		inline std::deque<bool>& GetFFlag() { return f_flag; }
		inline std::deque<bool>& GetVFalg() { return v_flag; }
		inline std::deque<bool>& GetCutV_Flag() { return cutv_flag; }
		inline std::deque<bool>& GetCutF_Flag() { return cutf_flag; }
		//int GetVidMap(int vid) { return vidmap[vid]; }
		inline double GetU(int vid) { return uv[0](vidmap[vid]); }
		inline double GetV(int vid) { return uv[1](vidmap[vid]); }

	//private:
		Mesh* mesh;
		std::vector<VertexHandle> vertex;
		std::vector<FaceHandle> face;
		std::vector<VertexHandle> cut;
		std::deque<bool> f_flag;
		std::deque<bool> v_flag;
		std::deque<bool> cutv_flag;
		std::deque<bool> cutf_flag;
		crossField* cf;
		std::vector<int> vidmap;
		Eigen::VectorXd uv[2];

		void run(VertexHandle v, int shift);
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
		std::vector<std::vector<InfoOnVertex*>> advancing_front[2];
		Eigen::VectorXd uv_para[2];
		std::vector<double> similarity_energy;
		std::vector<InfoOnVertex> InfoOnMesh;
		//std::vector<std::vector<InfoOnVertex*>> advancing_front[2];
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
		void ConstructSubRegion(InfoOnVertex* iov, std::vector<std::vector<InfoOnVertex*>> advancing_front[2]);
		bool SpreadSubRegion(std::vector<std::vector<InfoOnVertex*>> advancing_front[2], LocalParametrization& lp);
		void ConstructRegionCut(VertexHandle v, int shift, std::deque<bool>& visited, std::vector<VertexHandle> &cut);
		void OptimizeLoop();
		bool IsGood(std::vector<InfoOnVertex*>& af, std::deque<bool> &update_mark, LocalParametrization &lp, double threshold = 2.0);

		bool FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle> &loop, int shift = 0);
		double RefineLoopByPlanarity(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift);
		bool RefineLoopByParametrization(InfoOnVertex &iov, LocalParametrization& lp);
		void GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3]);
		double EvaluateSimilarity(Eigen::Matrix3Xd &loop0, Eigen::Matrix3Xd &loop1, double u, int begin_seg);
	};

	void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4]);
	double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4]);
}