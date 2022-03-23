#pragma once
/*
-- A method for generating quad-dominant meshes inspired from paper: Instant Field-Aligned Meshes
                                                                 and QuadriFlow: A Scalable and Robust Method for Quadrangulation
-- By yanyisheshou at GCL in USTC
-- Email: duweiyou@mail.ustc.edu.cn
-- 2021/09/18
*/

#include <vector>
#include <Eigen\Core>
#include "..\src\MeshViewer\MeshDefinition.h"
//#include "Tooler/flow/flow.hpp"
#include <tbb/tbb.h>

#define DefaultDepth 16
#define USING_SMALLMESH

#define RUN_MESHFIELD_PARALLEL

#ifdef RUN_MESHFIELD_PARALLEL
#define SerialParallelBlock(codeBody) tbb::parallel_for(static_cast<ui>(0), nv, [&](ui id){##codeBody##});
#else
#define SerialParallelBlock(codeBody) for (ui id = 0; id < nv; id++){##codeBody##}
#endif // RUN_MESHFIELD_PARALLEL

#if 0
class MeshField
{
public:
	MeshField() {};
	virtual ~MeshField() {};

protected:

	typedef unsigned int ui;
//#define invalid (ui)-1
	typedef Eigen::Vector3d Vec3d;
	typedef std::vector<Vec3d> vv3d;
	struct Link
	{
		ui id;
		double weight;
		Link() {}
		Link(ui v_id, double w) :id(v_id), weight(w) {}
		Link(ui v_id) :id(v_id), weight(1.0) {}
		inline bool operator<(const Link &link)const { return id < link.id; }
		inline bool operator==(const Link &link) const { return id == link.id; }
	};
	typedef std::vector<Link> vL;
	typedef std::vector<vL> vvL;

	struct edgeWeight
	{
		ui i; ui j;
		double w;
		edgeWeight(ui i_, ui j_, double w_) :i(i_), j(j_), w(w_) {}
		inline bool operator>(const edgeWeight &right) const { return w > right.w; }
		inline bool operator<(const edgeWeight &right) const { return w < right.w; }
	};

	struct hierarchy {
		vv3d v;//顶点坐标
		vv3d n;//法向
		std::vector<double> a;//面积
		std::vector<std::vector<Link>> adj;//邻接关系
		std::vector<ui> tocoarser;//寻找粗糙层指标
		std::vector<std::pair<ui, ui>> tofiner;//寻找精细层指标
		struct hierarchy* coarser_hierarchy;//粗糙层指针
		struct hierarchy* finer_hierarchy;//精细层指针
		std::vector<std::vector<ui>> phase;//网格边染色
		ui hierarchy_depth;//当前层数

		

		vv3d o;//orientation
		vv3d p;//position

		//给定z轴，计算x轴
		inline void calcXAxis(const Vec3d &z, Vec3d &x)
		{
			using Eigen::Vector3d;
			using std::fabs;
			if (fabs(z[0]) > fabs(z[1]))
			{
				x = Vec3d(z[2], 0.0, -z[0]).normalize();
			}
			else
			{
				x = Vec3d(0.0, z[2], -z[1]).normalize();
			}
		}

		//由mesh直接创建depth层的粗糙网格
		hierarchy(const TriMesh &mesh, ui depth);
		~hierarchy();
	private:
		//计算Voronoi面积和邻接关系，一个较为简单的版本
		double calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv);
		//由mesh创建当前层网格
		void init(const TriMesh &mesh);
		//计算并行块
		void graph_color();
		//由当前精细层网格创建下一粗糙层网格
		hierarchy(hierarchy* fh);
	};
	static struct hierarchy *h;

	virtual void randomInit(TriMesh &mesh) = 0;
	virtual void GSOptimize() = 0;
public:
	std::vector<std::vector<ui>> returnPhase() { return h->phase; }
};

#include <OpenMesh\Tools\Subdivider\Uniform\LoopT.hh>
#include <list>
class OrientationField :public MeshField
{
public:
	OrientationField(TriMesh &mesh)
	{
		randomInit(mesh);
	}
	~OrientationField()
	{
		//if (h) { delete h; h = nullptr; }
	};

public:
	void GSOptimize() override;
	vv3d returnOrientation() { return h->o; }

private:
	void loopDivision(TriMesh &mesh, ui divideTimes = 0);
	void randomInit(TriMesh &mesh) override;//初始化hierarchy h & orientation o
};

class PositionField :public MeshField
{
public:
	PositionField(TriMesh &mesh)
	{
		randomInit(mesh);
	}
	~PositionField()
	{
		if (h) { delete h; h = nullptr; }
	}

public:
	void GSOptimize() override;
	vv3d returnPosition() { return h->p; }
	struct VertexValue {
		OV vert;
		Vec3d normal;
		std::set<ui> adj;
		std::list<ui> orientedAdj;
		VertexValue(){}
		VertexValue(OV &v, Vec3d &n, std::set<ui> &a) :vert(v), normal(n), adj(a) {}
	};
	//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
	void extractMesh(PolyMesh &polymesh, std::vector<Vec3d> &pos, std::vector<VertexValue> &polyFrame);
	enum LinkType { SAME = 0, UNIT, DISCONNECTION };
	typedef std::vector<std::vector<LinkType>> vvlt;
	vvlt lP;

private:
	double length;
	double scale = 4.0;

	void randomInit(TriMesh &mesh) override;

	void computeLinkProperty(vvlt &linkProperty);
	void computePolyFrame(PolyMesh &polymesh, std::vector<Vec3d> &pos, std::vector<VertexValue> &polyFrame);
};

#else
//class MeshField
//{
//public:
//	MeshField() {};
//	virtual ~MeshField() {};
//
//protected:
//
//	typedef unsigned int ui;
//	typedef Eigen::Matrix<double, 1, -1> vxu;
//	typedef Eigen::Vector3d Eigen::Vector3d;
//	typedef Eigen::VectorXd vxd;
//
//	typedef Eigen::Matrix<ui, 2, -1> m2xu;
//	//typedef Eigen::Matrix2Xd m2xd;
//	typedef Eigen::Matrix3Xd m3xd;
//
//	typedef Eigen::Vector3d Vec3d;
//
//	inline static Eigen::Vector3d OpenMesh2EigenVector(const Vec3d &pos)
//	{
//		return Eigen::Vector3d(pos[0], pos[1], pos[2]);
//	}
//	inline static Vec3d Eigen2OpenMeshVector(const Eigen::Vector3d &pos)
//	{
//		return Vec3d(pos(0), pos(1), pos(2));
//	}
//	//#define invalid (ui)-1
//		/*typedef Eigen::Vector3d Vec3d;
//		typedef std::vector<Vec3d> vv3d;*/
//	struct Link
//	{
//		ui id;
//		double weight;
//		Link() {}
//		Link(ui v_id, double w) :id(v_id), weight(w) {}
//		Link(ui v_id) :id(v_id), weight(1.0) {}
//		inline bool operator<(const Link &link)const { return id < link.id; }
//		inline bool operator==(const Link &link) const { return id == link.id; }
//	};
//	typedef std::vector<Link> vL;
//	typedef std::vector<vL> vvL;
//
//	struct edgeWeight
//	{
//		ui i; ui j;
//		double w;
//		edgeWeight(ui i_, ui j_, double w_) :i(i_), j(j_), w(w_) {}
//		inline bool operator>(const edgeWeight &right) const { return w > right.w; }
//		inline bool operator<(const edgeWeight &right) const { return w < right.w; }
//	};
//
//	struct hierarchy {
//		m3xd v;//顶点坐标
//		m3xd n;//法向
//		vxd a;//面积
//		std::vector<std::vector<Link>> adj;//邻接关系
//		vxu tocoarser;//寻找粗糙层指标
//		m2xu tofiner;//寻找精细层指标
//		struct hierarchy* coarser_hierarchy;//粗糙层指针
//		struct hierarchy* finer_hierarchy;//精细层指针
//		std::vector<std::vector<ui>> phase;//网格边染色
//		ui hierarchy_depth;//当前层数
//
//		m3xd o;//orientation
//		m3xd p;//position
//
//		//给定z轴，计算x轴
//		inline Eigen::Vector3d calcXAxis(const Eigen::Vector3d &z)
//		{
//			Eigen::Vector3d x;
//			using std::fabs;
//			if (fabs(z[0]) > fabs(z[1]))
//			{
//				x = Eigen::Vector3d(z[2], 0.0, -z[0]); 
//			}
//			else
//			{
//				x = Eigen::Vector3d(0.0, z[2], -z[1]);
//			}
//			x.normalize();
//			return x;
//		}
//
//		//由mesh直接创建depth层的粗糙网格
//		hierarchy(const TriMesh &mesh, ui depth);
//		~hierarchy();
//	private:
//		//计算Voronoi面积和邻接关系，一个较为简单的版本
//		double calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv);
//		//由mesh创建当前层网格
//		void init(const TriMesh &mesh);
//		//计算并行块
//		void graph_color();
//		//由当前精细层网格创建下一粗糙层网格
//		hierarchy(hierarchy* fh);
//	};
//	static struct hierarchy *h;
//
//	virtual void randomInit(TriMesh &mesh) = 0;
//	virtual void GSOptimize() = 0;
//public:
//	std::vector<std::vector<ui>> returnPhase() { return h->phase; }
//};
//
//#include <OpenMesh\Tools\Subdivider\Uniform\LoopT.hh>
//#include <list>
//class OrientationField :public MeshField
//{
//public:
//	OrientationField(TriMesh &mesh)
//	{
//		randomInit(mesh);
//	}
//	~OrientationField()
//	{
//		//if (h) { delete h; h = nullptr; }
//	};
//
//public:
//	void GSOptimize() override;
//	m3xd returnOrientation() { return h->o; }
//
//private:
//	void loopDivision(TriMesh &mesh, ui divideTimes = 0);
//	void randomInit(TriMesh &mesh) override;//初始化hierarchy h & orientation o
//};
//
//class PositionField :public MeshField
//{
//public:
//	PositionField(TriMesh &mesh)
//	{
//		randomInit(mesh);
//	}
//	~PositionField()
//	{
//		if (h) { delete h; h = nullptr; }
//	}
//
//public:
//	void GSOptimize() override;
//	m3xd returnPosition() { return h->p; }
//	struct VertexValue {
//		OV vert;
//		Eigen::Vector3d normal;
//		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//		std::set<ui> adj;
//		std::list<ui> orientedAdj;
//		VertexValue() {}
//		//VertexValue(OV &v, Eigen::Vector3d &n, std::set<ui> &a) :vert(v), normal(n), adj(a) {}
//		void setValue(OV &v, Eigen::Vector3d &n)
//		{
//			vert = v; normal = n; normal.normalize();
//		}
//	};
//	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//	//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
//	void extractMesh(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame);
//	enum LinkType { SAME = 0, UNIT, DISCONNECTION };
//	typedef std::vector<std::vector<LinkType>> vvlt;
//	vvlt lP;
//
//private:
//	double length;
//	double scale = 4.0;
//
//	void randomInit(TriMesh &mesh) override;
//
//	void computeLinkProperty(vvlt &linkProperty);
//	void computePolyFrame(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame);
//};
#endif

class MeshField
{
public:
	MeshField() {};
	virtual ~MeshField() {};

protected:
#ifdef USING_SMALLMESH
	typedef unsigned int ui;
#else 
	typedef unsigned long long ui;
#endif
#define invalid (ui)-1
	typedef Eigen::Matrix<double, 1, -1> vxu;
	typedef Eigen::Vector3d v3d;
	typedef Eigen::VectorXd vxd;

	typedef Eigen::Matrix<ui, 2, -1> m2xu;
	//typedef Eigen::Matrix2Xd m2xd;
	typedef Eigen::Matrix3Xd m3xd;

	typedef OpenMesh::Vec3d Vec3d;

	inline static Eigen::Vector3d OpenMesh2EigenVector(const Vec3d &pos)
	{
		return Eigen::Vector3d(pos[0], pos[1], pos[2]);
	}
	inline static Vec3d Eigen2OpenMeshVector(const Eigen::Vector3d &pos)
	{
		return Vec3d(pos(0), pos(1), pos(2));
	}
		/*typedef Eigen::Vector3d Vec3d;
		typedef std::vector<Vec3d> vv3d;*/
	template<typename T0, typename T1, typename T2>
	struct triple {
		T0 first; T1 second; T2 third;
		triple(){}
		triple(T0 f,T1 s,T2 t):first(f),second(s),third(t){}
	};
	struct Link
	{
		ui id;
		double weight;
		Link() {}
		Link(ui v_id, double w) :id(v_id), weight(w) {}
		Link(ui v_id) :id(v_id), weight(1.0) {}
		inline bool operator<(const Link &link)const { return id < link.id; }
		inline bool operator==(const Link &link) const { return id == link.id; }
	};
	typedef std::vector<Link> vL;
	typedef std::vector<vL> vvL;

	struct edgeWeight
	{
		ui i; ui j;
		double w;
		edgeWeight(ui i_, ui j_, double w_) :i(i_), j(j_), w(w_) {}
		inline bool operator>(const edgeWeight &right) const { return w > right.w; }
		inline bool operator<(const edgeWeight &right) const { return w < right.w; }
	};

	struct hierarchy {
		m3xd v;//顶点坐标
		m3xd n;//法向
		vxd a;//面积
		std::vector<std::vector<Link>> adj;//邻接关系
		vxu tocoarser;//寻找粗糙层指标
		m2xu tofiner;//寻找精细层指标
		struct hierarchy* coarser_hierarchy;//粗糙层指针
		struct hierarchy* finer_hierarchy;//精细层指针
		std::vector<std::vector<ui>> phase;//网格边染色
		ui hierarchy_depth;//当前层数

		m3xd o;//orientation
		m3xd p;//position

		//给定z轴，计算x轴
		inline Eigen::Vector3d calcXAxis(const Eigen::Vector3d &z)
		{
			Eigen::Vector3d x;
			using std::fabs;
			if (fabs(z[0]) > fabs(z[1]))
			{
				x = Eigen::Vector3d(z[2], 0.0, -z[0]);
			}
			else
			{
				x = Eigen::Vector3d(0.0, z[2], -z[1]);
			}
			x.normalize();
			return x;
		}

		//由mesh直接创建depth层的粗糙网格
		hierarchy(const TriMesh &mesh, ui depth);
		~hierarchy();
	private:
		//计算Voronoi面积和邻接关系，一个较为简单的版本
		double calcVoronoiAndAdjacency(const TriMesh &mesh, OV tv);
		//由mesh创建当前层网格
		void init(const TriMesh &mesh);
		//计算并行块
		void graph_color();
		//由当前精细层网格创建下一粗糙层网格
		hierarchy(hierarchy* fh);
	};
	static struct hierarchy *h;
	std::vector<ui> singularities;
	virtual void randomInit(TriMesh &mesh) = 0;
	virtual void GSOptimize() = 0;
public:
	std::vector<std::vector<ui>> returnPhase() { return h->phase; }
	std::vector<ui> returnSingularities() { return singularities; }
};

#include <OpenMesh\Tools\Subdivider\Uniform\LoopT.hh>
#include <list>
class OrientationField :public MeshField
{
public:
	OrientationField(TriMesh &mesh)
	{
		randomInit(mesh); mesh_ptr = &mesh;
	}
	~OrientationField()
	{
		//if (h) { delete h; h = nullptr; }
	};

public:
	void GSOptimize() override;
	m3xd returnOrientation() { return h->o; }

private:
	TriMesh* mesh_ptr;
	void loopDivision(TriMesh &mesh, ui divideTimes = 0);
	void randomInit(TriMesh &mesh) override;//初始化hierarchy h & orientation o
};

class PositionField :public MeshField
{
public:
	PositionField(TriMesh &mesh)
	{
		randomInit(mesh);
	}
	~PositionField()
	{
		if (h) { delete h; h = nullptr; }
	}

public:
	void GSOptimize() override;
	m3xd returnPosition() { return h->p; }
	struct VertexValue {
		OV vert;
		Eigen::Vector3d normal;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		std::set<ui> adj;
		std::list<ui> orientedAdj;
		VertexValue() {}
		//VertexValue(OV &v, Eigen::Vector3d &n, std::set<ui> &a) :vert(v), normal(n), adj(a) {}
		void setValue(OV &v, Eigen::Vector3d &n)
		{
			vert = v; normal = n; normal.normalize();
		}
	};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	//typedef std::pair<std::pair<Vec3d, Vec3d>, std::set<ui>> VertexValue;
	void extractMesh(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame);
	enum LinkType { SAME = 0, UNIT, DISCONNECTION };
	typedef std::vector<std::vector<LinkType>> vvlt;
	vvlt lP;

private:
	double length;
	double scale = 4.0;

	void randomInit(TriMesh &mesh) override;

	void computeLinkProperty(vvlt &linkProperty);
	void computePolyFrame(PolyMesh &polymesh, m3xd &pos, VertexValue* &polyFrame);

	//remove singularity
	void searchBFSTree(TriMesh &trimesh);
};


