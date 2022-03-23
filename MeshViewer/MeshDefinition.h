#ifndef MESHDEFINITION_H
#define MESHDEFINITION_H

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

typedef OpenMesh::PolyConnectivity::FaceIter faceiter;
typedef OpenMesh::PolyConnectivity::EdgeIter edgeiter;
typedef OpenMesh::PolyConnectivity::HalfedgeIter halfedgeiter;
typedef OpenMesh::PolyConnectivity::VertexIter vertexiter;
typedef OpenMesh::VertexHandle OV;
typedef OpenMesh::SmartVertexHandle SOV;
typedef OpenMesh::EdgeHandle OE;
typedef OpenMesh::SmartEdgeHandle SOE;
typedef OpenMesh::HalfedgeHandle OH;
typedef OpenMesh::SmartHalfedgeHandle SOH;
typedef OpenMesh::FaceHandle OF;
typedef OpenMesh::SmartFaceHandle SOF;
typedef OpenMesh::Vec3d O3d;

struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	FaceTraits
	{
		FaceT(): area( 0.0 ),new_area(0.0), w(0), faceflag(false)
		{};
	private:
		double area; double new_area;
		double w;
		bool faceflag ;
	public:
		void set_area(const double& a){area = a;};
		double& get_area(){return area;};

		void set_new_area(const double& a){new_area = a;};
		double& get_new_area(){return new_area;};

		void set_weight(const double& w_){w = w_;};
		double& get_weight(){return w;};

		void set_faceflag(bool b) { faceflag = b; };
		bool get_faceflag() { return faceflag; };
	};

	EdgeTraits
	{
		EdgeT() : is_feature(false), weight(1.0), edgeflag(false)
		{
		};
	private:
		bool is_feature; double weight; bool edgeflag = false;
	public:
		void set_edge_feature(bool b){is_feature = b;};
		bool get_edge_feature(){return is_feature;};
		bool flag1 = false;
		bool flag2 = false;//曲率特征
		
		void set_edge_weight(double w){ weight = w; };
		double get_edge_weight(){ return weight; };

		void set_edgeflag(bool b) { edgeflag = b; };
		bool get_edgeflag() { return edgeflag; };

		bool offsetflag = false;
		std::vector<O3d> offsetpnt;   //存放offset点
		std::vector<OV> offsetvh;     //存放offset点对应的handle

	};

	HalfedgeTraits
	{
		HalfedgeT() :opposite_angle(0.0), face_he_var(-1)
		{
		};
	private:
		double opposite_angle; int face_he_var;
	public:
		void set_opposite_angle(const double& a){opposite_angle = a;};
		double& get_opposite_angle(){return opposite_angle;};
		int get_face_he_var() { return face_he_var; };
		void set_face_he_var(int fhe) { face_he_var = fhe; };
	};

	VertexTraits
	{
		VertexT() : node_type(0), new_pos_fixed(false), color(1), density(1.0), h(OpenMesh::Vec6d(1, 0, 0, 1, 0, 1)), vertflag(false), targetlength(0)
		{
		};
	private:
		OpenMesh::Vec6d h;
		OpenMesh::Vec3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
		int node_type;
		int color;
		double density;
		bool vertflag;
		double targetlength;
	public:
		Eigen::Matrix2d M;   //存放第一基本量

		void set_Hessian(const OpenMesh::Vec6d& h_){h = h_;};
		OpenMesh::Vec6d& get_Hessian(){return h;};

		void set_New_Pos(const OpenMesh::Vec3d& n_p){new_pos = n_p;};
		OpenMesh::Vec3d& get_New_Pos(){return new_pos;};

		void set_node_type(int n){node_type = n;};
		int get_node_type(){return node_type;};

		void set_new_pos_fixed(bool f){new_pos_fixed = f;};
		bool get_new_pos_fixed(){return new_pos_fixed;};

		void set_color(int c){ color = c; };
		int get_color(){ return color; };

		void set_density(double d){ density = d; };
		double get_density(){ return density; };

		void set_vertflag(bool b) { vertflag = b; };
		bool get_vertflag() { return vertflag; };

		void set_targetlength(double d) { targetlength = d; };
		double get_targetlength() { return targetlength; };

	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> PolyMesh;
typedef PolyMesh Mesh;


double calc_mesh_ave_edge_length(Mesh* mesh_);

bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);//just copy the code from openmesh
bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);

bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p);
bool baryCoord( const OpenMesh::Vec3d& _p, const OpenMesh::Vec3d& _u, const OpenMesh::Vec3d& _v, const OpenMesh::Vec3d& _w, OpenMesh::Vec3d&_result );

void compute_point_area(Mesh* mesh_, std::vector<std::map<int,double>>& cornerArea, std::vector<double>& pointArea , bool use_np = false);

void rot_coord_sys(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
				   const OpenMesh::Vec3d &new_norm,
				   OpenMesh::Vec3d &new_u, OpenMesh::Vec3d &new_v);

void proj_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
			   double old_ku, double old_kuv, double old_kv,
			   const OpenMesh::Vec3d &new_u, const OpenMesh::Vec3d &new_v,
			   double &new_ku, double &new_kuv, double &new_kv);

// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
void diagonalize_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
					  double ku, double kuv, double kv,
					  const OpenMesh::Vec3d &new_norm,
					  OpenMesh::Vec3d &pdir1, OpenMesh::Vec3d &pdir2, double &vk1, double &vk2);

void compute_principal_curvature(Mesh* mesh_, 
								 std::vector<double>& K1, std::vector<double>& K2, 
								 std::vector<OpenMesh::Vec3d>& dir1,std::vector<OpenMesh::Vec3d>& dir2);


#pragma region functions by yanyisheshou at GCL
double meshMinAngle(TriMesh &mesh);

void printMeshQuality(TriMesh &mesh);

template <typename T>
void initMeshStatusAndNormal(T& m)
{
	m.request_vertex_status();
	m.request_edge_status();
	m.request_face_status();

	m.request_face_normals();
	m.request_vertex_normals();

	m.update_face_normals();
	m.update_vertex_normals();
}

template <typename T>
bool isClosedMesh(T& mesh)
{
	for (auto tv : mesh.vertices())
	{
		if (mesh.is_boundary(tv))
			return false;
	}
	return true;
}

template <typename T>
double meshAverageLength(T &mesh)
{
	double le = 0;
	for (auto te : mesh.edges())
	{
		le += mesh.calc_edge_length(te);
	}
	return le / mesh.n_edges();
}
#pragma endregion

#endif