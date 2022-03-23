#ifndef MESHDEFINITION_H
#define MESHDEFINITION_H

//#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

using namespace Eigen;

struct MeshTraits : public OpenMesh::DefaultTraits
{
	using Point = Vector3d;
	using Normal = Vector3d;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	FaceTraits
	{
		FaceT()
		{};
	};

	EdgeTraits
	{
		EdgeT()
		{};
	};

	HalfedgeTraits
	{
		HalfedgeT()
		{};
	};

	VertexTraits
	{
		VertexT() : new_pos_fixed(false)
		{
		};
	private:
		Vector3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
	public:
		void set_New_Pos(const Vector3d& n_p){new_pos = n_p;};
		Vector3d& get_New_Pos(){return new_pos;};

		void set_new_pos_fixed(bool f){new_pos_fixed = f;};
		bool get_new_pos_fixed(){return new_pos_fixed;};
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> PolyMesh;
typedef PolyMesh Mesh;


//double calc_mesh_ave_edge_length(Mesh* mesh_);

bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);//just copy the code from openmesh
bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);

bool check_in_triangle_face(const Matrix3Xd& tri, const Vector3d& p);
//bool baryCoord( const OpenMesh::Vec3d& _p, const OpenMesh::Vec3d& _u, const OpenMesh::Vec3d& _v, const OpenMesh::Vec3d& _w, OpenMesh::Vec3d&_result );
//
//void compute_point_area(Mesh* mesh_, std::vector<std::map<int,double>>& cornerArea, std::vector<double>& pointArea , bool use_np = false);
//
//void rot_coord_sys(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
//				   const OpenMesh::Vec3d &new_norm,
//				   OpenMesh::Vec3d &new_u, OpenMesh::Vec3d &new_v);
//
//void proj_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
//			   double old_ku, double old_kuv, double old_kv,
//			   const OpenMesh::Vec3d &new_u, const OpenMesh::Vec3d &new_v,
//			   double &new_ku, double &new_kuv, double &new_kv);
//
//// Given a curvature tensor, find principal directions and curvatures
//// Makes sure that pdir1 and pdir2 are perpendicular to normal
//void diagonalize_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v,
//					  double ku, double kuv, double kv,
//					  const OpenMesh::Vec3d &new_norm,
//					  OpenMesh::Vec3d &pdir1, OpenMesh::Vec3d &pdir2, double &vk1, double &vk2);
//
//void compute_principal_curvature(Mesh* mesh_, 
//								 std::vector<double>& K1, std::vector<double>& K2, 
//								 std::vector<OpenMesh::Vec3d>& dir1,std::vector<OpenMesh::Vec3d>& dir2);
//
//
//#pragma region functions by yanyisheshou at GCL
//double meshMinAngle(TriMesh &mesh);
//
//void printMeshQuality(TriMesh &mesh);
//
//template <typename T>
//void initMeshStatusAndNormal(T& m)
//{
//	m.request_vertex_status();
//	m.request_edge_status();
//	m.request_face_status();
//
//	m.request_face_normals();
//	m.request_vertex_normals();
//
//	m.update_face_normals();
//	m.update_vertex_normals();
//}
//
//template <typename T>
//bool isClosedMesh(T& mesh)
//{
//	for (auto tv : mesh.vertices())
//	{
//		if (mesh.is_boundary(tv))
//			return false;
//	}
//	return true;
//}
//
//template <typename T>
//double meshAverageLength(T &mesh)
//{
//	double le = 0;
//	for (auto te : mesh.edges())
//	{
//		le += mesh.calc_edge_length(te);
//	}
//	return le / mesh.n_edges();
//}
//#pragma endregion

#endif