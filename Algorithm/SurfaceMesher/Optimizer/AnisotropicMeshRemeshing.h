#ifndef ANISOTROPICMESHREMESHING_H
#define ANISOTROPICMESHREMESHING_H

//#include <QObject>
#include "..\src\MeshViewer\MeshDefinition.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "..\src\Dependency\Common\CommonDefinitions.h"
#include "CGAL_Definition.h"
#include "..\Generator\basic_def.h"
#if 0
namespace CADMesher
{
	class AnisotropicMeshRemeshing
	{
		//Q_OBJECT
	public:
		AnisotropicMeshRemeshing();
		~AnisotropicMeshRemeshing();

		void SetMesh(Mesh* mesh)
		{
			reset_all_State();
			mesh_ = mesh;
		}
		void reset_all_State();

		void load_ref_mesh(Mesh* aniso_ref_mesh);
		void sample_mesh_anisotropic_edge_length(double ref_edge_len = 1.0, double a = 1.5, bool add_flip = true);
		void do_remeshing(double ref_edge_len = 1.0, double a = 1.5);
		void build_AABB_tree_using_Ref();
		void build_AABB_tree_feature_edge_using_Ref();
		void calc_tri_quality();

		void LCOT_Optimize(int iter_num, double step_length);
		bool reposition_LCOT(double step_length);
		bool reposition_exp_LCOT(double step_length);
		bool flip_based_energy();

		void uniform_optimize(int iter_num, double step_length);

		void exp_mips_optimize(int iter_num, double area_angle_ratio, double energy_power);
		void reposition_exp_mips(double area_angle_ratio, double energy_power);

		void reposition_particle(double step_length);
		void flip_based_particle_energy();

		void delete_boundary_small_tri();
		double get_ref_mesh_ave_anisotropic_edge_length() { return ref_mesh_ave_anisotropic_edge_length; }

	private:
		bool draw_small_tri_ok;

		void project_on_reference_mesh_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void project_on_reference_edge_with_metric(Mesh::VertexHandle vh, OpenMesh::Vec3d& p);
		void find_nearst_point_on_reference_mesh(OpenMesh::Vec3d& p, bool is_boundary);
		void project_on_reference();
		void project_on_reference_new_p(std::vector<OpenMesh::Vec3d>& np);
		void project_on_reference(OpenMesh::Vec3d& p, OpenMesh::Vec3d& sp, OpenMesh::Vec3d& dir, double& dis);

		bool split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p);
		double calc_flip_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec3d& p3, const OpenMesh::Vec6d& M, bool use_area);
		double calc_flip_particle_energy(const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2, const OpenMesh::Vec6d& M);

		Mesh* mesh_; Mesh* ref_mesh_;
		std::vector<CGAL_3_Triangle> triangle_vectors;
		CGAL_AABB_Tree* AABB_tree;
		std::vector<CGAL_3_Segment> segment_vectors;
		std::vector<unsigned> segment_edge_id;
		CGAL_AABB_Segment_Tree* AABB_Segment_tree;

		std::vector<int> below_30_tri;
		std::vector<double> below_30_tri_angle;
		double smallest_angle_th;

		void compute_src_mesh_ave_anisotropic_edge_length();
		double ref_mesh_ave_anisotropic_edge_length;

		double compute_exp_mips_area_energy(Mesh::VertexHandle vh, OpenMesh::Vec3d& np,
			const std::vector<OpenMesh::Vec6d>& vH, double area_angle_ratio, double energy_power);

		double least_angle = PI / 90.0;
		double largest_angle = PI * 0.975;
	};
}
#endif
#endif