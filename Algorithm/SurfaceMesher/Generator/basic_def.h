#pragma once

#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_OffsetSurface.hxx>

#include "..\src\MeshViewer\MeshDefinition.h"
#include "..\src\Dependency\CPS\CPS_AABBTree.h"
#include "..\src\Toolbox\dprinter\dprint.h"
#include"src\Dependency\BSpline\BSplineSurface.h"

namespace CADMesher 
{
#define PI 3.1415926535897932
#define epsilonerror 1.1e-15
	using namespace Eigen;
	using std::vector;
	
	struct ShapeFace
	{
		int id;
		TopoDS_Face face;
		vector<vector<int>> wires;
		BSplineSurface *Surface;
		ShapeFace(int id_, TopoDS_Face face_)
		{
			id = id_;
			face = face_;
		}
	};

	struct ShapeEdge
	{
		int id;
		bool if_merged;
		bool if_splitted;
		bool if_trimmed; 
		bool if_C0;
		bool if_curvature;
		TopoDS_Edge edge;
		int main_face;
		int secondary_face;
		int reversed_edge;
		int prev_edge;
		int next_reversed_edge;
		Matrix2Xd parameters;
		int begin_id;
		int end_id;

		int next_edge;
		bool if_visited;
		ShapeEdge(int id_, TopoDS_Edge edge_)
		{
			id = id_;
			edge = edge_;
			if_merged = false;
			if_splitted = false;
			if_trimmed = true;
			if_curvature = true;
			if_C0 = false;
			next_reversed_edge = -1;
			main_face = -1;
			secondary_face = -1;
			reversed_edge = -1;
			begin_id = -1;

			next_edge = -1;
			if_visited = false;
		}
	};

	struct GlobalGeometry {
		TopoDS_Shape aShape;
		vector<ShapeFace> faceshape;
		vector<ShapeEdge> edgeshape;
		vector<unsigned> triangle_surface_index;
		TriMesh initial_trimesh;
		TriMesh isotropic_trimesh;
		TriMesh anisotropic_trimesh;
		ClosestPointSearch::AABBTree* init_trimesh_tree = nullptr;
		ClosestPointSearch::AABBTree** init_surfacemesh_tree = nullptr;

		PolyMesh initial_polymesh;
		GlobalGeometry() {}
		void clear() {
			faceshape.clear();
			edgeshape.clear();
			triangle_surface_index.clear();
			initial_trimesh.clear();
			isotropic_trimesh.clear();
			anisotropic_trimesh.clear();

			if (init_trimesh_tree)
				init_trimesh_tree->clear();
			if (init_surfacemesh_tree)
			{
				delete[] init_surfacemesh_tree;
				init_surfacemesh_tree = nullptr;
			}
			initial_polymesh.clear();
		}
	};

	extern GlobalGeometry globalmodel;

	void MeshProjectToSurface(Mesh* mesh, vector<vector<unsigned>> &vertex_surface_index, GlobalGeometry* model);
}

