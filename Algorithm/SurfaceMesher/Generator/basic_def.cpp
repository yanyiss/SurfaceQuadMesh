#pragma once
#include "basic_def.h"

#include "..\src\Dependency\BSpline\TestClosestPoint.h"

#include <ctime>

namespace CADMesher
{
	GlobalGeometry globalmodel = GlobalGeometry();

	void MeshProjectToSurface(Mesh* mesh, vector<vector<unsigned>> &vertex_surface_index, GlobalGeometry* model) {
		vector<ShapeFace> &faceshape = globalmodel.faceshape;
		vector<vector<unsigned>>::iterator index_itr = vertex_surface_index.begin();
		vector<int> surface_type(10, 0);

		int itertimes = 0;
		clock_t start = clock();
		//针对不同类型的曲面进行不同方法的投影
		for (auto surface : faceshape)
			//auto surface = faceshape[67];
		{
			TopLoc_Location loca;
			opencascade::handle<Geom_Surface> geom_surface = BRep_Tool::Surface(surface.face, loca);
			opencascade::handle<Standard_Type> type = geom_surface->DynamicType();

			if (type == STANDARD_TYPE(Geom_Plane)) {
				surface_type[0]++;
				dprint("Plane", itertimes++);
				opencascade::handle<Geom_Plane> geom_plane = Handle(Geom_Plane)::DownCast(geom_surface);
				Standard_Real a, b, c, d;
				geom_plane->Coefficients(a, b, c, d);

				OpenMesh::Vec3d p;
				if (fabs(a) > epsilonerror) p = { -d / a,0,0 };
				else if (fabs(b) > epsilonerror) p = { 0,-d / b,0 };
				else p = { 0,0,-d / c };
				OpenMesh::Vec3d n = OpenMesh::Vec3d(a, b, c).normalize();

				//print("pointnum:", index_itr->size());
				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					double d = (p - pos).dot(n);
					/*if (d > 0.072458)
						print("dist: ", id, d);
					else*/
					mesh->set_point(tv, d*n + pos);
				}
			}
			else if (type == STANDARD_TYPE(Geom_CylindricalSurface)) {
				surface_type[1]++;
				dprint("Cylinder", itertimes++);
				opencascade::handle<Geom_CylindricalSurface> geom_cylindricalsurface = opencascade::handle<Geom_CylindricalSurface>::DownCast(geom_surface);
				/*Standard_Real a1, a2, a3, b1, b2, b3, c1, c2, c3, d;
				geom_cylindricalsurface->Coefficients(a1, a2, a3, b1, b2, b3, c1, c2, c3, d);*/
				gp_Pnt loc = geom_cylindricalsurface->Location();
				OpenMesh::Vec3d p(loc.X(), loc.Y(), loc.Z());
				gp_Dir z_axis = geom_cylindricalsurface->Axis().Direction();
				OpenMesh::Vec3d n = OpenMesh::Vec3d(z_axis.X(), z_axis.Y(), z_axis.Z()).normalized();
				double radius = static_cast<double>(geom_cylindricalsurface->Radius());
				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					OpenMesh::Vec3d d = (pos - p).dot(n)*n + p;
					double lamda = radius / (d - pos).norm();
					mesh->set_point(tv, (1 - lamda)*d + lamda * pos);
				}
			}
			else if (type == STANDARD_TYPE(Geom_ConicalSurface)) {
				surface_type[2]++;
				dprint("Cone", itertimes++);
				opencascade::handle<Geom_ConicalSurface> geom_conicalsurface = opencascade::handle<Geom_ConicalSurface>::DownCast(geom_surface);

				////???????????????????????????????????????????????????????????????????????????????????????待定
				////可能需要用Apex()
				//gp_Pnt loc = geom_conicalsurface->Location();
				//OpenMesh::Vec3d p(loc.X(), loc.Y(), loc.Z());
				////??????????????????????????????????????????????????????????????????????????????????????????????????

				//gp_Dir z_axis = geom_conicalsurface->Axis().Direction();
				//OpenMesh::Vec3d n = OpenMesh::Vec3d(z_axis.X(), z_axis.Y(), z_axis.Z()).normalized();
				//double angle = static_cast<double>(geom_conicalsurface->SemiAngle());
				//for (auto id : *index_itr) {
				//	OV tv = mesh->vertex_handle(id);
				//	if (mesh->data(tv).get_vertflag()) continue;
				//	OpenMesh::Vec3d pos = mesh->point(tv);
				//}
			}
			else if (type == STANDARD_TYPE(Geom_SphericalSurface)) {
				surface_type[3]++;
				dprint("Sphere", itertimes++);
				opencascade::handle<Geom_SphericalSurface> geom_sphericalsurface = opencascade::handle<Geom_SphericalSurface>::DownCast(geom_surface);
				gp_Pnt loc = geom_sphericalsurface->Location();
				OpenMesh::Vec3d p(loc.X(), loc.Y(), loc.Z());
				double radius = static_cast<double>(geom_sphericalsurface->Radius());
				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					double lamda = radius / (p - pos).norm();
					mesh->set_point(tv, (1 - lamda)*p + lamda * pos);
				}
			}
			else if (type == STANDARD_TYPE(Geom_ToroidalSurface)) {
				surface_type[4]++;
				dprint("Torus", itertimes++);
				opencascade::handle<Geom_ToroidalSurface> geom_toroidalsurface = opencascade::handle<Geom_ToroidalSurface>::DownCast(geom_surface);
				double R = static_cast<double>(geom_toroidalsurface->MajorRadius());
				double r = static_cast<double>(geom_toroidalsurface->MinorRadius());
				gp_Pnt loc = geom_toroidalsurface->Location();
				OpenMesh::Vec3d p(loc.X(), loc.Y(), loc.Z());
				gp_Dir z_axis = geom_toroidalsurface->Position().Axis().Direction();
				OpenMesh::Vec3d n = OpenMesh::Vec3d(z_axis.X(), z_axis.Y(), z_axis.Z()).normalized();
				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					OpenMesh::Vec3d s = p + (pos - (pos - p).dot(n)*n - p).normalize()*R;
					double lamda = r / (s - pos).norm();
					mesh->set_point(tv, (1 - lamda)*s + lamda * pos);
				}
			}
			else if (type == STANDARD_TYPE(Geom_BezierSurface)) {
				surface_type[5]++;
				dprint("Bezier Surface", itertimes++);
				BezierSurface beziersurface;
				opencascade::handle<Geom_BezierSurface> geom_beziersurface = opencascade::handle<Geom_BezierSurface>::DownCast(geom_surface);

				Standard_Real ui, ua, vi, va;
				geom_beziersurface->Bounds(ui, ua, vi, va);

				TColgp_Array2OfPnt controlpoints = geom_beziersurface->Poles();
				vector<vector<Point>> cp(controlpoints.NbRows());
				for (int r = 1; r <= controlpoints.NbRows(); r++)
					for (int c = 1; c <= controlpoints.NbColumns(); c++) {
						gp_Pnt pos = controlpoints.Value(r, c);
						cp[r - 1].emplace_back(pos.X(), pos.Y(), pos.Z());
					}

				const TColStd_Array2OfReal* weights = geom_beziersurface->Weights();
				if (weights) {
					dprint(weights->NbRows(), weights->NbColumns());
					vector<vector<double>> w(weights->NbRows());
					for (int r = 1; r <= weights->NbRows(); r++)
					{
						w[r - 1].reserve(weights->NbColumns());
						for (int c = 1; c <= weights->NbColumns(); c++)
							w[r - 1].push_back(weights->Value(r, c));
					}
					beziersurface = { geom_beziersurface->UDegree(), geom_beziersurface->VDegree(), ui, ua, vi, va, w, cp };
				}
				else
					beziersurface = { geom_beziersurface->UDegree(), geom_beziersurface->VDegree(), ui, ua, vi, va, cp };

				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					vector<ProjectionPointToSurface> ppts;
					TestClosestPoint tcp(beziersurface, Point(pos[0], pos[1], pos[2]), ppts);
					/*if (ppts.empty()) continue;

					if (ppts.size() >= 2) continue;
	*/
					if (ppts.size() != 1) continue;
					ProjectionPointToSurface newpos = ppts.front();
					Point p = beziersurface(newpos.u, newpos.v);
					mesh->set_point(tv, OpenMesh::Vec3d(p(0), p(1), p(2)));
				}

			}
			else if (type == STANDARD_TYPE(Geom_BSplineSurface)) {
				surface_type[6]++;
				dprint("BSpline Surface", itertimes++);
				BSplineSurface bsplinesurface;

				opencascade::handle<Geom_BSplineSurface> geom_bsplinesurface = Handle(Geom_BSplineSurface)::DownCast(geom_surface);
				TColStd_Array1OfReal uknotsequence = geom_bsplinesurface->UKnotSequence();
				TColStd_Array1OfReal vknotsequence = geom_bsplinesurface->VKnotSequence();
				vector<double> u;
				vector<double> v;
				for (auto itr = uknotsequence.begin(); itr != uknotsequence.end(); itr++)
					u.push_back(*itr);
				for (auto itr = vknotsequence.begin(); itr != vknotsequence.end(); itr++)
					v.push_back(*itr);

				TColgp_Array2OfPnt controlpoints = geom_bsplinesurface->Poles();
				vector<vector<Point>> cp(controlpoints.NbRows());
				for (int r = 1; r <= controlpoints.NbRows(); r++)
				{
					cp[r - 1].reserve(controlpoints.NbColumns());
					for (int c = 1; c <= controlpoints.NbColumns(); c++) {
						gp_Pnt pos = controlpoints.Value(r, c);
						cp[r - 1].emplace_back(pos.X(), pos.Y(), pos.Z());
					}
				}


				const TColStd_Array2OfReal* weights = geom_bsplinesurface->Weights();
				if (weights) {
					dprint(weights->NbRows(), weights->NbColumns());
					vector<vector<double>> w(weights->NbRows());
					for (int r = 1; r <= weights->NbRows(); r++)
					{
						w[r - 1].reserve(weights->NbColumns());
						for (int c = 1; c <= weights->NbColumns(); c++)
							w[r - 1].push_back(weights->Value(r, c));
					}
					bsplinesurface = { geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, w, cp };
				}
				else
					bsplinesurface = { geom_bsplinesurface->UDegree(), geom_bsplinesurface->VDegree(), u, v, cp };

				for (auto id : *index_itr) {
					OV tv = mesh->vertex_handle(id);
					if (mesh->data(tv).get_vertflag()) continue;
					OpenMesh::Vec3d pos = mesh->point(tv);
					vector<ProjectionPointToSurface> ppts;
					TestClosestPoint tcp(bsplinesurface, Point(pos[0], pos[1], pos[2]), ppts);
					/*if (ppts.empty()) continue;

					if (ppts.size() >= 2) continue;
	*/
					if (ppts.size() != 1) continue;
					ProjectionPointToSurface newpos = ppts.front();
					Point p = bsplinesurface.DeBoor(newpos.u, newpos.v);
					mesh->set_point(tv, OpenMesh::Vec3d(p(0), p(1), p(2)));
				}
			}
			else if (type == STANDARD_TYPE(Geom_SurfaceOfRevolution)) {
				surface_type[7]++;
				dprint("Surface of Revolution", itertimes++);
			}
			else if (type == STANDARD_TYPE(Geom_SurfaceOfLinearExtrusion)) {
				surface_type[8]++;
				dprint("Surface of Extrusion", itertimes++);
			}
			else if (type == STANDARD_TYPE(Geom_OffsetSurface)) {
				surface_type[9]++;
				dprint("Offset Surface", itertimes++);
			}
			else {
				dprint("new face", itertimes++);
				assert(false);
			}
			index_itr++;
		}
		dprint("project time:", clock() - start, "ms");
		dprint("num of different surfaces in the model:\nPlane:", surface_type[0],
			"\nCylinder:", surface_type[1],
			"\nCone:", surface_type[2],
			"\nSphere:", surface_type[3],
			"\nTorus:", surface_type[4],
			"\nBezier:", surface_type[5],
			"\nBSpline:", surface_type[6],
			"\nSurface of Revolution:", surface_type[7],
			"\nSurface of Linear Extrusion:", surface_type[8],
			"\nOffset Surface:", surface_type[9],
			"\n");
	}
}


