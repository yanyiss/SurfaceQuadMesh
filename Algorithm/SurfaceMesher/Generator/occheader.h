#pragma once

//OpenCascade Library
#pragma region OpenCascade Library
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_TangentialDeflection.hxx>

#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include <gp_Lin2d.hxx>
#include <Geom_Surface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>

#include <GCE2d_MakeSegment.hxx>

#include <TopoDS.hxx>
//#include <TopoDS_Edge.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>

#include <BRepLib.hxx>
#include <GeomConvert.hxx>
#include <Adaptor3d_HSurfaceTool.hxx>
#include <Geom2d_Geometry.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>

#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>


#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <BRepMesh_IncrementalMesh.hxx>


#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <GeomTools.hxx>
#include <GeomTools_CurveSet.hxx>
#include <GeomLProp_SLProps.hxx>
#pragma endregion

#include <QString.h>

//General Library
#include <stdlib.h>
#include <fstream>
#include <string>
#include <map>
#include <Eigen\Dense>

#include <random> 
#include <qstring.h>

#include "..\src\Toolbox\Math\TriangleInterface.h"

#include<limits.h>
using namespace Eigen;
