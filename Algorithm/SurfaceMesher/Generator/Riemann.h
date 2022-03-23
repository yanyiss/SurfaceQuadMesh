#include"src\Dependency\BSpline\BSplineSurface.h"
#include"src\MeshViewer\MeshDefinition.h"
class Riemannremesh
{
public:
	Riemannremesh(BSplineSurface *surface, TriMesh *paramesh)
	{
		B = surface;
		mesh = paramesh;
		uu = B->GetKnotsU();
		vv = B->GetKnotsV();
	}
	Eigen::Matrix2d Riemanndata(const TriMesh::Point &p);
	double Riemannlen(const TriMesh::VertexHandle &h1, const TriMesh::VertexHandle &h2);
	void split();
	void collapse();
	void flip();
	void updatepoint();
	void calulenth();
	void remesh();
	void curvature_feature(Eigen::Matrix2Xd &all_pnts, std::vector<bool> &curvature);   //ÇúÂÊÌØÕ÷¼ÆËã

	double highlenth;
	double lowlenth;
	std::vector<double> uu, vv;
	BSplineSurface *B;
	TriMesh *mesh;
};



