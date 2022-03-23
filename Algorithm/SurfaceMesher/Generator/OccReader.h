#ifndef OCCREADER_H
#define OCCREADER_H
#include "occheader.h"
#include "basic_def.h"
#include"Riemann.h"

namespace CADMesher
{
	class OccReader {
	public:
		explicit OccReader(QString &fileName)
		{
			std::string filetype;
			if (fileName.endsWith(".stp") || fileName.endsWith(".STP") || fileName.endsWith(".STEP")) {
				reader = new STEPControl_Reader();
				dprint("CAD model from STEP file");
				filetype = "STEP";
			}
			else if (fileName.endsWith(".igs") || fileName.endsWith(".IGS") || fileName.endsWith(".IGES")) {
				reader = new IGESControl_Reader();
				dprint("CAD model from IGES file");
				filetype = "IGES";
			}
			else
			{
				dprinterror("Is anything wrong? It can't be other file types");
				exit(0);
			}

			dprint(filetype + "file read beginning!\n");
			reader->ReadFile(fileName.toLatin1().data());
			Standard_Integer NbTrans = reader->TransferRoots();
			globalmodel.aShape = reader->OneShape();
			dprint(filetype + "file read finished\n");

			ComputeFaceAndEdge();
			Discrete_Edge();
			Face_type();
			Trim_Edge();
			C0_Feature();
			curvature_feature();
		}
		OccReader(const OccReader& or) = delete;
		~OccReader() {
			if (reader) { delete reader; reader = nullptr; }
		}

	protected:
		XSControl_Reader *reader;


	public:
		double expected_edge_length;
		double mu = 1.2;//三角形面积允许扩张系数
		vector<TriMesh> Surface_TriMeshes;
		vector<PolyMesh> Surface_PolyMeshes;

		void ComputeFaceAndEdge();
		void Discrete_Edge();
		void Face_type();
		void Trim_Edge();
		void C0_Feature();
		void curvature_feature();
		void Surface_delete();
		void Set_TriMesh();
		void Set_PolyMesh();
		Matrix2Xd Subdomain(Matrix2Xd &all_pnts, vector<Matrix2Xi> &bnd, int &pointsnumber);
	private:
		double initialRate = 0.003;
		double degeneratedRate = 0.01;
	};
}
#endif // !OCCREADER_H
