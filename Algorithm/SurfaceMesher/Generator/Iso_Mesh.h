#pragma once
#include "OccReader.h"

namespace CADMesher
{
	class Iso_Mesh
	{
	public:
		explicit Iso_Mesh(QString &fileName);
		Iso_Mesh(const Iso_Mesh& im) = delete;
		~Iso_Mesh()
		{
			if (occ_reader) { delete occ_reader; occ_reader = nullptr; }
		};

	private:
		OccReader *occ_reader = nullptr;
		void InitTree();


	public:
		void MergeModel();
		int EndId(vector<ShapeEdge> &edgeshape, int edge_id);
		void ResetFeature();
		void ResetFeature1();

		template <typename T>
		std::string Write_Obj(T &aMesh)
		{
			std::ofstream file_writer;
		    Open_File(file_writer);
			
			for (auto tv = aMesh.vertices_begin(); tv != aMesh.vertices_end(); tv++)
			{
				auto pos = aMesh.point(*tv);
				file_writer << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			}
			for (auto tf = aMesh.faces_begin(); tf != aMesh.faces_end(); tf++)
			{
				file_writer << "f";
				for (auto tfv = aMesh.fv_begin(*tf); tfv != aMesh.fv_end(*tf); tfv++)
				{
					file_writer << " " << tfv->idx() + 1;
				}
				file_writer << "\n";
			}
			file_writer.close();
			return "step_to_obj.obj"; 
		}
		void Open_File(std::ofstream &file_writer);
	};
}
