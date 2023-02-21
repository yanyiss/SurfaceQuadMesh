#pragma once
#include <string>
#include <vector>
#include <io.h>
void getFiles(std::string path, std::vector<std::string>& files);

#include <qstring.h>
void truncateFilePath(std::string& file);

void truncateFileExtension(std::string& file);

void truncateFileName(std::string& file);

#include <fstream>
#include "dprint.h"
#include "..//MeshViewer/MeshDefinition.h"
namespace LoopGen
{
	template <typename PlaneLoop>
	bool WritePlaneLoop(std::vector<PlaneLoop> &pls, std::string& model_name, Mesh *mesh)
	{
		std::ofstream file_writer;
		file_writer.open("..//resource//plane loop//" + model_name + ".pl");
		if (file_writer.fail()) {
			return false;
		}
		int lineNum = 0;
		for (int i = 0; i < pls.size(); ++i)
		{
			auto& pl = pls[i];
			file_writer << pl.size() << "\n";
			++lineNum;
			for (auto& plpl : pl)
			{
				file_writer << plpl.hl->id << " " << plpl.c << "\n";
			}
			lineNum += pl.size();
		}
		/*for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				auto& pl = InfoOnMesh[2 * i + j].pl;
				file_writer << pl.size() << "\n";
				++lineNum;
				for (auto& plpl : pl)
				{
					file_writer << plpl.hl->id << " " << plpl.c << "\n";
				}
				lineNum += pl.size();
			}
		}*/
		file_writer.close();

		file_writer.open("..//resource//plane loop//" + model_name + "_pl.txt");
		if (file_writer.fail()) {
			return false;
		}
		file_writer << lineNum;
		file_writer.close();
		return true;
	}

	template <typename M4, typename PlaneLoop>
	bool ReadPlaneLoop(M4 &m4, std::vector<PlaneLoop>& pls, std::string& model_name, Mesh *mesh)
	{
		std::ifstream file_reader;
		file_reader.open("..//resource//plane loop//" + model_name + "_pl.txt", std::ios::in);
		if (!file_reader.good())
			return false;
		int lineNum = 0;
		char line[1024] = { 0 };
		file_reader.getline(line, sizeof(line));
		std::stringstream num(line);
		num >> lineNum;
		file_reader.close();

		FILE* fp;
		int size = 0;
		char* ar;
		fp = fopen(("..//resource//plane loop//" + model_name + ".pl").c_str(), "r");
		if (NULL == fp)
		{
			return false;
		}
		fseek(fp, 0, SEEK_END);
		size = ftell(fp);
		rewind(fp);
		ar = (char*)malloc(sizeof(char) * size);
		fread(ar, 1, size, fp);//每次读一个，共读size次  
		std::stringstream data(ar);
		int ii = 0;
		int hid; double c;
		while (true)
		{
			int nn;
			data >> nn;
			--lineNum;
			pls[ii].clear();
			pls[ii].reserve(nn);
			for (int i = 0; i < nn; ++i)
			{
				data >> hid >> c;
				pls[ii].emplace_back(&m4.halfedgelayers[hid], c);
			}
			lineNum -= nn;
			++ii;
			if (lineNum <= 0)
				break;  
		}
		//while (true)
		//{
		//	int nn;
		//	data >> nn;
		//	--lineNum;
		//	InfoOnMesh[ii].pl.clear();
		//	InfoOnMesh[ii].pl.reserve(nn);
		//	for (int i = 0; i < nn; ++i)
		//	{
		//		data >> hid >> c;
		//		InfoOnMesh[ii].pl.emplace_back(&m4.halfedges[hid]/*mesh->halfedge_handle(hid)*/, c);
		//	}
		//	lineNum -= nn;
		//	++ii;
		//	if (lineNum <= 0)
		//		break;
		//}
		fclose(fp);
		return true;
	}

	bool WriteEnergy(std::vector<double>& energy, std::string& model_name);
	bool ReadEnergy(std::vector<double>& energy, std::string& model_name);

	template <typename cylinder, typename MatrixXd>
	bool WriteRegion(std::vector<cylinder> &cylinders, MatrixXd &crossfield, std::string& model_name)
	{
		std::ofstream file_writer;
		file_writer.open("..//resource//region//" + model_name + ".region");
		if (file_writer.fail()) {
			return false;
		}
		for (auto &cy : cylinders)
		{
			file_writer << cy.vertices.size() << "\n";
			int fid = cy.vertices.front()->hl->left;
			file_writer << crossfield(0, fid) << " " << crossfield(1, fid) << " " << crossfield(2, fid) << "\n";
			for (auto vl : cy.vertices)
			{
				file_writer << vl->v.idx() << "\n";
			}
		}
		file_writer.close();
		return true;
	}
	
	bool ReadRegion(std::vector<std::vector<int>> &vh_set, std::vector<OpenMesh::Vec3d> &dir, std::string& model_name);
}