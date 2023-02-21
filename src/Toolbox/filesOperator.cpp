#include "filesOperator.h"
void getFiles(std::string path, std::vector<std::string>& files)
{
	//文件句柄
	intptr_t  hFile = 0;
	//文件信息
	struct _finddata_t fileinfo;
	std::string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之
			//如果不是,加入列表
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

#include <qstring.h>
void truncateFilePath(std::string& file)
{
	QString fileName = QString::fromStdString(file);
	int id = fileName.lastIndexOf("/");
	if (id == -1)
	{
		id = fileName.lastIndexOf("\\");
		if (id == -1)
			return;
	}
	file = fileName.right(fileName.length() - id - 1).toLatin1().data();
}

void truncateFileExtension(std::string& file)
{
	QString fileName = QString::fromStdString(file);
	int id = fileName.lastIndexOf(".");
	if (id == -1)
		return;
	fileName.truncate(id);
	file = fileName.toLatin1().data();
}

void truncateFileName(std::string& file)
{
	truncateFilePath(file);
	truncateFileExtension(file);
}


namespace LoopGen
{
	bool WriteEnergy(std::vector<double>& energy, std::string& model_name)
	{
		std::ofstream file_writer;
		file_writer.open("..//resource//energy//" + model_name + ".energy");
		if (file_writer.fail()) {
			return false;
		}
		for (auto e : energy)
			file_writer << e << "\n";
		file_writer.close();
		return true;
	}

	bool ReadEnergy(std::vector<double>& energy, std::string& model_name)
	{
		std::ifstream file_reader;
		file_reader.open("..//resource//energy//" + model_name + ".energy", std::ios::in);
		if (!file_reader.good())
			return false;
		char line[1024] = { 0 };
		int ii = 0;
		while (file_reader.getline(line, sizeof(line)))
		{
			std::stringstream num(line);
			num >> energy[ii];
			++ii;
		}
		file_reader.close();
		return true;
	}

	bool ReadRegion(std::vector<std::vector<int>> &vh_set, std::vector<OpenMesh::Vec3d> &dir, std::string& model_name)
	{
		std::ifstream file_reader;
		file_reader.open("..//resource//region//" + model_name + ".region", std::ios::in);
		if (!file_reader.good())
			return false;
		char line[1024] = { 0 };
		int n = 0;
		vh_set.clear();
		while (file_reader.getline(line, sizeof(line)))
		{
			std::stringstream num0(line);
			num0 >> n;
			std::vector<int> vhs(n);

			file_reader.getline(line, sizeof(line));
			std::stringstream num1(line);
			dir.push_back(OpenMesh::Vec3d());
			num1 >> dir.back()[0] >> dir.back()[1] >> dir.back()[2];

			for (int i = 0; i < n; ++i)
			{
				file_reader.getline(line, sizeof(line));
				std::stringstream num2(line);
				num2 >> vhs[i];
			}
			vh_set.push_back(std::move(vhs));
		}
		file_reader.close();
		return true;
	}
}