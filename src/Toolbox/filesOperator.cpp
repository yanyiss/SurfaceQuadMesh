#include "filesOperator.h"
void getFiles(std::string path, std::vector<std::string>& files)
{
	//�ļ����
	intptr_t  hFile = 0;
	//�ļ���Ϣ
	struct _finddata_t fileinfo;
	std::string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//�����Ŀ¼,����֮
			//�������,�����б�
			if ((fileinfo.attrib &  _A_SUBDIR))
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

void WriteVector(std::string filename, std::vector<double> &data)
{
	std::fstream ofile(filename.c_str(), std::ios_base::out);
	int n = data.size();
	for (int i = 0; i < n; ++i)
	{
		ofile << data[i] << std::endl;
	}
	ofile.close();
}