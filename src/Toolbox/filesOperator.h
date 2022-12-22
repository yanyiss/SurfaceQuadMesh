#pragma once
#include <string>
#include <vector>
#include <io.h>
#include <fstream>
void getFiles(std::string path, std::vector<std::string>& files);

void WriteVector(std::string filename, std::vector<double> &data);