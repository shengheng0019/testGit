#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

vector<string> FindLinear(string filename)
{
	vector<string> namelist;
	std::ifstream file(filename);
	std::string line;
	std::getline(file, line);
	while (std::getline(file, line))
	{
		// 对每行字符串进行处理
		//std::cout << "Line: " << line << std::endl;

		// 获得第一个位置的元素
		std::string first_element;
		std::string second_element;
		std::istringstream ss(line);
		std::getline(ss, first_element, ',');
		std::getline(ss, second_element, ',');
		//std::cout << "First element: " << first_element<<" second element: "<<second_element << std::endl;
		//string totalname = Folder + "\\" + first_element + ".mps";
		//namelist.emplace_back(first_element, totalname);
		if (atoi(second_element.c_str())==0)
		{
			namelist.emplace_back(first_element);
			//std::cout << first_element <<endl;
		}
	}
	file.close();
	return namelist;
}