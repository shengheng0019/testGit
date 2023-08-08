#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>



std::vector<std::string> FindLinear(std::string filename)
{
	std::vector<std::string> namelist;
	std::ifstream file(filename);
	std::string line;
	std::getline(file, line);
    std::cout<<"open linear file success"<<std::endl;
	while (std::getline(file, line))
	{
		std::string first_element;
		std::string second_element;
		std::istringstream ss(line);
		std::getline(ss, first_element, ',');
		std::getline(ss, second_element, ',');
		//std::cout << "First element: " << first_element<<" second element: "<<second_element << std::endl;
		//string totalname = Folder + "\\" + first_element + ".mps";
		//namelist.emplace_back(first_element, totalname);
		if (atoi(second_element.c_str()) == 0)
		{
			namelist.emplace_back(first_element);
			//std::cout << first_element <<endl;
		}
	}
	file.close();

    // find the easy test example
    std::vector<std::string> integer_name_list;
    std::string integer_name;
    integer_name = "D:\\MIPsolver\\overall_algorithm_test_result\\test_example_0717.csv";

    std::ifstream file_inter(integer_name);
    while (std::getline(file_inter, line)) {
        std::string first_element;
        std::istringstream ss(line);
        std::getline(ss, first_element, ',');
        int n_count = std::count(namelist.begin(), namelist.end(), first_element);
        if (n_count == 1) {
            integer_name_list.push_back(first_element);
        }
    }
    file_inter.close();


    return integer_name_list;
}