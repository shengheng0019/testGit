#include <chrono>
#include "Parallelfuncs.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include "linear.hpp"
#include "ADMM.h"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/basic_file_sink.h>
#include "generate_cut.hpp"
#include <mutex>

#ifdef _WIN32 // Windows

#include <Windows.h>

#elif __linux__ // Linux
#include <unistd.h>
#endif
//key: 算法+index, value: sol
std::map<std::string, writeJson::Soljson> all_solsmap;
std::string Folder;
std::string testExample;

int main(int argc, char** argv)
{

#ifdef _WIN32
    static DWORD processId = GetCurrentProcessId();
#elif __linux__
    pid_t processId = getpid();
#endif

    // log
    std::cout << "start   !\n" << std::endl;
    auto my_logger = spdlog::basic_logger_mt("basic_logger", "log/log_time.txt");
    my_logger->set_level(spdlog::level::info);
    my_logger->set_pattern("%v");
    my_logger->flush_on(spdlog::level::trace);


    std::thread([]
        {
            std::this_thread::sleep_for(7200s);

            // save test example name : more than 2 hours
            ofstream outF("log/moreTwoHours.txt", std::ios::app);
            outF << testExample << std::endl;
            outF.close();

            std::cerr << "TLE" << std::endl;

            std::vector<std::pair<std::string, writeJson::Soljson>> vec(all_solsmap.begin(), all_solsmap.end());

            std::sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {
                return a.second.obj < b.second.obj;
                });

            cout << vec[0].first << vec[0].second.obj;

#ifdef _WIN32 // Windows
            HANDLE hProcess = OpenProcess(PROCESS_TERMINATE, FALSE, processId);
            if (hProcess == NULL) {
                std::cout << "Failed to open process. Error code: " << GetLastError() << std::endl;
                return 1;
            }

            if (!TerminateProcess(hProcess, 0)) {
                std::cout << "Failed to terminate process. Error code: " << GetLastError() << std::endl;
                CloseHandle(hProcess);
                return 1;
            }

            std::cout << "Process terminated successfully." << std::endl;

            CloseHandle(hProcess);
#elif __linux__
            if (kill(processId, SIGTERM) == -1) {
                std::cout << "Failed to kill process." << std::endl;
                return 1;
            }

            std::cout << "Process killed successfully." << std::endl;
#endif
        }).detach();

        // for exe
        //std::cout << "argv[1] : " << argv[1] << std::endl;
        //std::cout << "argv[2] : " << argv[2] << std::endl;
        //Folder = argv[1];
        //testExample = argv[2];
        Folder = "D:\\MIPsolver\\collectionunzip"; 
        testExample = "neos-1367061"; 
        std::string tempname = (Folder + "\\" + testExample + ".mps");

        auto start = std::chrono::system_clock::now();

        int k = 6;

        auto [block_ori, scipmodel, scipvars, scipcons, dualblock, Var_Ni] = Dection(tempname, k);
        auto block_ = GenerateCut(block_ori);

        //If the size of var_ni is different from scipvars, complete var_ni
        for (int zz = 0; zz < scipvars.size(); zz++) {
            string name = SCIPvarGetName(scipvars[zz]);
            auto it = Var_Ni.find(name);
            if (it == Var_Ni.end()) {
                Var_Ni.insert(std::pair<std::string, int>(name, Var_Ni.size()));
            }
        }

        auto [vars, conss] = getAllVarsandCons(block_);
        SetVarsBounds(vars, conss);
        //update blocks
        for (int blkIdx = 0; blkIdx < block_.size(); blkIdx++) {
            for (int varIdx = 0; varIdx < block_[blkIdx].bVars.size(); varIdx++) {
                block_[blkIdx].bVars[varIdx].Lb = vars[block_[blkIdx].bVars[varIdx].Id].Lb;
                block_[blkIdx].bVars[varIdx].Ub = vars[block_[blkIdx].bVars[varIdx].Id].Ub;
            }
        }

        //    if (scipvars.size() == 0) {
        //        continue;
        //    }

            //@wqy add, diving using
        Objective master_obj;
        for (int j = 0; j < block_.size(); j++) {
            for (auto v : block_[j].bobj.coef) {
                master_obj.coef[v.first] = v.second;
            }
        }
        //the end of @wqy add

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = double(std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count()) / 1000;

        auto num_of_block = block_.size() - 1;
        auto coupling1 = 100 * double(block_[num_of_block].bCons.size()) / scipcons.size();
        //std::cout << Linear[i] << "," << duration << "," << scipcons.size() << "," << scipvars.size() << "," << block_[num_of_block].bCons.size() << "," << coupling1 << std::endl;
        //outfile << Linear[i] << "," << duration << "," << scipcons.size() << "," << scipvars.size() << "," << block_[num_of_block].bCons.size() << "," << coupling1<<"," << block_.size()<< std::endl;

        double IndependentNum = 0;
        for (auto tempBlock : dualblock) {
            IndependentNum += tempBlock.Independent_vars.size();
        }
        auto coupling2 = 100 * (scipvars.size() - IndependentNum) / scipvars.size();

        ofstream out_file("log/hypergraph.txt", std::ios::app);
        //算例名、超图划分算法时间、约束个数、变量个数、原始角块个数、对偶角块个数、耦合约束比例、耦合变量比例
        out_file << testExample << "\t" << duration << "\t" << scipcons.size() << "\t" << scipvars.size() << "\t"
            << block_.size() - 1 << "\t" << dualblock.size() << "\t" << coupling1 << "\t" << coupling2 << std::endl;
        out_file.close();

        ///
        ///=====================hypergraph partition end=====================
        ///
        //variables in ADMM will be changed after LG, must do it first

        /// ADMM start
        ADMMFunc(scipmodel, scipvars, dualblock, tempname, my_logger, start, Var_Ni, all_solsmap);
        /// ADMM end

        /// Lagrangean start
        LagrangeanFunc(block_, scipmodel, scipvars, scipcons, Var_Ni, tempname, testExample, start, my_logger);
        /// Lagrangean end

        ///
        ///=====================DW=====================
        ///
        DWFunc(block_, scipmodel, master_obj, all_solsmap, testExample, start);

        ///
        ///=====================Lagrangean end=====================
        ///
        std::cout << "algorithm all run time: " << double(std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count()) / 1000 << std::endl;

        my_logger->info("test_name:{}:run time:{}", testExample, double(std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start).count()) / 1000);

        for (auto iter = all_solsmap.begin(); iter != all_solsmap.end(); iter++) {
            std::cout << iter->first << " " << iter->second.obj << std::endl;
        }

}
