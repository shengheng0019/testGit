//#pragma once
//#include "Partition.hpp"
//#include <map>
//#include <set>
//#include "runmain.hpp"
//
//std::map<int, block> DualBlock(std::vector<part::Partition> parts, std::vector<std::vector<int64_t>> matrix_dual)
//{
//    std::map<int, block> Blocks;
//    for (int i = 0; i < parts.size(); ++i)
//    {
//        auto temp_block = std::move(parts[i]);
//        std::set<int64_t> nodelist(temp_block.getNodes().begin(), temp_block.getNodes().end());
//        //cout << "第" << temp_block.getId() << "块:有" << nodelist.size() << "个变量" << endl;
//        block temp_b;
//        temp_b.Cons = nodelist;
//        Blocks.emplace(temp_block.getId(), temp_b);
//    }
//    block temp_b;
//    Blocks.emplace(parts.size(), temp_b);
//    //block检验
//    //cout << "block检验" << endl;
//    //for (int i = 0; i < Blocks.size() - 1; ++i)
//    //{
//    //    for (int j = i + 1; j < Blocks.size(); ++j)
//    //    {
//    //        set<int64_t> diff;
//    //        set_intersection(Blocks[i].Vars.begin(), Blocks[i].Vars.end(), Blocks[j].Vars.begin(), Blocks[j].Vars.end(), inserter(diff, diff.begin()));
//    //        if (!diff.empty())
//    //        {
//    //            cout << i << " " << j << "存在交集" << endl;
//    //        }
//    //    }
//    //}
//    //约束归类
//    //cout << "------------变量归类------------" << endl;
//    for (int i = 0; i < matrix_dual.size(); ++i)
//    {
//        std::set<int64_t> temp_cons(matrix_dual[i].begin(), matrix_dual[i].end());
//        //cout << temp_cons.size() << endl;
//        //cout << endl;
//        int mark = 0;
//        std::set<int64_t> diff_all;
//        for (int b = 0; b < parts.size(); ++b)
//        {
//            //bool issubset = includes(Blocks[b].Vars.begin(), Blocks[b].Vars.end(), temp_cons.begin(), temp_cons.end());
//            //cout << Blocks[b].Vars.size() << " ";
//            //if (issubset)
//            //{
//            //    Blocks[b].Cons.emplace(i);
//            //    mark = 1;
//            //    break;
//            //}
//            std::set<int64_t> diff;
//            set_difference(temp_cons.begin(), temp_cons.end(), Blocks[b].Cons.begin(), Blocks[b].Cons.end(), inserter(diff, diff.begin()));
//            if (diff.empty())
//            {
//                Blocks[b].Vars.emplace(i);
//                mark = 1;
//                break;
//            }
//        }
//        //cout << endl;
//        //for (auto all : diff_all)
//        //{
//        //    cout << all << " ";
//        //}
//        if (mark == 0)
//        {
//            Blocks[parts.size()].Vars.emplace(i);
//        }
//    }
//    //处理最后一个block的节点信息
//    for (auto it = Blocks[parts.size()].Vars.begin(); it != Blocks[parts.size()].Vars.end(); ++it)
//    {
//        auto nodelist = matrix_dual[*it];
//        for (auto value : nodelist)
//        {
//            Blocks[parts.size()].Cons.emplace(value);
//        }
//    }
//    //for (int i = 0; i < Blocks.size(); ++i)
//    //{
//    //    std::cout << endl << "block " << i << ":" << endl;
//    //    cout << "变量数:" << Blocks[i].Vars.size() << endl;
//    //    cout << "约束数:" << Blocks[i].Cons.size() << endl;
//    //    //for (auto var : Blocks[i].Vars)
//    //    //{
//    // //       cout << var << " ";
//    //    //}
//    // //   cout << endl;
//    //    //for (auto con : Blocks[i].Cons)
//    //    //{
//    // //       cout << con << " ";
//    //    //}
//    //    cout << endl;
//    //}
//    return Blocks;
//}


