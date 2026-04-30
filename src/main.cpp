// Copyright 2020 The Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//
// The main file for running experiments
//

// compile this as C++14 (or later)

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <chrono>
#include <limits>

#include <algorithm>
#include <sys/_intsup.h>
#include "utilities.h"
#include "point.h"
#include "dataset.h"
#include "LDPMHQ.h"
#include "LDP_TLaplace.h"
#include "Skyline_fre.h"
#include "Skyline_dis.h"
#include "LDP_SPM.h"
#include "LDP_SW.h"
#include "Skyline_prior.h"
#include "sphere.h"

using std::ifstream;
using std::ostream;
using std::pair;
using std::unordered_set;
using std::vector;

template<class T, class R>
ostream &operator<<(ostream &s, const pair<T, R> &p) {
    return s << p.first << ' ' << p.second;
}


std::vector<double> Allocate_privacy_budgets(double epi, double Alpha, std::vector<double> weights) {
    double low = (1 - Alpha) * epi;
    double remain = weights.size() * Alpha * epi;

    std::vector<double> budgets;

    double min_weight = *std::min_element(weights.begin(), weights.end());
    // 调整权重
    std::vector<double> adjusted_weights;
    for (double w: weights) {
        adjusted_weights.push_back(w - min_weight);
    }
    double sum = std::accumulate(adjusted_weights.begin(), adjusted_weights.end(), 0.0);
    // 归一化
    std::vector<double> normalized_weights;
    for (double w: adjusted_weights) {
        if (w == 0) {
            normalized_weights.push_back(0.0);
        } else {
            normalized_weights.push_back(w / sum);
        }
    }
    double sum_weight = 0.0;
    for (double w: adjusted_weights) {
        sum_weight += w;
    }
    if (sum_weight == 0.0) {
        for (int i = 0; i < weights.size(); ++i) {
            budgets.push_back(epi);
        }
    } else {
        for (int i = 0; i < weights.size(); ++i) {
            budgets.push_back(low + normalized_weights[i] * remain);
        }
    }
    return budgets;
}

std::vector<double> getWeightsByDataset(const std::string &datasetname) {
    std::vector<double> weights;

    if (datasetname == "Compas") {
        // dimension = 9
        // weights = {
        //   0.05632,  // A1
        //   0.02492,  // A2
        //   0.27442,  // A3
        //   0.06112,  // A4
        //   0.13240,  // A5
        //   0.11205,  // A6
        //   0.05632,  // A7
        //   0.02492,  // A8
        //   0.25750   // A9
        // };
        weights = {
                0.2,
                0.2,
                0.2,
                0.2,
                0.2
        };
    } else if (datasetname == "Credit" || datasetname == "German" || datasetname == "GermanTest") {
        weights = {
                0.35734, 0.16446, 0.07403, 0.11401, 0.04864, 0.03329, 0.20819
        };
    } else if (datasetname == "Lawschs" || datasetname == "anti_2_10000") {
        weights = {
                0.42376, 0.57624
        };
    } else if (datasetname == "Anti-Cor_3_10000" || datasetname == "IND_3_10000" || datasetname == "Cor_3_10000") {
        weights = {
                0.63334572, 0.26049796, 0.10615632
        };
    } else if (datasetname == "Anti-Cor_4_10000" || datasetname == "IND_4_10000" || datasetname == "Cor_4_10000" ||
               datasetname == "Cor_4_5000") {
        weights = {
                0.12345, 0.28761, 0.19876, 0.39018
        };
    } else if (datasetname == "Anti-Cor_5_10000" || datasetname == "Cor_5_10000") {
        weights = {
                0.08765, 0.21345, 0.17654, 0.29876, 0.2236
        };
    } else if (datasetname == "Anti-Cor_6_10000" || datasetname == "House") {
        weights = {
                0.05432, 0.18765, 0.12345, 0.23456, 0.19876, 0.20126
        };
    } else if (datasetname == "Anti-Cor_7_10000") {
        weights = {
                0.04321, 0.12345, 0.09876, 0.18765, 0.15432, 0.21345, 0.17916
        };
    } else if (datasetname == "anti_4_10000") {
        weights = {
                0.27315, 0.18642, 0.34108, 0.19935
        };
    } else if (datasetname == "Adult" || datasetname == "Adult_age") {
        weights = {
                0.28376, 0.12654, 0.22018, 0.28376, 0.08576
        };
    } else {
        // 默认返回空向量或抛出异常，根据需求选择
        // throw std::invalid_argument("Unknown dataset name: " + datasetname);
    }

    return weights;
}

void applyLDPProtection(TLaplace &tlaplace,
                        SPM &spm,
                        SquareWaveMechanism sw,
                        const std::string &LDP,
                        int dim,
                        LDPHMQ &func_ldpmhq) {

    if (LDP == "TLaplace") {
        func_ldpmhq.dataset_.LDP_points_.clear();
        for (Point p: func_ldpmhq.dataset_.points_) {
            vector<double> values = p.get_coordinates();
            // for (int i = 0; i < dim; i++) {
            //   values.push_back(p.coordinates[i]);
            // }
            vector<double> ldp_values = tlaplace.randomize_batch(values);
            Point ldp_point = p;
            ldp_point.coordinates = ldp_values;
            func_ldpmhq.dataset_.LDP_points_.push_back(ldp_point);
        }
        func_ldpmhq.data_to_LDP();
    } else if (LDP == "SPM") {
        func_ldpmhq.dataset_.LDP_points_.clear();
        for (Point p: func_ldpmhq.dataset_.points_) {
            vector<double> values;
            for (int i = 0; i < dim; i++) {
                values.push_back(p.coordinates[i]);
            }
            vector<double> ldp_values = spm.randomize_batch(values);
            Point ldp_point = p;
            ldp_point.coordinates = ldp_values;
            func_ldpmhq.dataset_.LDP_points_.push_back(ldp_point);
        }
        func_ldpmhq.data_to_LDP();
    } else if (LDP == "SW") {
        func_ldpmhq.dataset_.LDP_points_.clear();
        for (Point p: func_ldpmhq.dataset_.points_) {
            vector<double> values;
            for (int i = 0; i < dim; i++) {
                values.push_back(p.coordinates[i]);
            }
            vector<double> ldp_values = sw.randomize_batch(values);
            Point ldp_point = p;
            ldp_point.coordinates = ldp_values;
            func_ldpmhq.dataset_.LDP_points_.push_back(ldp_point);
        }
        func_ldpmhq.data_to_LDP();
    } else {
        std::cout << "Unknown LDP type: " << LDP << std::endl;
    }
}

vector<double> Find_Max_attribute(const vector<Point> &Points) {
    if (Points.empty()) {
        return vector<double>();
    }

    // 获取第一个点的维度数
    size_t dimension = Points[0].get_dimension();

    // 初始化最大值向量，初始值为最小可能的double
    vector<double> max_values(dimension, 1);

    // 遍历所有点，更新每个维度的最大值
    for (const auto &point: Points) {
        // 确保所有点的维度一致
        if (point.get_dimension() != dimension) {
            // 如果维度不一致，返回空向量
            return vector<double>();
        }
        vector<double> coords = point.get_coordinates();
        for (size_t i = 0; i < dimension; ++i) {
            if (coords[i] > max_values[i]) {
                max_values[i] = coords[i];
            }
        }
    }
    return max_values;
}


double euclideanDistance(const Point& a, const Point& b) {
    const vector<double>& ca = a.coordinates;
    const vector<double>& cb = b.coordinates;

    if (ca.size() != cb.size()) {
        return std::numeric_limits<double>::infinity();
    }

    double sum = 0.0;
    for (size_t i = 0; i < ca.size(); ++i) {
        double diff = ca[i] - cb[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

vector<Point> getSphereSolutionPoints(const vector<Point>& points, int k) {
    vector<Point> solution_points;

    if (points.empty()) return solution_points;

    point_set_t* points_convert = PointVectorToPointSet(points);
    point_set_t* points_sphere = sphereWSImpLP(points_convert, k);

    // 先把原 points 建一个 id -> Point 的映射
    std::unordered_map<int, Point> id_to_point;
    for (const auto& p : points) {
        id_to_point[p.id] = p;
    }

    for (int i = 0; i < points_sphere->numberOfPoints; i++) {
        if (points_sphere->points[i]) {
            int id = points_sphere->points[i]->id;
            if (id_to_point.find(id) != id_to_point.end()) {
                solution_points.push_back(id_to_point[id]);
            }
        }
    }

    // 如果你的工程里 release_point_set 可用，建议打开
    // release_point_set(points_sphere, false);
    // release_point_set(points_convert, true);

    return solution_points;
}

vector<Point> mapToOriginalPoints(const vector<Point>& selected_points,
                                  const vector<Point>& original_points) {
    vector<Point> mapped_points;

    std::unordered_map<int, Point> id_to_original;
    for (const auto& p : original_points) {
        id_to_original[p.id] = p;
    }

    for (const auto& p : selected_points) {
        auto it = id_to_original.find(p.id);
        if (it != id_to_original.end()) {
            mapped_points.push_back(it->second);
        }
    }

    return mapped_points;
}

//double calculateAMD(const vector<Point>& S, const vector<Point>& S_mapped) {
//    if (S.empty() || S_mapped.empty()) {
//        return std::numeric_limits<double>::infinity();
//    }
//
//    double total = 0.0;
//
//    for (const auto& q : S_mapped) {
//        double min_dist = std::numeric_limits<double>::infinity();
//        for (const auto& p : S) {
//            double dist = euclideanDistance(p, q);
//            if (dist < min_dist) {
//                min_dist = dist;
//            }
//        }
//        total += min_dist;
//    }
//
//    // 按你的定义分母是 k
//    return total / static_cast<double>(S_mapped.size());
//}
double calculateAMD(const vector<Point>& S, const vector<Point>& S_mapped) {
    if (S.empty() || S_mapped.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    double total_1 = 0.0; // S'' -> S
    for (const auto& q : S_mapped) {
        double min_dist = std::numeric_limits<double>::infinity();
        for (const auto& p : S) {
            double dist = euclideanDistance(p, q);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
        total_1 += min_dist;
    }

    double total_2 = 0.0; // S -> S''
    for (const auto& p : S) {
        double min_dist = std::numeric_limits<double>::infinity();
        for (const auto& q : S_mapped) {
            double dist = euclideanDistance(p, q);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
        total_2 += min_dist;
    }

    double avg_1 = total_1 / static_cast<double>(S_mapped.size());
    double avg_2 = total_2 / static_cast<double>(S.size());

    return (avg_1 + avg_2) / 2.0;
}

// 计算points的mhr
double calculateMHR(LDPHMQ & ldpmhq, const vector<Point>& points, int k) {
    // cout << "the id of points:" << endl;
    // for (auto point : points) {
    //   cout << point.id << " ";
    // }
    // cout << endl;
    point_set_t* points_convert = PointVectorToPointSet(points);
    // cout << "the id of points_convert：" << endl;
    // for (int i = 0; i < points_convert->numberOfPoints; i++) cout << points_convert->points[i]->id << " ";
    // cout << endl;

    point_set_t* points_sphere = sphereWSImpLP(points_convert, k);
    vector<pair<int, int>> solutionVector;
    solutionVector.clear();
    cout << "the id of sphere_points :" << endl;
    for (int i = 0; i < points_sphere->numberOfPoints; i++) {
        if (points_sphere->points[i]) {
            // 第一个值为索引(从1开始)，第二个值为点的id
            cout << points_sphere->points[i]->id << " ";
            solutionVector.push_back(make_pair(i + 1, points_sphere->points[i]->id));
        }
    }
    // cout << endl;
    double mhr = ldpmhq.get_realmhr(solutionVector);
    // release_point_set(points_sphere, false);
    // release_point_set(points_convert, true);
    return mhr;
}


// 先Loop再Skyline三种方法
int main() {
    RandomHandler::CheckRandomNumberGenerator();

//      LDPHMQ func_ldpmhq("Anti-Cor_4_10000", 100);
//       LDPHMQ func_ldpmhq("smart-home",100);
//       LDPHMQ func_ldpmhq("smart-home-skyline",100);
       LDPHMQ func_ldpmhq("wearable",100);
//    LDPHMQ func_ldpmhq("healthcare", 100);


    //================================================================================ ==========================================================================================
    std::string algorithm = "sphere";

//    double Alpha = 0.0; //个性化参数设置
    int loops = 50; //循环次数设置
//    string choose = "percent";
//    string choose = "epsilon";
    string choose = "k";


    //找出每个属性的最大值
    vector<double> max_atr;
    max_atr = Find_Max_attribute(func_ldpmhq.dataset_.points_);
    func_ldpmhq.Max_attr = max_atr;


    vector<int> percent_vector; vector<int> epsilon_vector; vector<int> k_vector;
    if (choose == "percent"){
        percent_vector = {10, 20, 40, 80, 100, 111};
        epsilon_vector = {7};  // 3 5 7
        k_vector = {20};
    }else if (choose == "epsilon"){
        epsilon_vector = {1, 3, 5, 7, 9};
        percent_vector = {20};
        k_vector = {25};     // 15 20 25
    }else {
        k_vector = {10, 15, 20, 25, 30};
        percent_vector = {80};   // 20 40 80
        epsilon_vector = {5};
    }


//    const vector<std::string> skyline_name = {"Pri","Fre", "Dis"};
    const vector<std::string> skyline_name = {"Dis","Pri", "Fre"};
//    const vector<std::string> LDP_name = {"SPM", "SW", "TLaplace"};
    const vector<std::string> LDP_name = {"TLaplace", "SPM", "SW"};


    //各个数据集的权重
    std::string datasetname = func_ldpmhq.get_name(); // 记录数据集名称
    int dim = func_ldpmhq.dataset_.dim_;

    std::string outfilePath;
    if (datasetname == "wearable") {
        if (choose == "percent")
//            outfilePath = "../results/CD/wearable/percent/from[5,10,20,50,100,total]_epsilon3_k20.txt";
//            outfilePath = "../results/CD/wearable/percent/from[5,10,20,50,100,total]_epsilon5_k20.txt";
            outfilePath = "../results/CD/wearable/percent/from[5,10,20,50,100,total]_epsilon7_k20.txt";
        else if (choose == "epsilon")
//            outfilePath = "../results/CD/wearable/epsilon/from[1,3,5,7,9]_percent20_k15.txt";
//            outfilePath = "../results/CD/wearable/epsilon/from[1,3,5,7,9]_percent20_k20.txt";
            outfilePath = "../results/CD/wearable/epsilon/from[1,3,5,7,9]_percent20_k25.txt";
//        else outfilePath = "../results/CD/wearable/k/from[10,15,20,25,30]_percent20_epi5.txt";
//        else outfilePath = "../results/CD/wearable/k/from[10,15,20,25,30]_percent40_epi5.txt";
        else outfilePath = "../results/CD/wearable/k/from[10,15,20,25,30]_percent80_epi5.txt";
    } else if (datasetname == "smart-home") {
        if (choose == "percent")
            outfilePath = "../results/CD/smart-home/percent/from[5,10,20,50,100,total]_epsilon3_k20.txt";
//            outfilePath = "../results/CD/smart-home/percent/from[5,10,20,50,100,total]_epsilon5_k20.txt";
//            outfilePath = "../results/CD/smart-home/percent/from[5,10,20,50,100,total]_epsilon7_k20.txt";
        else if (choose == "epsilon")
//            outfilePath = "../results/CD/smart-home/epsilon/from[1,3,5,7,9]_percent20_k15.txt";
//            outfilePath = "../results/CD/smart-home/epsilon/from[1,3,5,7,9]_percent20_k20.txt";
            outfilePath = "../results/CD/smart-home/epsilon/from[1,3,5,7,9]_percent20_k25.txt";
//        else outfilePath = "../results/CD/smart-home/k/from[10,15,20,25,30]_percent20_epi5.txt";
//        else outfilePath = "../results/CD/smart-home/k/from[10,15,20,25,30]_percent40_epi5.txt";
        else outfilePath = "../results/CD/smart-home/k/from[10,15,20,25,30]_percent80_epi5.txt";
    } else if (datasetname == "Anti-Cor_4_10000") {
        if (choose == "percent")
//            outfilePath = "../results/CD/ANTI/percent/from[5,10,20,50,100,total]_epsilon3_k20.txt";
//            outfilePath = "../results/CD/ANTI/percent/from[5,10,20,50,100,total]_epsilon5_k20.txt";
            outfilePath = "../results/CD/ANTI/percent/from[5,10,20,50,100,total]_epsilon7_k20.txt";
        else if (choose == "epsilon")
//            outfilePath = "../results/CD/ANTI/epsilon/from[1,3,5,7,9]_percent20_k15.txt";
//            outfilePath = "../results/CD/ANTI/epsilon/from[1,3,5,7,9]_percent20_k20.txt";
            outfilePath = "../results/CD/ANTI/epsilon/from[1,3,5,7,9]_percent20_k25.txt";
        else outfilePath = "../results/CD/ANTI/k/from[10,15,20,25,30]_percent20_epi5.txt";
//        else outfilePath = "../results/CD/ANTI/k/from[10,15,20,25,30]_percent40_epi5.txt";
//        else outfilePath = "../results/CD/ANTI/k/from[10,15,20,25,30]_percent80_epi5.txt";
    }  else if (datasetname == "healthcare") {
        if (choose == "percent")
//            outfilePath = "../results/CD/healthcare/percent/from[5,10,20,50,100,total]_epsilon3_k20.txt";
//            outfilePath = "../results/CD/healthcare/percent/from[5,10,20,50,100,total]_epsilon5_k20.txt";
            outfilePath = "../results/CD/healthcare/percent/from[5,10,20,50,100,total]_epsilon7_k20.txt";
        else if (choose == "epsilon")
//            outfilePath = "../results/CD/healthcare/epsilon/from[1,3,5,7,9]_percent20_k15.txt";
//            outfilePath = "../results/CD/healthcare/epsilon/from[1,3,5,7,9]_percent20_k20.txt";
            outfilePath = "../results/CD/healthcare/epsilon/from[1,3,5,7,9]_percent20_k25.txt";
//        else outfilePath = "../results/CD/healthcare/k/from[10,15,20,25,30]_percent20_epi5.txt";
//        else outfilePath = "../results/CD/healthcare/k/from[10,15,20,25,30]_percent40_epi5.txt";
        else outfilePath = "../results/CD/healthcare/k/from[10,15,20,25,30]_percent80_epi5.txt";
    }

//    outfilePath = "../results/CD/test.txt";
    std::ofstream fCD(outfilePath, std::ios::out);
    if (!fCD.is_open()) {
        cout << outfilePath + "文件无法打开！" << std::endl;
        return 0;
    }


    func_ldpmhq.LDP_to_data();

    vector<vector<vector<double>>> CD_total;     // Loop * percent_from1 * skyline_name 的vector


    size_t maxSkylineNum = 0;
    // 每一次LDP都需要循环Loop次，CD_total保存Loop次的结果
    for (auto ldp: LDP_name) {
        cout << "ldp: " << ldp << endl;
        CD_total.clear();
        CD_total.resize(loops);

        if (choose == "epsilon") {
            for (size_t i = 0; i < loops; i++) {
                CD_total[i].resize(epsilon_vector.size());
                for (size_t j = 0; j < epsilon_vector.size(); j++) {
                    CD_total[i][j].resize(skyline_name.size());
                }
            }
        }else if (choose == "percent"){
            for (size_t i = 0; i < loops; ++i) {
                CD_total[i].resize(percent_vector.size());
                for (size_t j = 0; j < percent_vector.size(); j++) {
                    CD_total[i][j].resize(skyline_name.size());
                }
            }
        }else{
            for (size_t i = 0; i < loops; i++) {
                CD_total[i].resize(k_vector.size());
                for (size_t j = 0; j < k_vector.size(); j++) {
                    CD_total[i][j].resize(skyline_name.size());
                }
            }
        }


        for (int loop = 0; loop < loops; loop++) {
            //清空LDP点数据集
            func_ldpmhq.dataset_.LDP_points_.clear();
            //情况skyline数据集
            func_ldpmhq.dataset_.Sky_points_.clear();
            vector<double> eps;  //每个属性的隐私预算
            for (int e = 0; e < epsilon_vector.size(); e++) {
                eps.clear(); //清空上一次循环的数据
                eps.assign(dim, epsilon_vector[e]);
//            eps = Allocate_privacy_budgets(epi_fix[e], Alpha,weights);
                // std::cout<< "===============================================================================" << std::endl;
                TLaplace tlaplace(eps, max_atr);
                SPM spm(eps, max_atr);
                SquareWaveMechanism sw(eps, max_atr);
                applyLDPProtection(tlaplace, spm, sw, ldp, dim, func_ldpmhq);

                for (int k = 0; k < k_vector.size(); k++) {
//                    double real_mhr = calculateMHR(func_ldpmhq, func_ldpmhq.dataset_.points_, k_vector[k]);
//                    cout << "real_mhr: " << real_mhr << endl;
                    vector<Point> S = getSphereSolutionPoints(func_ldpmhq.dataset_.points_, k_vector[k]);
                    cout << "real solution size: " << S.size() << endl;


                    // 如果choose == percent 则计算percent=100 的skyline解集数组，然后再依次取前百分之几
                    vector<int> Fre, Pri, Dis;
                    if (choose != "percent")
                        computerSkylineBNL(func_ldpmhq.dataset_.LDP_points_, dim, Fre, Pri, Dis, percent_vector[0]);
                    else computerSkylineBNL(func_ldpmhq.dataset_.LDP_points_, dim, Fre, Pri, Dis, 100);
                    maxSkylineNum = max(Fre.size(), maxSkylineNum);
                    for (int i = 0; i < percent_vector.size(); i++) {
                        for (int j = 0; j < skyline_name.size(); j++) {
//                            double LDP_mhr = 0.0;
                            vector<Point> S_prime;
                            vector<Point> S_double_prime;   // mapped back to original dataset
                            double amd = 0.0;
                            int fre_candidates = 0, pri_candidates = 0, dis_candidates = 0;
                            std::cout
                                    << "==============================================================================================="
                                    << endl;
                            std::cout << "loop: " << loop << " LDP: " << ldp << " epi: " << epsilon_vector[e] << " k: "
                                      << k_vector[k] << " percent: " << percent_vector[i] << " skyline: " << skyline_name[j]
                                      << std::endl;

                            if (percent_vector[i] != 111) {
                                func_ldpmhq.dataset_.Sky_points_.clear();
                                vector<int> indices;
                                indices.clear();

                                if (skyline_name[j] == "Fre") {
                                    fre_candidates = static_cast<int>(Fre.size() * percent_vector[i] / 100.0);
                                    indices.assign(Fre.begin(), Fre.begin() + fre_candidates);
                                } else if (skyline_name[j] == "Pri") {
                                    pri_candidates = static_cast<int>(Pri.size() * percent_vector[i] / 100.0);
                                    indices.assign(Pri.begin(), Pri.begin() + pri_candidates);
                                } else {
                                    dis_candidates = static_cast<int>(Dis.size() * percent_vector[i] / 100.0);
                                    indices.assign(Dis.begin(), Dis.begin() + dis_candidates);
                                }

                                cout << "the num of skyline_points: " << indices.size() << endl;
                                for (auto x: indices) {
                                    cout << x << " ";
                                    Point sky_point = func_ldpmhq.dataset_.LDP_points_.at(x);
                                    func_ldpmhq.dataset_.Sky_points_.push_back(sky_point);
                                }
                                cout << endl;
                            }

                            // total: 在整个加噪数据集上执行 HMQ
                            // 其他 percent: 在 skyline / representative skyline 上执行 HMQ
                            if (percent_vector[i] == 111) {
//                                LDP_mhr = calculateMHR(func_ldpmhq, func_ldpmhq.dataset_.LDP_points_, k_vector[k]);
                                S_prime = getSphereSolutionPoints(func_ldpmhq.dataset_.LDP_points_, k_vector[k]);
                            } else {
//                                LDP_mhr = calculateMHR(func_ldpmhq, func_ldpmhq.dataset_.Sky_points_, k_vector[k]);
                                S_prime = getSphereSolutionPoints(func_ldpmhq.dataset_.Sky_points_, k_vector[k]);
                            }

                            // 按 id 映射回原始数据中的 S''
                            S_double_prime = mapToOriginalPoints(S_prime, func_ldpmhq.dataset_.points_);
                            amd = calculateAMD(S, S_double_prime);

                            cout << "amd: " << amd << endl;
                            if (choose == "percent") CD_total[loop][i][j] = amd;            // percent: i    skyline: j
                            else if (choose == "k") CD_total[loop][k][j] = amd;
                            else CD_total[loop][e][j] = amd;
                        }
                    }
                }
            }
        } // Loop结束，输出

        // --------------输出percent数组--------
        if (choose == "percent") {
            vector<double> Dis(percent_vector.size(), 0);  // 大小5，初始值0
            vector<double> Pri(percent_vector.size(), 0);
            vector<double> Fre(percent_vector.size(), 0);
            for (auto item: CD_total) {
                for (int i = 0; i < percent_vector.size(); i++) {
                    Dis[i] += item[i][0];
                    Pri[i] += item[i][1];
                    Fre[i] += item[i][2];
                }
            }
            cout << "Distance (CD)-------------" << endl;
            cout << "percent 10:" << Dis[0] / loops << " percent 20:" << Dis[1] / loops << " percent 40:"
                 << Dis[2] / loops << " percent 80:" << Dis[3] / loops << " percent 100:" << Dis[4] / loops
                 << " percent total:" << Dis[5] / loops << endl;
            fCD << Dis[0] / loops << "  " << Dis[1] / loops << "  " << Dis[2] / loops << "  " << Dis[3] / loops
                     << "  "
                     << Dis[4] / loops << "  " << Dis[5] / loops << endl;
            cout << "Priority (CD)-------------" << endl;
            cout << "percent 10:" << Pri[0] / loops << " percent 20:" << Pri[1] / loops << " percent 40:"
                 << Pri[2] / loops << " percent 80:" << Pri[3] / loops << " percent 100:" << Pri[4] / loops
                 << " percent total:" << Pri[5] / loops << endl;
            fCD << Pri[0] / loops << "  " << Pri[1] / loops << "  " << Pri[2] / loops << "  " << Pri[3] / loops
                << "  "
                << Pri[4] / loops << "  " << Pri[5] / loops << endl;
            cout << "Frequency (CD)-------------" << endl;
            cout << "percent 10:" << Fre[0] / loops << " percent 20:" << Fre[1] / loops << " percent 40:"
                 << Fre[2] / loops << " percent 80:" << Fre[3] / loops << " percent 100:" << Fre[4] / loops
                 << " percent total:" << Fre[5] / loops << endl;
            fCD << Fre[0] / loops << "  " << Fre[1] / loops << "  " << Fre[2] / loops << "  " << Fre[3] / loops
                << "  "
                << Fre[4] / loops << "  " << Fre[5] / loops << endl;
        } else if (choose == "epsilon") {                              // 输出epsilon数组
            vector<double> Dis(epsilon_vector.size(), 0);  // 大小5，初始值0
            vector<double> Pri(epsilon_vector.size(), 0);
            vector<double> Fre(epsilon_vector.size(), 0);
            for (auto item: CD_total) {
                for (int e = 0; e < epsilon_vector.size(); e++) {
                    Dis[e] += item[e][0];
                    Pri[e] += item[e][1];
                    Fre[e] += item[e][2];
                }
            }
            cout << "Distance (CD)-------------" << endl;
            cout << "epsilon 1:" << Dis[0] / loops << " epsilon 3:" << Dis[1] / loops << " epsilon 5:" << Dis[2] / loops
                 << " epsilon 7:" << Dis[3] / loops << " epsilon 9:" << Dis[4] / loops << endl;
            fCD << Dis[0] / loops << "  " << Dis[1] / loops << "  " << Dis[2] / loops << "  " << Dis[3] / loops
                << "  " << Dis[4] / loops << endl;
            cout << "Priority (CD)-------------" << endl;
            cout << "epsilon 1:" << Pri[0] / loops << " epsilon 3:" << Pri[1] / loops << " epsilon 5:" << Pri[2] / loops
                 << " epsilon 7:" << Pri[3] / loops << " epsilon 9:" << Pri[4] / loops << endl;
            fCD << Pri[0] / loops << "  " << Pri[1] / loops << "  " << Pri[2] / loops << "  " << Pri[3] / loops
                     << "  " << Pri[4] / loops << endl;
            cout << "Frequency (CD)-------------" << endl;
            cout << "epsilon 1:" << Fre[0] / loops << " epsilon 3:" << Fre[1] / loops << " epsilon 5:" << Fre[2] / loops
                 << " epsilon 7:" << Fre[3] / loops << " epsilon 9:" << Fre[4] / loops << endl;
            fCD << Fre[0] / loops << "  " << Fre[1] / loops << "  " << Fre[2] / loops << "  " << Fre[3] / loops
                << "  " << Fre[4] / loops << endl;

        } else {
            vector<double> Dis(k_vector.size(), 0);  // 大小5，初始值0
            vector<double> Pri(k_vector.size(), 0);
            vector<double> Fre(k_vector.size(), 0);
            for (auto item: CD_total) {
                for (int k = 0; k < k_vector.size(); k++) {
                    Dis[k] += item[k][0];
                    Pri[k] += item[k][1];
                    Fre[k] += item[k][2];
                }
            }
            cout << "Distance (CD)-------------" << endl;
            cout << "k 10:" << Dis[0] / loops << " k 15:" << Dis[1] / loops << " k 20:" << Dis[2] / loops << " k 25:"
                 << Dis[3] / loops << " k 30:" << Dis[4] / loops << endl;
            fCD << Dis[0] / loops << "  " << Dis[1] / loops << "  " << Dis[2] / loops << "  " << Dis[3] / loops
                << "  " << Dis[4] / loops << endl;
            cout << "Priority (CD)-------------" << endl;
            cout << "k 10:" << Pri[0] / loops << " k 15:" << Pri[1] / loops << " k 20:" << Pri[2] / loops << " k 25:"
                 << Pri[3] / loops << " k 30:" << Pri[4] / loops << endl;
            fCD << Pri[0] / loops << "  " << Pri[1] / loops << "  " << Pri[2] / loops << "  " << Pri[3] / loops
                     << "  " << Pri[4] / loops << endl;
            cout << "Frequency (CD)-------------" << endl;
            cout << "k 10:" << Fre[0] / loops << " k 15:" << Fre[1] / loops << " k 20:" << Fre[2] / loops << " k 25:"
                 << Fre[3] / loops << " k 30:" << Fre[4] / loops << endl;
            fCD << Fre[0] / loops << "  " << Fre[1] / loops << "  " << Fre[2] / loops << "  " << Fre[3] / loops
                << "  " << Fre[4] / loops << endl;
        }
    }
//    fCD << "maxSkylineNum:" << maxSkylineNum << endl;
    fCD.close();
}




