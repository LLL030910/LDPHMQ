#include "Skyline_prior.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <bitset>

#include "point.h"

using namespace std;

// 从Point向量读取数据
vector<Point> readDataFromPoints_prio(const vector<Point>& points, int& n, int& d_num) {
    if (points.empty()) {
        n = 0;
        d_num = 0;
        return vector<Point>();
    }

    n = points.size();
    d_num = points[0].coordinates.size(); // 数值维度数

    return points;
}

// 使用BNL算法计算子空间的Skyline点
vector<int> computeSkylineBNL_prio(const vector<Point>& points, const vector<int>& dimensions) {
    vector<int> skyline;
    int n = points.size();

    if (n == 0) return skyline;

    // 窗口算法
    vector<int> window;
    window.push_back(0); // 第一个点的索引

    if (dimensions.size() == 1) {
        // -------- 单维度处理 --------
        int dim = dimensions[0];
        double max_val = points[0].coordinates[dim];    // int ---> double
        vector<int> skyline_ids = { static_cast<int>(points[0].id) };

        for (int i = 1; i < n; i++) {
            double val = points[i].coordinates[dim];
            if (val > max_val) {
                // 找到更大的值，清空之前的结果
                max_val = val;
                skyline_ids.clear();
                skyline_ids.push_back(points[i].id);
            } else if (val == max_val) {
                // 相等则并列最大
                skyline_ids.push_back(points[i].id);
            }
        }
        window = skyline_ids;  // Skyline 结果
    } else {
        // -------- 多维度处理（保持原逻辑） --------
        for (int i = 1; i < n; i++) {
            bool dominated = false;

            for (int j : window) {
                bool j_dominates_i = true;
                bool strictly_better = false;
                for (int dim : dimensions) {
                    if (points[j].coordinates[dim] < points[i].coordinates[dim]) {
                        j_dominates_i = false;
                        break;
                    } else if (points[j].coordinates[dim] > points[i].coordinates[dim]) {
                        strictly_better = true;
                    }
                }
                if (j_dominates_i && strictly_better) {
                    dominated = true;
                    break; // 在这里 break 是安全的，因为我们只关心 dominated 标志
                }
            }

            // 步骤 2: 只有当 i 没有被支配时，才去构建新的 window
            if (!dominated) {
                vector<int> new_window;
                // 遍历旧 window，淘汰掉被 i 支配的点
                for (int j : window) {
                    bool i_dominates_j = true;
                    bool strictly_better = false;
                    for (int dim : dimensions) {
                        if (points[i].coordinates[dim] < points[j].coordinates[dim]) {
                            i_dominates_j = false;
                            break;
                        } else if (points[i].coordinates[dim] > points[j].coordinates[dim]) {
                            strictly_better = true;
                        }
                    }

                    // 如果 i 不支配 j，那么 j 幸存下来
                    if (!(i_dominates_j && strictly_better)) {
                        new_window.push_back(j);
                    }
                }

                // 把新的 Skyline 点 i 也加进去 (注意：这里要用索引 i，而不是 id)
                new_window.push_back(i);

                // 用最终构建好的 new_window 替换旧的
                window = new_window;
            }
        }
    }

    return window;
}

// 计算点的维度数值之和
double calculateSum_prio(const Point& point) {
    return accumulate(point.coordinates.begin(), point.coordinates.end(), 0.0);
}

// 主计算函数 - 计算Skyline优先级
std::vector<int> computeSkylinePriority(const vector<Point>& points, double x_percent) {
    int n, d_num;

    // 从Point向量读取数据
    vector<Point> data = readDataFromPoints_prio(points, n, d_num);

    if (n == 0) {
        return vector<int>();
    }

    // 初始化优先级数组为总维度数d_num
    // vector<int> priority(n, (1 << d_num) );
    unordered_map<int, int> priority;
    vector<double> point_sums(n);

    for (int i = 0; i < n; i++) {
        point_sums[i] = calculateSum_prio(data[i]);
    }

    // 枚举所有非空子空间 (1 to 2^d_num - 1)
    int total_subspaces = (1 << d_num) - 1;

    // 处理每个子空间
    for (int bitmap = 1; bitmap <= total_subspaces; bitmap++) {
        // 获取当前子空间的维度
        vector<int> dimensions;
        for (int i = 0; i < d_num; i++) {
            if (bitmap & (1 << i)) {
                dimensions.push_back(i);
            }
        }

        // 当前子空间的维度大小
        int current_dimension_count = dimensions.size();

        // 使用BNL算法计算Skyline
        vector<int> skyline_indices = computeSkylineBNL_prio(data, dimensions);

        // 更新优先级
        for (int idx : skyline_indices) {
            if (priority[idx] == 0) priority[idx] = current_dimension_count;
            // 如果当前子空间维度数小于该点当前优先级，则更新优先级
            if (current_dimension_count < priority[idx]) {
                priority[idx] = current_dimension_count;
            }
        }
    }


    // 创建索引数组并排序
    // vector<int> indices(n);
    // iota(indices.begin(), indices.end(), 0);
    vector<int> indices;
    for (const auto& pair : priority) {
        indices.push_back(pair.first); // pair.first 是点id，pair.second 是频率
    }

    // 自定义排序：先按优先级升序（优先级越小越好），再按维度数值之和降序
    sort(indices.begin(), indices.end(), [&](int a, int b) {
        if (priority[a] != priority[b]) {
            return priority[a] < priority[b]; // 优先级升序
        }
        return point_sums[a] > point_sums[b]; // 数值和降序
    });

    // 设置x为指定百分比
    int num_candidates = static_cast<int>(priority.size() * x_percent / 100.0);
    // 取前m个点
    vector<int> result_indices(indices.begin(), indices.begin() + num_candidates);
    return result_indices;
}

// 辅助函数：打印结果
void printResults(const vector<Point>& points, const vector<int>& indices, const vector<int>& priority) {
    cout << "\n==========================================" << endl;
    cout << "SKYLINE PRIORITY RESULTS" << endl;
    cout << "==========================================" << endl;
    cout << "Index\tPriority\tSum\t\tCoordinates" << endl;
    cout << "------------------------------------------" << endl;

    for (int idx : indices) {
        cout << idx << "\t" << priority[idx] << "\t\t"
             << accumulate(points[idx].coordinates.begin(), points[idx].coordinates.end(), 0.0)
             << "\t[";
        for (size_t j = 0; j < points[idx].coordinates.size(); j++) {
            cout << points[idx].coordinates[j];
            if (j < points[idx].coordinates.size() - 1) cout << ", ";
        }
        cout << "]" << endl;
    }
}