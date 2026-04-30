//
// Created by Administrator on 2026/4/9.
//
#include "Skyline_dis.h"
#include <random>
#include <chrono>

using namespace std;

//// ==================== DistanceCache 实现 ====================
//DistanceCache::DistanceCache(const vector<Point>& points) : N(points.size()) {
//    // 预计算所有点对之间的距离（只计算一次）
//    distMatrix.resize(N, vector<double>(N, 0.0));
//
//#pragma omp parallel for
//    for (int i = 0; i < N; i++) {
//        for (int j = i + 1; j < N; j++) {
//            double dist = 0;
//            for (size_t d = 0; d < points[i].coordinates.size(); d++) {
//                double diff = points[i].coordinates[d] - points[j].coordinates[d];
//                dist += diff * diff;
//            }
//            distMatrix[i][j] = dist;
//            distMatrix[j][i] = dist;
//        }
//    }
//}

//double DistanceCache::getDistance(int i, int j) const {
//    return distMatrix[i][j];
//}
//
//// ==================== CenterDistanceCache 实现 ====================
//CenterDistanceCache::CenterDistanceCache(const vector<Point>& pts, const DistanceCache& cache)
//        : points(pts), distCache(cache), N(pts.size()) {
//    minDistToCenters.assign(N, numeric_limits<double>::max());
//    nearestCenter.assign(N, -1);
//    isCenter.assign(N, false);
//}
//
//void CenterDistanceCache::addCenter(int centerIdx) {
//    isCenter[centerIdx] = true;
//    minDistToCenters[centerIdx] = 0;
//    nearestCenter[centerIdx] = centerIdx;
//
//    // 更新所有点到新中心的距离
//#pragma omp parallel for
//    for (int i = 0; i < N; i++) {
//        if (!isCenter[i]) {
//            double dist = distCache.getDistance(i, centerIdx);
//            if (dist < minDistToCenters[i]) {
//                minDistToCenters[i] = dist;
//                nearestCenter[i] = centerIdx;
//            }
//        }
//    }
//}
//
//double CenterDistanceCache::getMinDistance(int idx) const {
//    return minDistToCenters[idx];
//}
//
//int CenterDistanceCache::getNearestCenter(int idx) const {
//    return nearestCenter[idx];
//}
//
//int CenterDistanceCache::findFarthestPoint() const {
//    int farthestIdx = -1;
//    double maxDist = -1;
//    for (int i = 0; i < N; i++) {
//        if (!isCenter[i] && minDistToCenters[i] > maxDist) {
//            maxDist = minDistToCenters[i];
//            farthestIdx = i;
//        }
//    }
//    return farthestIdx;
//}

// ==================== computeSkylineSFS 实现 ====================
//vector<int> computeSkylineSFS(const vector<Point>& points, const vector<int>& dimensions) {
//    int n = points.size();
//    if (n == 0) return {};
//
//    int d = dimensions.size();
//    if (d == 1) {
//        // 单维度：直接找最大值
//        int dim = dimensions[0];
//        double maxVal = -1e100;
//        vector<int> result;
//        for (int i = 0; i < n; i++) {
//            double val = points[i].coordinates[dim];
//            if (val > maxVal + 1e-9) {
//                maxVal = val;
//                result.clear();
//                result.push_back(i);
//            } else if (abs(val - maxVal) < 1e-9) {
//                result.push_back(i);
//            }
//        }
//        return result;
//    }
//
//    // 多维度：先按第一维排序
//    vector<int> indices(n);
//    iota(indices.begin(), indices.end(), 0);
//
//    sort(indices.begin(), indices.end(), [&](int a, int b) {
//        return points[a].coordinates[dimensions[0]] > points[b].coordinates[dimensions[0]];
//    });
//
//    vector<int> skyline;
//    for (int idx : indices) {
//        bool dominated = false;
//        for (int s : skyline) {
//            bool isDominated = true;
//            bool strictBetter = false;
//            for (int dim : dimensions) {
//                if (points[s].coordinates[dim] < points[idx].coordinates[dim]) {
//                    isDominated = false;
//                    break;
//                } else if (points[s].coordinates[dim] > points[idx].coordinates[dim]) {
//                    strictBetter = true;
//                }
//            }
//            if (isDominated && strictBetter) {
//                dominated = true;
//                break;
//            }
//        }
//        if (!dominated) {
//            // 移除被当前点支配的点
//            vector<int> newSkyline;
//            for (int s : skyline) {
//                bool isDominated = true;
//                bool strictBetter = false;
//                for (int dim : dimensions) {
//                    if (points[idx].coordinates[dim] < points[s].coordinates[dim]) {
//                        isDominated = false;
//                        break;
//                    } else if (points[idx].coordinates[dim] > points[s].coordinates[dim]) {
//                        strictBetter = true;
//                    }
//                }
//                if (!(isDominated && strictBetter)) {
//                    newSkyline.push_back(s);
//                }
//            }
//            newSkyline.push_back(idx);
//            skyline = move(newSkyline);
//        }
//    }
//    return skyline;
//}
//
//// ==================== computeSkylineDistance2_Optimized 实现 ====================
//vector<int> computeSkylineDistance2_Optimized(vector<Point>& p, int k) {
//    int N = p.size();
//    if (k < 1 || N == 0) return {};
//
//    vector<int> indices;
//    indices.reserve(k);
//
//    // 预计算距离矩阵（空间换时间）
//    DistanceCache distCache(p);
//
//    // 优化初始点：使用Skyline中的最佳点
//    vector<int> skyline = computeSkylineSFS(p, {0});
//    if (skyline.empty()) return {};
//
//    // 选择Skyline中离原点最远的点作为初始点
//    int startIdx = skyline[0];
//    double maxNorm = -1;
//    for (int idx : skyline) {
//        double norm = 0;
//        for (double coord : p[idx].coordinates) {
//            norm += coord * coord;
//        }
//        if (norm > maxNorm) {
//            maxNorm = norm;
//            startIdx = idx;
//        }
//    }
//
//    // 初始化距离缓存
//    CenterDistanceCache centerCache(p, distCache);
//    centerCache.addCenter(startIdx);
//    indices.push_back(p[startIdx].id);
//
//    // 使用优先队列加速最远点查找
//    auto cmp = [](const pair<double, int>& a, const pair<double, int>& b) {
//        return a.first < b.first;
//    };
//    priority_queue<pair<double, int>, vector<pair<double, int>>, decltype(cmp)> pq(cmp);
//
//    // 初始化堆
//    for (int i = 0; i < N; i++) {
//        if (i != startIdx) {
//            pq.push({centerCache.getMinDistance(i), i});
//        }
//    }
//
//    const double epsilon = 1e-7;
//
//    // 贪心选择
//    while (indices.size() < (size_t)k && !pq.empty()) {
//        auto [dist, idx] = pq.top();
//        pq.pop();
//
//        // 验证距离是否仍然有效
//        if (abs(dist - centerCache.getMinDistance(idx)) > epsilon) continue;
//        if (dist < epsilon) break;
//
//        // 添加新中心点
//        centerCache.addCenter(idx);
//        indices.push_back(p[idx].id);
//    }
//
//    // 确保返回k个点（如果不够，用第一个点填充）
//    while (indices.size() < (size_t)k) {
//        indices.push_back(indices[0]);
//    }
//
//    return indices;
//}



//double L2Distance(Point& p,  Point& q)
////Calculates the L_2 distance between two points
//{
//    double sum=0;
//    int i;
//    size_t d = p.dimension;
//    for (i=0; i<d; ++i)
//        sum+=(p.coordinates[i]-q.coordinates[i])*(p.coordinates[i]-q.coordinates[i]);
//    return sum;  // no need for sqrt
//}


//double CenterDistance(int index, vector<Point>& sol, Point& p)
////Calculates the distance from a Point p to the closest center in sol
//{
//    double mindist = 10000;
//    int i;
//    for (i = 0; i < index; i++) {
//        double dist = L2Distance(sol[i],p);
//        if (dist < mindist)
//            mindist = dist;
//    }
//    return mindist;
//}

double L2DistanceSquared(const Point& a, const Point& b)
{
    double sum = 0.0;
    size_t d = a.dimension;

    for (size_t i = 0; i < d; ++i) {
        double diff = a.coordinates[i] - b.coordinates[i];
        sum += diff * diff;
    }

    return sum;
}

double distanceToNearestCenter(int point_idx, const vector<int>& selected_indices, const vector<Point>& points)
{
    double min_dist = numeric_limits<double>::infinity();

    for (int center_idx : selected_indices) {
        double dist = L2DistanceSquared(points[point_idx], points[center_idx]);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }

    return min_dist;
}

vector<int> computeSkylineDistance(const vector<Point>& points, int k, const vector<int>& skyline_indices ) {
    vector<int> selected_ids;

    if (k <= 0 || skyline_indices.empty()) {
        return selected_ids;
    }

    int m = static_cast<int>(skyline_indices.size());
    k = min(k, m);

    // selected_pos[t] 表示 skyline_indices[t] 是否已经被选入 DS
    vector<char> selected_pos(m, 0);

    // nearest_dist[t] 表示 skyline_indices[t] 到当前 DS 的最近距离
    vector<double> nearest_dist(
            m,
            numeric_limits<double>::infinity()
    );

    // selected_indices 存储的是 points 中的数组下标
    vector<int> selected_indices;
    selected_indices.reserve(k);
    selected_ids.reserve(k);

    // 1. 初始点：选择 skyline_indices[0]
    int first_pos = 0;
    int first_idx = skyline_indices[first_pos];

    selected_pos[first_pos] = 1;
    selected_indices.push_back(first_idx);

    int last_center_idx = first_idx;

    // 2. 每轮只用新加入的中心更新 nearest_dist
    while (static_cast<int>(selected_indices.size()) < k) {

        // Step 1: 更新每个 skyline 点到 DS 的最近距离
        for (int t = 0; t < m; ++t) {
            if (selected_pos[t]) {
                nearest_dist[t] = 0.0;
                continue;
            }

            int point_idx = skyline_indices[t];

            double dist = L2DistanceSquared(
                    points[point_idx],
                    points[last_center_idx]
            );

            if (dist < nearest_dist[t]) {
                nearest_dist[t] = dist;
            }
        }

        // Step 2: 找到当前离 DS 最远的 skyline 点
        double max_dist = -1.0;
        int farthest_pos = -1;

        for (int t = 0; t < m; ++t) {
            if (!selected_pos[t] && nearest_dist[t] > max_dist) {
                max_dist = nearest_dist[t];
                farthest_pos = t;
            }
        }

        // 所有剩余点都与已选点重合，提前停止
        if (farthest_pos == -1 || max_dist <= 1e-12) {
            break;
        }

        selected_pos[farthest_pos] = 1;

        last_center_idx = skyline_indices[farthest_pos];
        selected_indices.push_back(last_center_idx);
    }

    // 3. index 转换为 id
    for (int idx : selected_indices) {
        selected_ids.push_back(static_cast<int>(points[idx].id));
    }

    return selected_ids;
}