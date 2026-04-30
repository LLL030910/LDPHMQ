//
// Created by Administrator on 25-9-7.
//
#ifndef SKYLINE_DIS_H
#define SKYLINE_DIS_H

#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include "point.h"

using namespace std;

// ==================== 优化1: 预计算距离矩阵 ====================
//class DistanceCache {
//private:
//    vector<vector<double>> distMatrix;
//    int N;
//
//public:
//    DistanceCache(const vector<Point>& points);
//    double getDistance(int i, int j) const;
//};
//
//// ==================== 优化2: 优化Skyline计算 ====================
//vector<int> computeSkylineSFS(const vector<Point>& points, const vector<int>& dimensions);
//
//// ==================== 优化3: 优化CenterDistance ====================
//class CenterDistanceCache {
//private:
//    const vector<Point>& points;
//    const DistanceCache& distCache;
//    vector<double> minDistToCenters;
//    vector<int> nearestCenter;
//    vector<bool> isCenter;
//    int N;
//
//public:
//    CenterDistanceCache(const vector<Point>& pts, const DistanceCache& cache);
//    void addCenter(int centerIdx);
//    double getMinDistance(int idx) const;
//    int getNearestCenter(int idx) const;
//    int findFarthestPoint() const;
//};


//int pointcmp(Point& a, Point& b);
//double L2Distance(Point& p, Point& q);
//double CenterDistance(int index, vector<Point>& sol, Point& p);
//std::vector<int> computeSkylineDistance(vector<Point>& p, int k, vector<int>& skyline_indices);

double L2DistanceSquared(const Point& a, const Point& b);
double distanceToNearestCenter(int point_idx, const vector<int>& selected_indices, const vector<Point>& points);
vector<int> computeSkylineDistance(const vector<Point>& points, int k, const vector<int>& skyline_indices);

#endif //SKYLINE_DIS_H