//
// Created by Administrator on 25-9-4.
//

#ifndef SKYLINE_PRIO_H
#define SKYLINE_PRIO_H

#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <iomanip>
#include "point.h"
// Skyline优先级计算相关函数声明
std::vector<Point> readDataFromPoints_prio(const std::vector<Point>& points, int& n, int& d_num);
std::vector<int> computeSkylineBNL_prio(const std::vector<Point>& points, const std::vector<int>& dimensions);
double calculateSum_prio(const Point& point);
std::vector<int> computeSkylinePriority(const std::vector<Point>& points, double x_percent = 20.0);

#endif //SKYLINE_PRIO_H
