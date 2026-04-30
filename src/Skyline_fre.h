#ifndef SKYLINE_FRE_H
#define SKYLINE_FRE_H

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
#include "Skyline_dis.h"
// Skyline频率计算相关函数声明
std::vector<Point> readDataFromPoints(const std::vector<Point>& points, int& n, int& d_num);
std::vector<int> computeSkylineBNL_Fre(const std::vector<Point>& points, const std::vector<int>& dimensions);
double calculateSum(const Point& point);
std::vector<int> computeSkylineFrequency(const std::vector<Point>& points, double x_percent = 20.0);

void computerSkylineBNL(std::vector<Point>& points, int dim, vector<int>& Fre, vector<int>& Pri, vector<int>& Dis, double percent);

#endif // SKYLINE_FRE_H