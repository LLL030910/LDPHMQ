#ifndef __POINT_H__
#define __POINT_H__

#include <vector>
#include <map>
#include <set>
#include <functional>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <sstream>

#include <vector>
#include <string>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

class Point
{
public:
    // 构造函数
    Point();
    Point(size_t dimension);
    Point(size_t dimension, size_t id);
    Point(size_t dimension, vector<double> &coordinates);
    Point(size_t id, size_t dimension, vector<double> &coordinates);
    Point(const Point &other_point);

    // 基本getter方法
    size_t get_dimension() const { return dimension; }
    size_t getId() const { return id; }
    vector<double> get_coordinates() const { return coordinates; }
    double get_coordinate(size_t idx) const;

    // 几何运算方法
    double distance_to(const Point &other_point) const;
    double length() const;
    double dot_product(const Point &other_point, vector<double> &Max_attr) const;
    double dotP(const Point &other_point) const;
    double dotP(vector<double> coord) const;
    void scale_to_length(double len);

    // 关系运算方法
    bool dominates(const Point &other_point) const;
    bool operator==(const Point &other_point) const;

    // 运算符重载
    Point operator-(const Point &other_point) const;
    Point operator*(const double factor) const;
    Point& operator=(const Point &other_point);

    // 设置方法
    void set(size_t dimension, vector<double> &coordinates);
    void set_id(size_t new_id) { id = new_id; }

    // 静态方法
    static Point abs(const Point &other_point);
    static boost::numeric::ublas::vector<double> to_ublas(const Point& p);
    static Point from_ublas(const boost::numeric::ublas::vector<double>& ublasp);

    // 输出方法
    void print() const;
    friend std::ostream &operator<<(std::ostream &s, const Point &p);

    // 数据成员
    size_t dimension;
    size_t id;
    vector<double> coordinates;
};


class UtilityFunction
{
private:


public:
    Point utility_function_;
    vector<double> utility_values_;  // utility values of whole dataset
    vector<double> ldp_utility_values_;  // 加噪数据集的效用值

    double real_max;    // 存储真实数据集效用值的最大值
    double ldp_max;     // 存储加噪数据集效用值的最大值

    UtilityFunction();
    UtilityFunction(const Point &e);
    UtilityFunction(const UtilityFunction &uf);
    //计算LDP后的函数值
    void calculate_LDP_utility_values(const vector<Point> &data, vector<double> &Max_attr);
    void calculate_utility_values(const vector<Point> &data, vector<double> &Max_attr);
};

class UtilityFunctionSet
{
public:
    UtilityFunctionSet();
    UtilityFunctionSet(size_t dim, size_t num_utility_function);
    void generate_random_function_set(double sphere_radius, bool first_orthant);

    size_t dim_;
    size_t num_utility_function_;
    vector<UtilityFunction> utility_functions_;
};

#endif /* __POINT_H__ */

