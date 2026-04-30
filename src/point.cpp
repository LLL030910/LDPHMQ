#include <cassert>
#include <cmath>
#include <iostream>
#include "point.h"

#include "rand_util.h"

using std::vector;

Point::Point() {

}

Point::Point(size_t dimension)
{
    this->dimension = dimension;
    for (size_t i = 0; i < dimension; i++)
        coordinates.push_back(0);
}

Point::Point(size_t dimension, size_t id)
{
    this->dimension = dimension;
    this->id=id;
    for(size_t i = 0; i < dimension; i++)
        coordinates.push_back(0);
}

Point::Point(size_t dimension, vector<double> &coordinates)
{
    assert(dimension == coordinates.size());
    this->dimension = dimension;
    this->coordinates = coordinates;
}


Point::Point(size_t id, size_t dimension, vector<double> &coordinates)
{
    assert(dimension == coordinates.size());
    this->dimension = dimension;
    this->coordinates = coordinates;
    this->id = id;
}

Point::Point(const Point &other_point)
{
    this->dimension = other_point.dimension;
    this->coordinates = other_point.coordinates;
    this->id = other_point.id;
}

Point Point::operator-(const Point &other_point) const
{
    assert(other_point.dimension == dimension);
    Point diff(dimension);
    for (size_t i = 0; i < dimension; i++)
        diff.coordinates[i] = coordinates[i] - other_point.coordinates[i];
    return diff;
}

Point Point::operator*(double factor) const
{
    Point mult(dimension);
    for (size_t i = 0; i < dimension; i++)
        mult.coordinates[i] = coordinates[i] * factor;
    return mult;
}

Point& Point::operator=(const Point &p)
{
    if(this != &p)
    {
        dimension = p.dimension;
        id = p.id;
        coordinates = p.coordinates;
    }
    return *this;
}

void Point::set(size_t dimension, vector<double> &coordinates)
{
    assert(dimension == coordinates.size());
    this->dimension = dimension;
    this->coordinates = coordinates;
}

double Point::distance_to(const Point &other_point) const
{
    Point diff = (*this) - other_point;
    return diff.length();
}

double Point::length() const
{
    double length = 0;
    for (size_t i = 0; i < dimension; i++)
        length += coordinates[i] * coordinates[i];
    return sqrt(length);
}

double Point::dotP(const Point& other_point) const
{
    assert(dimension == other_point.dimension);


    double dotp = 0.0;
    for(size_t i = 0; i < dimension; i++)
        dotp += coordinates[i] * other_point.coordinates[i];

    return dotp;
}

double Point::dotP(vector<double> coord) const
{
    double dotp = 0.0;
    for(size_t i = 0; i < dimension; i++)
        dotp += coordinates[i] * coord[i];

    return dotp;
}

double Point::dot_product(const Point &other_point, vector<double> &Max_attr) const
{
    assert(dimension == other_point.dimension);
    double dotp = 0.0;
    for (size_t i = 0; i < dimension; i++)
    {
        if (Max_attr[i] <= 1)
        {
            dotp += coordinates[i] * other_point.coordinates[i];
        }
        else
        {
            dotp += coordinates[i] * other_point.coordinates[i] / Max_attr[i];
        }
    }
    return dotp;
}

double Point::get_coordinate(size_t idx) const
{
    assert(idx < dimension);
    return coordinates[idx];
}

bool Point::dominates(const Point &other_point) const
{
    assert(dimension == other_point.dimension);
    bool at_least_one = false;
    for (size_t i = 0; i < dimension; i++)
    {
        if (coordinates[i] < other_point.coordinates[i])
        {
            return false;
        }
        if (coordinates[i] > other_point.coordinates[i])
        {
            at_least_one = true;
        }
    }
    return at_least_one;
}


Point Point::abs(const Point &other_point)  //函数用于计算并返回另一个 Point 对象的坐标的绝对值。
{
    Point p(other_point.dimension, other_point.id);
    for (size_t i = 0; i < other_point.dimension; i++)
    {
        p.coordinates[i] = std::abs(other_point.coordinates[i]);
    }
    return p;
}

boost::numeric::ublas::vector<double> Point::to_ublas(const Point& p)
{
    boost::numeric::ublas::vector<double> ublasp(p.get_dimension());
    for(size_t i = 0; i < p.get_dimension(); i++)
        ublasp[i] = p.get_coordinate(i);

    return ublasp;
}

Point Point::from_ublas(const boost::numeric::ublas::vector<double>& ublasp)
{
    std::vector<double> coords;
    for(size_t i = 0; i < ublasp.size(); i++)
        coords.push_back(ublasp[i]);

    return Point(ublasp.size(), coords);
}

bool Point::operator==(const Point &other_point) const
{
    return dimension == other_point.dimension && coordinates == other_point.coordinates;
}

void Point::scale_to_length(double len)
{
    assert(len >= 0);
    double my_len = length();
    assert(len == 0 || my_len > 0);
    double factor = (my_len > 0) ? (len / my_len) : 0;
    for (size_t i = 0; i < dimension; i++)
    {
        coordinates[i] *= factor;
        coordinates[i] *= coordinates[i];
    }
}

UtilityFunctionSet::UtilityFunctionSet()
{
}

UtilityFunctionSet::UtilityFunctionSet(size_t dim, size_t num_utility_function) : dim_(dim), num_utility_function_(num_utility_function)
{
    utility_functions_.reserve(num_utility_function_);
}

void UtilityFunctionSet::generate_random_function_set(double sphere_radius, bool first_orthant)
{
    //随机生成效用函数

    for (size_t i = 0; i < num_utility_function_; ++i)
    {
        Point p(dim_, i);
        RandUtil::get_random_direction(dim_, p);
        if (first_orthant)
        {
            p = Point::abs(p);
        }
        p.scale_to_length(sphere_radius);
        //TODO: check (in theory will call constructior UtilityFunction(const Point &e))
        utility_functions_.emplace_back(p);
    }

    std::cout << "Already sampled " << num_utility_function_ << " utility functions." << std::endl;
}

std::ostream &operator<<(std::ostream &s, const Point &p)
{
    s << "Dimension: " << p.dimension << ", Coordinates: [";
    for (size_t i = 0; i < p.coordinates.size(); ++i)
    {
        s << p.coordinates[i];
        if (i < p.coordinates.size() - 1)
        {
            s << ", ";
        }
    }
    s << "]";
    return s;
}

UtilityFunction::UtilityFunction()
{
}

UtilityFunction::UtilityFunction(const Point &e) : utility_function_(e)
{
    utility_values_.clear();
    // rank_table_.clear();
    // sorted_utility_values_.clear();
}

UtilityFunction::UtilityFunction(const UtilityFunction &uf) : utility_function_(uf.utility_function_), utility_values_(uf.utility_values_),real_max(uf.real_max), ldp_max(uf.ldp_max)
{
}


void UtilityFunction:: calculate_LDP_utility_values(const vector<Point> &data, vector<double> &Max_attr) {
    ldp_utility_values_.clear();
    ldp_utility_values_.reserve(data.size());

    double max = 0;
    vector<pair<double, int>> utility_values;  //第一个是函数值，第二个是元素
    for(int i = 0; i < data.size(); ++i) {
        double uv = utility_function_.dot_product(data[i], Max_attr);
        if (uv > max) {
            max = uv;
        }
        //保证尺度不变性，判断Max_attr是否大于1

        ldp_utility_values_.push_back(uv);
        utility_values.emplace_back(uv, i);
    }
    ldp_max = max;
}

void UtilityFunction:: calculate_utility_values(const vector<Point> &data, vector<double> &Max_attr) {
    utility_values_.clear();
    utility_values_.reserve(data.size());

    double max = 0;
    vector<pair<double, int>> utility_values;  //第一个是函数值，第二个是元素
    for(int i = 0; i < data.size(); ++i) {
        double uv = utility_function_.dot_product(data[i], Max_attr);
        if (uv > max) {
            max = uv;
        }
        //保证尺度不变性，判断Max_attr是否大于1

        utility_values_.push_back(uv);
        utility_values.emplace_back(uv, i);
    }
    real_max = max;
}

