#ifndef __LDPMHQ_H__
#define __LDPMHQ_H__

#include <algorithm>
#include <memory>
#include <vector>
#include <set>
#include <utility>
#include <execution>
#include <iterator>
#include "utilities.h"
#include "dataset.h"
#include "point.h"

using std::set;
using std::vector;

class LDPHMQ {
public:

    LDPHMQ(const std::string& dataset_name, const size_t m);

    double get_realmhr(const vector<std::pair<int,int>>& set) const;


    void data_to_LDP();

    void LDP_to_data();

    // Return the name.
    std::string get_name() const;

    vector<double> Max_attr;
    Dataset& dataset_; // whole dataset, including universe points
    UtilityFunctionSet utility_function_set_; // sampled utility functions

private:

    vector<int> dataset_points_; // indexes of whole dataset

};

#endif /* __LDPMHQ_H__ */