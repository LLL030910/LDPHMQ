#ifndef __DATASET_H__
#define __DATASET_H__

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>

#include "utilities.h"
#include "point.h"

using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::vector;

class Dataset {
public:
    // Forbids copying.
    Dataset(const Dataset&) = delete;
    Dataset& operator=(const Dataset&) = delete;

    Dataset(const std::string& dataset_name);

    static Dataset& get_dataset(const std::string dataset_name)
    {
        static std::unordered_map<std::string, Dataset> name_to_dataset;

        if (!name_to_dataset.count(dataset_name))
        {
            name_to_dataset.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(dataset_name),
                                    std::forward_as_tuple(dataset_name));
        }

        return name_to_dataset.at(dataset_name);
    }


    const vector<Point>& get_universe_points() const;

    // Returns the dataset name.
    const std::string& get_name() const;

    const std::string name_;
    size_t dim_;
    size_t D_size_;
    // size_t att_num_;
    int64_t num_points_;
    vector<Point> points_;
    vector<Point> LDP_points_; //添加差分隐私后的数据集
    vector<Point> Sky_points_; //根据skyline选择的数据集

};

#endif /* __DATASET_H__ */
