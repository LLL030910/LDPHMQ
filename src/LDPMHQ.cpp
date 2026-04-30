#include <omp.h>
#include "LDPMHQ.h"
#include "rand_util.h"

LDPHMQ::LDPHMQ(const std::string &dataset_name, const size_t m) : dataset_(Dataset::get_dataset(dataset_name)), utility_function_set_(dataset_.dim_, m)
{
    for (int i = 0; i < dataset_.points_.size(); ++i)
    {
        dataset_points_.push_back(i);
    }

    utility_function_set_.generate_random_function_set(1.0, true);
}

void LDPHMQ::data_to_LDP() {
    vector<double> max(dataset_.dim_, 0);
    for (UtilityFunction &uf : utility_function_set_.utility_functions_)
    {
        uf.calculate_LDP_utility_values(dataset_.LDP_points_, max);
    }
}

void LDPHMQ::LDP_to_data() {
    for (UtilityFunction &uf : utility_function_set_.utility_functions_)
    {
        uf.calculate_utility_values(dataset_.points_,this->Max_attr);
    }
}

double LDPHMQ::get_realmhr(const vector<std::pair<int, int>> &set) const {
    double mhr_val = 1;
    int n = utility_function_set_.num_utility_function_;
    for (size_t j = 0; j < n; ++j) {
        double happy_ratio, happy_max = 0, happy_tmp;
        for (size_t i = 0; i < set.size(); ++i)
        {
            //happy_tmp = utility_function_set_.utility_functions_[j].utility_function_.dot_product(dataset_.points_[set[i].first]);
            happy_tmp = utility_function_set_.utility_functions_[j].utility_values_[set[i].second];
            happy_max = happy_max > happy_tmp ? happy_max : happy_tmp;
        }
        // std::cout << happy_max << "  ";
        happy_ratio = happy_max / utility_function_set_.utility_functions_[j].real_max;
        if (mhr_val > happy_ratio) {
            mhr_val = happy_ratio;
        }
    }
    // std::cout << std::endl;
    return mhr_val;
}


std::string LDPHMQ::get_name() const
{
    return std::string(dataset_.get_name());
}


