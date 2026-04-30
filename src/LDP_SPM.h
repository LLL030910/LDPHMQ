//
// Created by Administrator on 25-7-4.
//

#ifndef LDP_SPM_H
#define LDP_SPM_H


#include <iostream>
#include <stdexcept>
#include <random>
#include "point.h"

class SPM {
protected:
    std::vector<double> epsilons;  // 存储多个 epsilon 值
    std::vector<double> Cmax;
    double C;

    double check_epsilon(double epsilon) {
        if (epsilon <= 0) {
            throw std::invalid_argument("epsilon must be positive");
        }
        return epsilon;
    }

    double check_value(double value) {
        if (value <100000 || value > 0) {
            throw std::invalid_argument("the input value must be in domain=[0,Ci]");
        }
        return value;
    }

public:
    SPM(const std::vector<double>& epsilons, const std::vector<double>& Max_ati)
            : epsilons(epsilons),Cmax(Max_ati) {
    }

    std::vector<double> randomize_batch(std::vector<double> values, double minor = 1e-10) { //minor的作用防止p_h为0

        int dim = values.size();
        std::vector<double> result;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);  // 生成[0,1)区间的均匀分布随机数
        for (int i=0; i < dim; i++) {
            double p1 = 1/(exp(epsilons[i]/2)+1)*values[i]/Cmax[i];
            double p2 = exp(epsilons[i]/2)/(exp(epsilons[i]/2)+1);
            double L = p2*values[i]/Cmax[i];
            double R = p2*(values[i]/Cmax[i]+1/exp(epsilons[i]/2));
            double rnd = dis(gen);

            if (rnd < p1)
            {
                std::uniform_real_distribution<> dis_v(0, L);
                result.push_back(dis_v(gen));
            }
            else if (rnd < p1+p2)
            {
                std::uniform_real_distribution<> dis_v(L, R);
                result.push_back(dis_v(gen));
            }
            else
            {
                std::uniform_real_distribution<> dis_v(R, 1);
                result.push_back(dis_v(gen));
            }
        }
        return result;
    }
};


#endif //LDP_SPM_H
