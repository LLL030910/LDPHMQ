//
// Created by Administrator on 25-7-1.
//

#ifndef LDP_TLAPULAS_H
#define LDP_TLAPULAS_H

#include <iostream>
#include <stdexcept>
#include <random>

class TLaplace {
private:
    std::vector<double> epsilons;  // 存储多个 epsilon 值
    std::vector<double> sensitivity;
    std::vector<double> lap_scales;  // 预计算的尺度参数（sensitivity / epsilon）
    std::vector<double> Cmax;

    // 随机数生成器（避免重复初始化）
    std::random_device rd;
    std::mt19937 gen;
    std::exponential_distribution<> exp_dist;  // 创建指数分布对象，默认lambda=1.0 f(x) = λ * e^(-λx)

    // 参数检查
    void check_params(double eps, double del) {
        if (eps <= 0) throw std::invalid_argument("Epsilon 必须为正数");
        if (del < 0 || del >= 1) throw std::invalid_argument("Delta 必须在 [0, 1) 范围内");
    }

public:
    // 构造函数（支持单个或多个 epsilon）
    TLaplace(const std::vector<double>& epsilons,const std::vector<double>& Max_ati, double delta = 0.0)
        : epsilons(epsilons), gen(rd()), exp_dist(1.0),Cmax(Max_ati), sensitivity(Max_ati) {

        // 检查每个 epsilon 和 delta
        for (int i = 0; i < Max_ati.size(); i++) {
            check_params(epsilons[i], delta);
            lap_scales.push_back(sensitivity[i] / epsilons[i]);  // 预计算尺度参数
        }
    }

    // 单值加噪（指定使用第 idx 个 epsilon）
    double randomize(double value, size_t idx = 0) {
        if (!std::isfinite(value)) throw std::invalid_argument("输入值必须为有限数");
        if (idx >= lap_scales.size()) throw std::out_of_range("Epsilon 索引越界");

        // 生成拉普拉斯噪声：u1 - u2 ~ Laplace(0, 1)
        double u1 = exp_dist(gen);
        double u2 = exp_dist(gen);
        double noise = (u1 - u2) * lap_scales[idx];  // 缩放至 Laplace(0, sensitivity/epsilon)

        return value + noise;
    }

    // 批量加噪（为多个值使用同一个 epsilon）
    std::vector<double> randomize_batch(const std::vector<double>& values) {

        std::vector<double> noisy_values;
        noisy_values.reserve(values.size());
        int dim = values.size();

        for (int i = 0; i < dim; i++) {
            if (!std::isfinite(values[i])) throw std::invalid_argument("输入值必须为有限数");
            double u1 = exp_dist(gen);
            double u2 = exp_dist(gen);
            double noise = (u1 - u2) * lap_scales[i];
            double noise2 = (values[i] + noise)/Cmax[i];
            if (noise2 > 1)
            {
                noise2 = 1 - generateTinyRandomDouble();
                // noise2 = 1;
            }
            else if (noise2 < 0)
            {
                noise2 = 0;
            }
            else
            {
                //do nothing
            }
            noisy_values.push_back(noise2);
        }

        return noisy_values;
    }

    // 获取隐私参数
    const std::vector<double>& get_epsilons() const { return epsilons; }

    double generateTinyRandomDouble() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<double> dist(1e-15, 1e-9);

        return dist(gen);
    }
};


#endif //LDP_TLAPULAS_H
