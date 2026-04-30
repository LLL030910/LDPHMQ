//
// Created by 1348942332 on 25-11-27.
//

#ifndef LDP_SW_H
#define LDP_SW_H
#include <random>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

class SquareWaveMechanism {
private:
    std::vector<double> epsilons;    // 隐私预算
    std::vector<double> Cmax;          // 输入数据的最大范围 [0, C]
    std::vector<double> b;           // 波宽参数
    std::vector<double> p;           // 高概率密度 (接近真实值区域)
    std::vector<double> q;           // 低概率密度 (远离真实值区域)
    std::vector<double> high_prob_mass; // 在 [v-b, v+b] 区间内采样的总概率 (面积)

    // 随机数生成器
    std::mt19937 gen;
    std::uniform_real_distribution<double> uniform_dist;

    // 计算最优的 b 值 (根据论文公式 5.3节)
    // 目标是最大化输入和输出之间的互信息
    double calculateOptimalB(double eps) {
        double exp_eps = std::exp(eps);
        // 分子: epsilon * e^eps - e^eps + 1
        double numerator = eps * exp_eps - exp_eps + 1.0;
        // 分母: 2 * e^eps * (e^eps - 1 - epsilon)
        double denominator = 2.0 * exp_eps * (exp_eps - 1.0 - eps);

        // 防止分母为0 (虽然 eps>0 时一般不会发生)
        if (std::abs(denominator) < 1e-9) return 0.5;
        double b = numerator / denominator;
        return b;
    }

public:
    SquareWaveMechanism(std::vector<double> eps, std::vector<double> max_val)
        : epsilons(eps), Cmax(max_val), uniform_dist(0.0, 1.0) {
        b.resize(max_val.size());
        p.resize(max_val.size());
        q.resize(max_val.size());
        high_prob_mass.resize(max_val.size());
        for (int i = 0; i < eps.size(); i++) {
            // 初始化随机数种子
            std::random_device rd;
            gen.seed(rd());

            // 1. 计算参数 b
            b[i] = calculateOptimalB(eps[i]);

            // 2. 计算概率密度 p 和 q
            double exp_eps = std::exp(eps[i]);
            double term = 2.0 * b[i] * exp_eps + 1.0;

            p[i] = exp_eps / term;
            q[i] = 1.0 / term;

            // 3. 计算落在高概率区间 [v-b, v+b] 的总概率质量 (面积 = 宽 * 高)
            // 宽度是 2b, 高度是 p
            high_prob_mass[i] = 2.0 * b[i] * p[i];

        }

        // 调试输出参数
        // std::cout << "--- SW Mechanism Initialized ---" << std::endl;
        // std::cout << "Input Range: [0, " << C << "]" << std::endl;
        // std::cout << "Epsilon: " << epsilon << std::endl;
        // std::cout << "Optimal b: " << b << std::endl;
        // std::cout << "High Prob Mass (p*2b): " << high_prob_mass << std::endl;
        // std::cout << "Low Prob Mass (q*1): " << (1.0 - high_prob_mass) << std::endl;
        // std::cout << "--------------------------------" << std::endl;
    }


    // 核心扰动函数
    std::vector<double> randomize_batch(std::vector<double> values) {
        std::vector<double> result;
        for (int i = 0; i < values.size(); i++) {
            // Step 1: 归一化 (Normalize)
            // 将输入 x 从 [0, C] 映射到 v ∈ [0, 1]
            // 这里做一个边界保护，防止 x 超出 [0, C]
            double v = std::clamp(values[i], 0.0, Cmax[i]) / Cmax[i];

            double noisy_v;

            // 生成一个 [0, 1] 的随机数决定采样区域
            double u = uniform_dist(gen);

            if (u < high_prob_mass[i]) {
                // Case A: 在高概率区间 [v-b, v+b] 内采样
                // 重新归一化随机数 u 到 [0, 1] 用于在区间内均匀分布
                // 或者直接生成一个新的随机数
                double u_local = uniform_dist(gen);

                // 映射到 [v-b, v+b]
                noisy_v = (v - b[i]) + u_local * (2.0 * b[i]);

            } else {
                // Case B: 在低概率区间 [-b, 1+b] \ [v-b, v+b] 内采样
                // 低概率区间的总长度是 1 (因为总长度 1+2b，减去高概率部分 2b)
                // 低概率区间由两部分组成:
                // 左边: [-b, v-b)  长度为 v
                // 右边: (v+b, 1+b] 长度为 1-v

                double u_local = uniform_dist(gen); // 生成 [0, 1]

                if (u_local < v) {
                    // 落在左边部分 (根据长度比例 v)
                    // 将 u_local 映射到 [-b, v-b]
                    // 此时 u_local 已经在 [0, v] 范围内，只需要平移 -b
                    noisy_v = -b[i] + u_local;
                } else {
                    // 落在右边部分
                    // 将 u_local 映射到 [v+b, 1+b]
                    // 偏移量计算: u_local - v 得到 [0, 1-v] 的偏移
                    noisy_v = (v + b[i]) + (u_local - v);
                }
            }

            // Step 3: 标准化
            // 原始 SW 输出范围是 [-b, 1+b]，标准化到[0, 1]
            noisy_v = (noisy_v + b[i]) / (2 * b[i] + 1);
            result.push_back(noisy_v);
        }
        return result;
    }

};
#endif //LDP_SW_H
