# Local Differential Privacy Preservation in Happiness Maximization Queries

This repository contains the source code, datasets, and experimental results for the paper **"Local Differential Privacy Preservation in Happiness Maximization Queries for Data-Centric Consumer Electronics"**.

The project studies how to answer **Happiness Maximization Queries (HMQ)** over multidimensional data while preserving **Local Differential Privacy (LDP)**. It implements three privacy mechanisms for HMQ:

- `T-Laplace` (Truncated Laplace Mechanism)
- `SPM` (Scalable Piecewise Mechanism)
- `SW` (square-wave-based mechanism)

To improve efficiency, the repository also includes skyline-based data reduction strategies:

- full skyline selection
- representative skyline by distance
- representative skyline by priority
- representative skyline by frequency

## Repository Structure

```text
.
|- datasets/              # Input datasets used in experiments
|- results/               # Output results reported by experiments
|- src/                   # C++ source code
```

## Datasets

This repository currently includes the following datasets under `datasets/`:

- `Biometric-wearable`: https://www.kaggle.com/datasets/ziya07/forensic-biometric-wearable-dataset
- `smart home`: https://zenodo.org/records/10925475
- `healthcare-IoT`: https://www.kaggle.com/datasets/ziya07/healthcare-iot-data
- `Anti-Cor_4_10000`: synthetic anti-correlated dataset

## Requirements

- CMake 3.20 or later
- A C++17-compatible compiler
- GLPK : https://ftp.gnu.org/gnu/glpk/
- Boost : https://www.boost.org/

## Code Organization

### Core HMQ and dataset components

- [src/LDPMHQ.h](src/LDPMHQ.h), [src/LDPMHQ.cpp](src/LDPMHQ.cpp): HMQ-related core logic
- [src/dataset.h](src/dataset.h), [src/dataset.cpp](src/dataset.cpp): dataset loading and dataset management
- [src/point.h](src/point.h), [src/point.cpp](src/point.cpp): multidimensional point representation
- [src/utilities.h](src/utilities.h), [src/utilities.cpp](src/utilities.cpp): shared utility functions
- [src/rand_util.h](src/rand_util.h), [src/rand_util.cpp](src/rand_util.cpp): random support utilities

### Local differential privacy mechanisms

- [src/LDP_TLaplace.h](src/LDP_TLaplace.h): truncated Laplace mechanism
- [src/LDP_SPM.h](src/LDP_SPM.h): scalable piecewise mechanism
- [src/LDP_SW.h](src/LDP_SW.h): square-wave-based mechanism

### Skyline-based data reduction

- [src/Skyline_dis.h](src/Skyline_dis.h), [src/Skyline_dis.cpp](src/Skyline_dis.cpp): representative skyline by distance
- [src/Skyline_prior.h](src/Skyline_prior.h), [src/Skyline_prior.cpp](src/Skyline_prior.cpp): representative skyline by priority
- [src/Skyline_fre.h](src/Skyline_fre.h), [src/Skyline_fre.cpp](src/Skyline_fre.cpp): representative skyline by frequency

### Sphere / optimization-related code

- [src/sphere.h](src/sphere.h), [src/sphere.cpp](src/sphere.cpp)
- [src/sphere_lp.h](src/sphere_lp.h), [src/sphere_lp.cpp](src/sphere_lp.cpp)
- [src/sphere_MVE.h](src/sphere_MVE.h), [src/sphere_MVE.cpp](src/sphere_MVE.cpp)
- [src/sphere_operation.h](src/sphere_operation.h), [src/sphere_operation.cpp](src/sphere_operation.cpp)
- [src/sphere_data_utility.h](src/sphere_data_utility.h), [src/sphere_data_utility.cpp](src/sphere_data_utility.cpp)
- [src/sphere_RMSUtils.h](src/sphere_RMSUtils.h), [src/sphere_RMSUtils.cpp](src/sphere_RMSUtils.cpp)

## How to Reproduce the Experiments

1. Build the project with CMake.
2. Open [src/main.cpp](src/main.cpp).
3. Select the dataset by editing the `LDPHMQ func_ldpmhq(...)` line.
4. Set the experiment mode through `choose`: `percent`, `epsilon`, or `k`.
5. Adjust `loops`, `percent_vector`, `epsilon_vector`, `k_vector`, and the output path under `results/`.
6. Run the executable and collect the generated result files.
