//
// Created by Charles Du on 5/3/21.
//

#include "optimization_util.h"

#include <iostream>
#include <fstream>

bool importData(const char* filename,
                Eigen::MatrixXd &restV,
                Eigen::Matrix2Xd &initV,
                Eigen::Matrix3Xi &F,
                Eigen::VectorXi &handles
)
{
    std::ifstream in_file(filename);

    if (! in_file.is_open()) {
        std::cerr << "Failed to open " << filename << "!" << std::endl;
        return false;
    }

    //read the file
    size_t n, ndim;
    // restV
    in_file >> n >> ndim;
    restV.resize(ndim, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < ndim; ++j) {
            in_file >> restV(j,i);
        }
    }

    //initV
    in_file >> n >> ndim;
    initV.resize(ndim,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < ndim; ++j) {
            in_file >> initV(j,i);
        }

    }

    //F
    size_t simplexSize;
    in_file >> n >> simplexSize;
    F.resize(simplexSize, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < simplexSize; ++j)
        {
            in_file >> F(j,i);
        }
    }

    //handles
    in_file >> n;
    handles.resize(n);
    for (int i = 0; i < n; ++i)
    {
        in_file >> handles(i);
    }

    in_file.close();

    return true;
}