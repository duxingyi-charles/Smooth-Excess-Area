//
// Created by Charles Du on 5/3/21.
//

#ifndef TLC_OPTIMIZATION_UTIL_H
#define TLC_OPTIMIZATION_UTIL_H

#include <Eigen/Core>

bool importData(const char* filename,
                Eigen::MatrixXd &restV,
                Eigen::Matrix2Xd &initV,
                Eigen::Matrix3Xi &F,
                Eigen::VectorXi &handles
);

#endif //TLC_OPTIMIZATION_UTIL_H
