//
// Created by Charles Du on 4/29/21.
//

#include "Total_Lifted_Content.h"

#include <utility>

#include <Eigen/Eigenvalues>

void Total_Lifted_Content::initialize(const Eigen::MatrixXd &rest_vertices, Eigen::Matrix3Xi faces,
                                           const std::string &form, double alpha)
{
    //
    F = faces;

    // compute restD
    if (form == "harmonic") {
        compute_squared_edge_Length(rest_vertices, F, restD);
        restD *= alpha;
    } else // Tutte form
    {
        restD = Eigen::MatrixXd::Constant(3, F.cols(), alpha);
    }
}



double Total_Lifted_Content::compute_total_lifted_content(const Matrix2Xd &vertices) const {

    VectorXd energyList(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (auto i = 0; i < F.cols(); ++i) {
        Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, i));
        vert.col(1) = vertices.col(F(1, i));
        vert.col(2) = vertices.col(F(2, i));

        const Vector3d &r = restD.col(i);

        energyList(i) = compute_lifted_TriArea(vert, r);
    }

    return energyList.sum();
}

double Total_Lifted_Content::compute_lifted_TriArea(const Matrix2Xd &vert, const Vector3d &r) {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;
    double d1 = e1.squaredNorm() + r(0);
    double d2 = e2.squaredNorm() + r(1);
    double d3 = e3.squaredNorm() + r(2);
    return compute_Heron_tri_area(d1, d2, d3);
}

double
Total_Lifted_Content::compute_total_lifted_content_with_gradient(const Matrix2Xd &vertices, Matrix2Xd &grad) const
{
    int vDim = 2;
    double energy = 0.0;
    grad = Matrix2Xd::Zero(2, vertices.cols());

    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Vector3d r = restD.col(i);


        Matrix2Xd g;
        energy += compute_lifted_TriArea_with_gradient(vert,r,g);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);
    }

    return energy;
}

double
Total_Lifted_Content::compute_lifted_TriArea_with_gradient(const Matrix2Xd &vert, const Vector3d &r,
                                                           Matrix2Xd &grad)
{
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;
    double d1 = e1.squaredNorm() + r(0);
    double d2 = e2.squaredNorm() + r(1);
    double d3 = e3.squaredNorm() + r(2);

    //
    double area = compute_Heron_tri_area(d1, d2, d3);

    //
    double g1 = d2 + d3 - d1;
    double g2 = d3 + d1 - d2;
    double g3 = d1 + d2 - d3;

    //
    auto ge1 = g1 * e1;
    auto ge2 = g2 * e2;
    auto ge3 = g3 * e3;

    //note: grad has the same dimension as vert
    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = ge3 - ge2;
    grad.col(1) = ge1 - ge3;
    grad.col(2) = ge2 - ge1;

    double s = 1 / (8 * area);
    grad *= s;

    return area;
}

double
Total_Lifted_Content::compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(const Matrix2Xd &vertices,
                                                                                           const VectorXi &freeI,
                                                                                           const Matrix3Xi  &F_free,
                                                                                           Matrix2Xd &grad,
                                                                                           SpMat &Hess) const {
    int vDim = 2;
    double energy = 0.0;
    grad = Matrix2Xd::Zero(2, vertices.cols());

    std::vector<Eigen::Triplet<double>> tripletList(3 * 3 * vDim * vDim * F.cols());

    // triangle-wise Hessian of signed area
    // this is used later in the PSD projection step
    MatrixXd signedHess(3 * 2, 3 * 2);
    signedHess << 0.0, 0.0, 0.0, 0.5, 0.0, -0.5,
            0.0, 0.0, -0.5, 0.0, 0.5, 0.0,
            0.0, -0.5, 0.0, 0.0, 0.0, 0.5,
            0.5, 0.0, 0.0, 0.0, -0.5, 0.0,
            0.0, 0.5, 0.0, -0.5, 0.0, 0.0,
            -0.5, 0.0, 0.5, 0.0, 0.0, 0.0;
    //

    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Vector3d r = restD.col(i);


        Matrix2Xd g;
        MatrixXd  hess;
        energy += compute_lifted_TriArea_with_gradient_Hessian(vert,r,g, hess);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);


        // subtract Hessian of signed triangle area
        hess -= signedHess;

        //project hess to PSD
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(hess);
        VectorXd eigenVals = eigenSolver.eigenvalues();
        for (auto j = 0; j < eigenVals.size(); ++j) {
            if (eigenVals(j) < 0.0) {
                eigenVals(j) = 0.0;
            }
        }
        MatrixXd eigenVecs = eigenSolver.eigenvectors();
        hess = eigenVecs * (eigenVals.asDiagonal()) * eigenVecs.transpose();
        //end project hess to PSD

        // update Hessian of free vertices
        int current_index = i * 3 * 3 * vDim * vDim;
        Vector3i indices = F_free.col(i);
        for (int j = 0; j < 3; ++j) {
            int idx_j = indices(j);
            for (int k = 0; k < 3; ++k) {
                int idx_k = indices(k);
                if (idx_j != -1 && idx_k != -1) {
                    for (int l = 0; l < vDim; ++l) {
                        for (int n = 0; n < vDim; ++n) {
                            tripletList[current_index] = Eigen::Triplet<double>(idx_j * vDim + l, idx_k * vDim + n,
                                                                hess(j * vDim + l, k * vDim + n));
                            ++current_index;
                        }
                    }
                }
            }
        }

    }

    // add small positive values to the diagonal of Hessian
    for (auto i = 0; i < vDim * freeI.size(); ++i) {
        tripletList.emplace_back(i, i, 1e-8);
    }

    // get Hessian on free vertices
    Hess.resize(vDim * freeI.size(), vDim * freeI.size());
    Hess.setFromTriplets(tripletList.begin(), tripletList.end());

    return energy;
}

double Total_Lifted_Content::compute_lifted_TriArea_with_gradient_Hessian(const Matrix2Xd &vert, const Vector3d &r,
                                                                          Matrix2Xd &grad, MatrixXd &Hess) {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;
    double d1 = e1.squaredNorm() + r(0);
    double d2 = e2.squaredNorm() + r(1);
    double d3 = e3.squaredNorm() + r(2);

    int vDim = v1.size();

    //
    double area = compute_Heron_tri_area(d1, d2, d3);

    //
    double g1 = d2 + d3 - d1;
    double g2 = d3 + d1 - d2;
    double g3 = d1 + d2 - d3;

    //
    auto ge1 = g1 * e1;
    auto ge2 = g2 * e2;
    auto ge3 = g3 * e3;

    //
    auto av1 = ge3 - ge2;
    auto av2 = ge1 - ge3;
    auto av3 = ge2 - ge1;

    //note: grad has the same dimension as vert
    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = av1;
    grad.col(1) = av2;
    grad.col(2) = av3;
    double s = 1 / (8 * area);
    grad *= s;

    // Hess 1: Laplacian
    Matrix3d Lap;
    Lap << g2 + g3, -g3, -g2,
            -g3, g1 + g3, -g1,
            -g2, -g1, g1 + g2;
    Lap *= s;

    // Kronecker product
    MatrixXd Hess1(3 * vDim, 3 * vDim);
    MatrixXd I = MatrixXd::Identity(vDim, vDim);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Hess1.block(i * vDim, j * vDim, vDim, vDim) = Lap(i, j) * I;
        }
    }

    // Hess 2
    MatrixXd E11(vDim, vDim);
    MatrixXd E22(vDim, vDim);
    MatrixXd E33(vDim, vDim);
    MatrixXd E13(vDim, vDim);
    MatrixXd E12(vDim, vDim);
    MatrixXd E231(vDim, vDim);
    MatrixXd E312(vDim, vDim);

    E11 = e1 * e1.transpose();
    E12 = e1 * e2.transpose();
    E13 = e1 * e3.transpose();
    E22 = e2 * e2.transpose();
    E33 = e3 * e3.transpose();
    E231 = (e2 - e3) * e1.transpose();
    E312 = (e3 - e1) * e2.transpose();

    MatrixXd Hess2(3 * vDim, 3 * vDim);
    Hess2.block(0, 0, vDim, vDim) = E11;
    Hess2.block(0, vDim, vDim, vDim) = E13 + E231;
    Hess2.block(0, 2 * vDim, vDim, vDim) = E12 - E231;
    Hess2.block(vDim, vDim, vDim, vDim) = E22;
    Hess2.block(vDim, 2 * vDim, vDim, vDim) = E12.transpose() + E312;
    Hess2.block(2 * vDim, 2 * vDim, vDim, vDim) = E33;

    Hess2.block(vDim, 0, vDim, vDim) = Hess2.block(0, vDim, vDim, vDim).transpose();
    Hess2.block(2 * vDim, 0, vDim, vDim) = Hess2.block(0, 2 * vDim, vDim, vDim).transpose();
    Hess2.block(2 * vDim, vDim, vDim, vDim) = Hess2.block(vDim, 2 * vDim, vDim, vDim).transpose();

    s *= 2; // 1/4area
    Hess2 *= s;

    // Hess 3
    MatrixXd Hess3(3 * vDim, 3 * vDim);
    Hess3.block(0, 0, vDim, vDim) = av1 * av1.transpose();
    Hess3.block(0, vDim, vDim, vDim) = av1 * av2.transpose();
    Hess3.block(0, 2 * vDim, vDim, vDim) = av1 * av3.transpose();
    Hess3.block(vDim, vDim, vDim, vDim) = av2 * av2.transpose();
    Hess3.block(vDim, 2 * vDim, vDim, vDim) = av2 * av3.transpose();
    Hess3.block(2 * vDim, 2 * vDim, vDim, vDim) = av3 * av3.transpose();

    Hess3.block(vDim, 0, vDim, vDim) = Hess3.block(0, vDim, vDim, vDim).transpose();
    Hess3.block(2 * vDim, 0, vDim, vDim) = Hess3.block(0, 2 * vDim, vDim, vDim).transpose();
    Hess3.block(2 * vDim, vDim, vDim, vDim) = Hess3.block(vDim, 2 * vDim, vDim, vDim).transpose();

    s = s * s / (4 * area);  // 1/64area^3
    Hess3 *= s;

    // Hessian
    Hess.resize(3 * vDim, 3 * vDim);
    Hess = Hess1 - Hess2 - Hess3;

    return area;
}

double Total_Lifted_Content::compute_total_lifted_content(const Matrix2Xd &vertices, VectorXd &energyList) const {
    energyList.resize(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (auto i = 0; i < F.cols(); ++i) {
        Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, i));
        vert.col(1) = vertices.col(F(1, i));
        vert.col(2) = vertices.col(F(2, i));

        const Vector3d &r = restD.col(i);

        energyList(i) = compute_lifted_TriArea(vert, r);
    }

    return energyList.sum();
}

double
Total_Lifted_Content::compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(const Matrix2Xd &vertices,
                                                                                           const VectorXi &freeI,
                                                                                           const Matrix3Xi &F_free,
                                                                                           VectorXd &lifted_content_list,
                                                                                           Matrix2Xd &grad,
                                                                                           SpMat &Hess) const {
    int vDim = 2;
    lifted_content_list.resize(F.cols());
    grad = Matrix2Xd::Zero(2, vertices.cols());

    std::vector<Eigen::Triplet<double>> tripletList(3 * 3 * vDim * vDim * F.cols());

    // triangle-wise Hessian of signed area
    // this is used later in the PSD projection step
    MatrixXd signedHess(3 * 2, 3 * 2);
    signedHess << 0.0, 0.0, 0.0, 0.5, 0.0, -0.5,
            0.0, 0.0, -0.5, 0.0, 0.5, 0.0,
            0.0, -0.5, 0.0, 0.0, 0.0, 0.5,
            0.5, 0.0, 0.0, 0.0, -0.5, 0.0,
            0.0, 0.5, 0.0, -0.5, 0.0, 0.0,
            -0.5, 0.0, 0.5, 0.0, 0.0, 0.0;
    //
#pragma omp parallel
#pragma omp for
    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Vector3d r = restD.col(i);


        Matrix2Xd g;
        MatrixXd  hess;
        lifted_content_list(i) = compute_lifted_TriArea_with_gradient_Hessian(vert,r,g, hess);

#pragma omp critical
        {
            grad.col(i1) += g.col(0);
            grad.col(i2) += g.col(1);
            grad.col(i3) += g.col(2);
        }
        

        // subtract Hessian of signed triangle area
        hess -= signedHess;

        //project hess to PSD
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(hess);
        VectorXd eigenVals = eigenSolver.eigenvalues();
        for (auto j = 0; j < eigenVals.size(); ++j) {
            if (eigenVals(j) < 0.0) {
                eigenVals(j) = 0.0;
            }
        }
        MatrixXd eigenVecs = eigenSolver.eigenvectors();
        hess = eigenVecs * (eigenVals.asDiagonal()) * eigenVecs.transpose();
        //end project hess to PSD

        // update Hessian of free vertices
        int current_index = i * 3 * 3 * vDim * vDim;
        Vector3i indices = F_free.col(i);
        for (int j = 0; j < 3; ++j) {
            int idx_j = indices(j);
            for (int k = 0; k < 3; ++k) {
                int idx_k = indices(k);
                if (idx_j != -1 && idx_k != -1) {
                    for (int l = 0; l < vDim; ++l) {
                        for (int n = 0; n < vDim; ++n) {
                            tripletList[current_index] = Eigen::Triplet<double>(idx_j * vDim + l, idx_k * vDim + n,
                                                                                hess(j * vDim + l, k * vDim + n));
                            ++current_index;
                        }
                    }
                }
            }
        }

    }

    // add small positive values to the diagonal of Hessian
    for (auto i = 0; i < vDim * freeI.size(); ++i) {
        tripletList.emplace_back(i, i, 1e-8);
    }

    // get Hessian on free vertices
    Hess.resize(vDim * freeI.size(), vDim * freeI.size());
    Hess.setFromTriplets(tripletList.begin(), tripletList.end());

    return lifted_content_list.sum();
}


