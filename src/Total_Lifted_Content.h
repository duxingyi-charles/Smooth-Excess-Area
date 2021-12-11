//
// Created by Charles Du on 4/29/21.
// This class computes Total lifted content energy and its gradient
//

#ifndef TLC_TOTAL_LIFTED_CONTENT_H
#define TLC_TOTAL_LIFTED_CONTENT_H

#include "geo_util.h"

using namespace Eigen;

//#include <Eigen/CholmodSupport>
#include <Eigen/Sparse>
typedef SparseMatrix<double> SpMat;

class Total_Lifted_Content {
public:

    Total_Lifted_Content() = default;

    // initialize from rest mesh and parameters: get ready for computing TLC energy and derivatives
    void initialize(const Eigen::MatrixXd &rest_vertices, Eigen::Matrix3Xi faces,
                         const std::string &form, double alpha);

    ~Total_Lifted_Content() = default;

    // compute total lifted content
    double compute_total_lifted_content(const Matrix2Xd &vertices) const;

    // compute total lifted content, record lifted content of each triangle
    double compute_total_lifted_content(const Matrix2Xd &vertices, VectorXd &lifted_content_list) const;

    // compute total lifted content and its gradient
    double compute_total_lifted_content_with_gradient(const Matrix2Xd &vertices,
                                        // output
                                        Eigen::Matrix2Xd &grad) const;

    // compute total lifted content and its gradient,
    // and PSD projected Hessian of (TLC - total signed area) on free vertices
    double compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(const Matrix2Xd &vertices,
                                                                                const VectorXi &freeI,
                                                                                const Matrix3Xi  &F_free,
                                                                                // output
                                                                                Matrix2Xd &grad,
                                                                                SpMat &Hess) const;

    // compute total lifted content and its gradient,
    // and PSD projected Hessian of (TLC - total signed area) on free vertices
    double compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(const Matrix2Xd &vertices,
                                                                                const VectorXi &freeI,
                                                                                const Matrix3Xi  &F_free,
                                                                                // output
                                                                                VectorXd &lifted_content_list,
                                                                                Matrix2Xd &grad,
                                                                                SpMat &Hess) const;


private:

    // compute lifted triangle area
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    static double compute_lifted_TriArea(const Matrix2Xd &vert, const Vector3d &r);

    // compute lifted triangle area with gradient wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    static double compute_lifted_TriArea_with_gradient(const Matrix2Xd &vert, const Vector3d &r,
                                                       Matrix2Xd &grad);

    // compute lifted triangle area with gradient and Hessian wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    static double compute_lifted_TriArea_with_gradient_Hessian(const Matrix2Xd &vert, const Vector3d &r,
                                                       Matrix2Xd &grad, MatrixXd &Hess);

private:
    // alpha parameter
//    double param_alpha;
    // dimension of target vertices
//    size_t vDim;
    // indices of triangle vertices
    Eigen::Matrix3Xi F;
    // squared edge lengths of auxiliary triangles
    Eigen::Matrix3Xd restD;

};


#endif //TLC_TOTAL_LIFTED_CONTENT_H
