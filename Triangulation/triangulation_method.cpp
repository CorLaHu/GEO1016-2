/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */

void Normalise2DPoints(const std::vector<Vector2D>& input, // input
                       std::vector<Vector2D>& output_points, // output
                       Vector2D& translation, double& scaling // output
                       )
                       {
    // Translate the points
    double total_x = 0, total_y = 0;
    for (const auto& point : input){
        total_x += point.x();
        total_y += point.y();
    }
    translation.x() = total_x / input.size();
    translation.y() = total_y / input.size();

    for (const auto& point : input){
        output_points.emplace_back(Vector2D(point.x() - translation.x(), point.y() - translation.y()));
    };

    // Calculate the scaling factor
    double total_distance = 0;
    for (const auto& point : output_points){
        total_distance += sqrt(pow(point.x(), 2) + pow(point.y(), 2));
    }
    double avg_distance = total_distance / input.size();
    scaling = sqrt(2) / avg_distance;

    // Apply the scaling to the translated points
    for (auto& point : output_points){
        point = point * scaling;
    }
};

Matrix ConstructWMatrix(const std::vector<Vector2D> &points_normalized_0,
                        const std::vector<Vector2D> &points_normalized_1) {
    Matrix wmatrix = Matrix(points_normalized_0.size(), 9, 0.0);

    for (int i = 0; i < points_normalized_0.size(); i++){
        wmatrix(i, 0) = points_normalized_0[i].x() * points_normalized_1[i].x();
        wmatrix(i, 1) = points_normalized_0[i].y() * points_normalized_1[i].x();
        wmatrix(i, 2) = points_normalized_1[i].x();
        wmatrix(i, 3) = points_normalized_0[i].x() * points_normalized_1[i].y();
        wmatrix(i, 4) = points_normalized_0[i].y() * points_normalized_1[i].y();
        wmatrix(i, 5) = points_normalized_1[i].y();
        wmatrix(i, 6) = points_normalized_0[i].x();
        wmatrix(i, 7) = points_normalized_0[i].y();
        wmatrix(i, 8) = 1;
    }
    return wmatrix;
}

std::vector<Vector3D> project_points(std::vector<Vector2D> points_0, std::vector<Vector2D> points_1,
                                     Matrix33 R,Vector3D t, Matrix K, Matrix &p_m0, Matrix &p_m1){
    Matrix34 Rt_0(1.0, 0.0, 0.0, 0.0,
                  0.0, 1.0, 0.0, 0.0,
                  0.0, 0.0, 1.0, 0.0);
    Matrix projection_matrix_0 = K * Rt_0;

    p_m0 = projection_matrix_0;


    Matrix34 Rt_1(R(0, 0), R(0, 1), R(0, 2), t.x(),
                  R(1, 0), R(1, 1), R(1, 2), t.y(),
                  R(2, 0), R(2, 1), R(2, 2), t.z());
    Matrix projection_matrix_1 = K * Rt_1;

    p_m1 = projection_matrix_1;

    // create matrix A
    Vector m1_0 = projection_matrix_0.get_row(0);
    Vector m2_0 = projection_matrix_0.get_row(1);
    Vector m3_0 = projection_matrix_0.get_row(2);
    Vector m1_1 = projection_matrix_1.get_row(0);
    Vector m2_1 = projection_matrix_1.get_row(1);
    Vector m3_1 = projection_matrix_1.get_row(2);

    std::vector<Vector3D> returnpoints;

    for (int i = 0; i < points_0.size(); ++i){
        Matrix A(4, 4, 0.0);
        A.set_row(0, points_0[i].x()* m3_0 - m1_0);
        A.set_row(1, points_0[i].y()* m3_0 - m2_0);
        A.set_row(2, points_1[i].x()* m3_1 - m1_1);
        A.set_row(3, points_1[i].y()* m3_1 - m2_1);
        Matrix U4(4, 4, 0.0);   // initialized with 0s
        Matrix D4(4, 4, 0.0);   // initialized with 0s
        Matrix V4(4, 4, 0.0);   // initialized with 0s
        svd_decompose(A, U4, D4, V4);
        Vector P_inhomogeneous = V4.get_column(V4.cols() - 1);
        Vector P_homogeneous = P_inhomogeneous / P_inhomogeneous[3];
        std::cout << "point " << i << ": " << P_homogeneous << std::endl;
        Vector3D P(P_homogeneous[0], P_homogeneous[1], P_homogeneous[2]);
        returnpoints.push_back(P);
    }
    return returnpoints;
}

bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.

    // TODO: think about whether this is actually correct for an input check.


    std::vector<Vector2D> points_normalized_0;
    Vector2D translation0(0.0, 0.0);
    double scaling0= 0.0;
    Normalise2DPoints(points_0, points_normalized_0, translation0, scaling0);
    std::cout << "scaling 0:" << scaling0 << std::endl;
    std::vector<Vector2D> points_normalized_1;
    Vector2D translation1(0.0, 0.0);
    double scaling1= 0.0;
    Normalise2DPoints(points_1, points_normalized_1, translation1, scaling1);

    Matrix WMatrix = ConstructWMatrix(points_normalized_0, points_normalized_1);

    // Single value decomposition
    int m = WMatrix.rows();
    int n = WMatrix.cols();
    Matrix U(m, m, 0.0);   // initialized with 0s
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0);   // initialized with 0s
    svd_decompose(WMatrix, U, S, V);
    Vector fVector = V.get_column(V.cols() - 1);
    Matrix33 F_estimate_Matrix(fVector[0], fVector[1], fVector[2],
                        fVector[3], fVector[4], fVector[5],
                        fVector[6], fVector[7], fVector[8]);
    std::cout << "F^:" << std::endl << F_estimate_Matrix << std::endl;

    // Constraint enforcement
    Matrix U2(3, 3, 0.0);   // initialized with 0s
    Matrix D(3, 3, 0.0);   // initialized with 0s
    Matrix V2(3, 3, 0.0);   // initialized with 0s
    svd_decompose(F_estimate_Matrix, U2, D, V2);
    Matrix33 D2(D(0, 0), 0.0, 0.0,
                0.0, D(1, 1), 0.0,
                0.0, 0.0, 0.0);
    Matrix FqMatrix = U2 * D2 * V2.transpose();
    std::cout << "Fq:" << std::endl << FqMatrix << std::endl;

    // Denormalization
    Matrix33 T0(scaling0, 0.0, -translation0.x() * scaling0,
                0.0, scaling0, -translation0.y() * scaling0,
                0.0, 0.0, 1.0);
    Matrix33 T1(scaling1, 0.0, -translation1.x() * scaling1,
                0.0, scaling1, -translation1.y() * scaling1,
                0.0, 0.0, 1.0);
    Matrix FMatrix = T1.transpose() * FqMatrix * T0;
    std::cout << "F:" << std::endl << FMatrix << std::endl;


    // Calibration matrices
    // Ask about whether the skew should be 1 or 0 in the K matrix
    Matrix33 K(fx, 0, cx,
               0, fy, cy,
               0, 0, 1);

    Matrix33 EMatrix = K.transpose() * FMatrix * K;

    // Single value decomposition of EMatrix
    Matrix U3(3, 3, 0.0);   // initialized with 0s
    Matrix D3(3, 3, 0.0);   // initialized with 0s
    Matrix V3(3, 3, 0.0);   // initialized with 0s
    svd_decompose(EMatrix, U3, D3, V3);

    Matrix33 WhyMatrix(0.0, -1.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0);

    Matrix33 ZMatrix(0.0, 1.0, 0.0,
                     -1.0, 0.0, 0.0,
                     0.0, 0.0, 0.0);

    struct Rt_options{
        Matrix33 R;
        Vector3D t;
    } option0, option1, option2, option3;

    option0.R = determinant(U3 * WhyMatrix * V3.transpose()) * U3 * WhyMatrix * V3.transpose();
    option0.t = U3.get_column(U3.cols() - 1);

    option1.R = determinant(U3 * WhyMatrix * V3.transpose()) * U3 * WhyMatrix * V3.transpose();
    option1.t = -U3.get_column(U3.cols() - 1);

    option2.R = determinant(U3 * WhyMatrix.transpose() * V3.transpose()) * U3 * WhyMatrix.transpose() * V3.transpose();
    option2.t = U3.get_column(U3.cols() - 1);

    option3.R = determinant(U3 * WhyMatrix.transpose() * V3.transpose()) * U3 * WhyMatrix.transpose() * V3.transpose();
    option3.t = -U3.get_column(U3.cols() - 1);

    std::vector<Rt_options> options = {option0, option1, option2, option3};

    int best_count = 0;
    int index_counter = 0;
    Matrix projection_matrix_0;
    Matrix projection_matrix_1;
    for (const auto& option : options){
        index_counter += 1;
        std::cout << "option: " << index_counter << std::endl;
        std::cout << "R: " << std::endl << option.R << std::endl;
        std::cout << "t: " << std::endl << option.t << std::endl;
        Matrix candidate_projection_matrix_0;
        Matrix candidate_projection_matrix_1;
        // count the amount of correct points, if higher, accept it as current best option
        std::vector<Vector3D> projected_points = project_points(points_0, points_1, option.R, option.t, K, candidate_projection_matrix_0, candidate_projection_matrix_1);
        int count = 0;
        for (const auto& point : projected_points){
            if (point.z() > 0){
                count += 1;
            }
        }
        if (count > best_count){
            best_count = count;
            projection_matrix_0 = candidate_projection_matrix_0;
            projection_matrix_1 = candidate_projection_matrix_1;
            R = option.R;
            t = option.t;
            points_3d = projected_points;
        }
    }


    if (points_3d.size() < 8){
        return false;
    };



    return true;
}



