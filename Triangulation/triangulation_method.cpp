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

std::vector<Vector2D> Normalise2DPoints(const std::vector<Vector2D>& input){
    // Translate the points
    double avg_x = 0, avg_y = 0;
    for (const auto& point : input){
        avg_x += point.x();
        avg_y += point.y();
    }
    avg_x = avg_x / input.size();
    avg_y = avg_y / input.size();

    std::vector<Vector2D> translated_points;
    for (const auto& point : input){
        translated_points.emplace_back(Vector2D(point.x() - avg_x, point.y() - avg_y));
    };

    // Calculate the scaling factor
    double total_distance = 0;
    for (const auto& point : translated_points){
        total_distance += sqrt(pow(point.x(), 2) + pow(point.y(), 2));
    }
    double avg_distance = total_distance / input.size();
    double scale = sqrt(2) / avg_distance;

    // Apply the scaling to the translated points
    for (auto& point : translated_points){
        point = point * scale;
    }
    return translated_points;
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
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

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

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
//    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
//    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
//
//    /// get the number of rows.
//    int num_rows = W.rows();
//
//    /// get the number of columns.
//    int num_cols = W.cols();
//
//    /// get the the element at row 1 and column 2
//    double value = W(1, 2);
//
//    /// get the last column of a matrix
//    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...


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


    std::vector<Vector2D> points_normalized_0 = Normalise2DPoints(points_0);
    std::vector<Vector2D> points_normalized_1 = Normalise2DPoints(points_1);

    Matrix WMatrix = ConstructWMatrix(points_normalized_0, points_normalized_1);

    // single value decomposition
    int m = WMatrix.rows();
    int n = WMatrix.cols();
    Matrix U(m, m, 0.0);   // initialized with 0s
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0);   // initialized with 0s
    svd_decompose(WMatrix, U, S, V);
    Vector fVector = V.get_column(V.cols() - 1);
    Matrix33 FMatrix(fVector[0], fVector[1], fVector[2],
                        fVector[3], fVector[4], fVector[5],
                        fVector[6], fVector[7], fVector[8]);
    std::cout << "F:" << std::endl << FMatrix << std::endl;


    if (points_3d.size() < 8){
        return false;
    }



    return true;
}



