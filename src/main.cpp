#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double degrees) {
    return degrees * 0.0174532925;
}

double rad_to_deg(double radians) {
    return radians * 57.2957795;
}

Eigen::Matrix3d skew_symmetric(Eigen::Vector3d vec) {
    // Equation (3.30) page 75, MR 3rd print 2019
    return Eigen::Matrix3d{ {0 , -vec[2], vec[1]},
                            {vec[2], 0, -vec[0]},
                            {-vec[1], vec[0], 0} };
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d& x, const Eigen::Vector3d& y,
    const Eigen::Vector3d& z) {
    Eigen::Matrix3d matrix;
    // Equation (3.16) page 65, MR 3rd print 2019
    matrix << x, y, z;
    return matrix;
}

Eigen::Matrix3d rotate_x(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations at page 72, MR 3rd print 2019
    matrix <<
        1, 0, 0,
        0, cos(radians), -sin(radians),
        0, sin(radians), cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_y(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations at page 72, MR 3rd print 2019
    matrix <<
        cos(radians), 0, sin(radians),
        0, 1, 0,
        -sin(radians), 0, cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_z(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations at page 72, MR 3rd print 2019
    matrix <<
        cos(radians), -sin(radians), 0,
        sin(radians), cos(radians), 0,
        0, 0, 1;
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d& axis, double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations at page 72 and Equation (3.52) at page 84, MR 3rd print 2019
    double minus_cos = (1 - cos(radians));
    matrix <<
        cos(radians) + pow(axis[0], 2) * minus_cos,
        axis[0] * axis[1] * minus_cos - axis[2] * sin(radians),
        axis[0] * axis[2] * minus_cos + axis[1] * sin(radians),

        axis[0] * axis[1] * minus_cos + axis[2] * sin(radians),
        cos(radians) + pow(axis[1], 2) * minus_cos,
        axis[1] * axis[2] * minus_cos - axis[0] * sin(radians),

        axis[0] * axis[2] * minus_cos - axis[1] * sin(radians),
        axis[1] * axis[2] * minus_cos + axis[0] * sin(radians),
        cos(radians) + pow(axis[2], 2) * minus_cos;
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d& e) {
    // Appendix B.1 at page 577, MR 3rd print 2019
    Eigen::Matrix3d r_z = rotate_z(e[0]);
    Eigen::Matrix3d r_y = rotate_y(e[1]);
    Eigen::Matrix3d r_x = rotate_x(e[2]);
    return Eigen::Matrix3d::Identity() * r_z * r_y * r_x;
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d& r, const Eigen::Vector3d& p) {
    Eigen::Matrix4d matrix;
    // Equation (3.62) at page 87, MR 3rd print 2019
    matrix << r(0, 0), r(0, 1), r(0, 2), p[0],
        r(1, 0), r(1, 1), r(1, 2), p[1],
        r(2, 0), r(2, 1), r(2, 2), p[2],
        0, 0, 0, 1;
    return matrix;
}

void transform_vector() {
    Eigen::Matrix3d e = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 60.0, 45.0, 0.0 });
    Eigen::Vector3d p{ 0.0, 0.0, 10.0 };
    Eigen::Matrix4d T = transformation_matrix(e, p);
    Eigen::Vector3d v_a{ 2.5, 3.0, -10.0 };

    // Equation (3.65) at page 88, MR 3rd print 2019
    Eigen::Vector4d  x{ v_a[0], v_a[1], v_a[2], 1.0 };
    Eigen::Vector4d Tx = T * x;
    Eigen::Vector3d v_w = { Tx[0], Tx[1], Tx[2] };

    std::cout << "Transform vector: " << std::endl;
    std::cout << v_w.transpose() << std::endl;
}

void skew_symmetric_test() {
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{ 0.5, 0.5, 0.707107 });
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

void rotation_matrix_test() {
    Eigen::Matrix3d rot = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 45.0, -45.0, 90.0 });
    Eigen::Matrix3d rot_aa = rotation_matrix_from_axis_angle(Eigen::Vector3d{ 0.8164966, 0.0,
     0.5773503 }, 120.0);
    Eigen::Matrix3d rot_fa = rotation_matrix_from_frame_axes(Eigen::Vector3d{ 0.5, 0.5, 0.707107 },
        Eigen::Vector3d{ -0.5, -0.5, 0.707107 }, Eigen::Vector3d{ 0.707107, -0.707107, 0.0 });
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}

void transformation_matrix_test() {
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 45, -45.0, 90.0 });
    Eigen::Vector3d v{ 1.0, -2.0, 3.0 };
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}

int main() {
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    return 0;
}
