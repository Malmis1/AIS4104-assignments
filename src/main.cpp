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
    // Equation (3.16) on page 65, MR 3rd print 2019
    matrix << x, y, z;
    return matrix;
}

Eigen::Matrix3d rotate_x(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations on page 72, MR 3rd print 2019
    matrix <<
        1, 0, 0,
        0, cos(radians), -sin(radians),
        0, sin(radians), cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_y(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations on page 72, MR 3rd print 2019
    matrix <<
        cos(radians), 0, sin(radians),
        0, 1, 0,
        -sin(radians), 0, cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_z(double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations on page 72, MR 3rd print 2019
    matrix <<
        cos(radians), -sin(radians), 0,
        sin(radians), cos(radians), 0,
        0, 0, 1;
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d& axis, double degrees) {
    Eigen::Matrix3d matrix;
    double radians = deg_to_rad(degrees);
    // Equations on page 72 and Equation (3.52) on page 84, MR 3rd print 2019
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
    // Appendix B.1 on page 577, MR 3rd print 2019
    Eigen::Matrix3d r_z = rotate_z(e[0]);
    Eigen::Matrix3d r_y = rotate_y(e[1]);
    Eigen::Matrix3d r_x = rotate_x(e[2]);
    return Eigen::Matrix3d::Identity() * r_z * r_y * r_x;
}

Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d& e) {
    // See rotation_matrix_from_euler_zyx function
    Eigen::Matrix3d r_y = rotate_y(e[0]);
    Eigen::Matrix3d r_z = rotate_z(e[1]);
    Eigen::Matrix3d r_x = rotate_x(e[2]);
    return Eigen::Matrix3d::Identity() * r_x * r_y * r_z;
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d& r, const Eigen::Vector3d& p) {
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Zero();
    // Equation (3.62) on page 87, MR 3rd print 2019
    matrix.block<3, 3>(0, 0) = r;
    matrix.block<3, 1>(0, 3) = p;
    matrix(3, 3) = 1;
    return matrix;
}

Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d& r) {
    // Equations on page 579 Section B.1.1, MR 3rd print 2019
    double beta;
    if (r(2, 0) == -1) {
        // Equal to pi/2
        beta = acos(0.0);
    }
    else if (r(2, 0) == 1) {
        // Equal to -pi/2
        beta = -acos(0.0);
    }
    else {
        beta = atan2(-r(2, 0), sqrt(pow(r(0, 0), 2) + pow(r(1, 0), 2)));
    }
    double alpha = atan2(r(1, 0), r(0, 0));
    double gamma = atan2(r(2, 1), r(2, 2));
    return Eigen::Vector3d{ rad_to_deg(alpha),rad_to_deg(beta),rad_to_deg(gamma) };
}

Eigen::VectorXd twist(const Eigen::Vector3d& w, const Eigen::Vector3d& v) {
    // Equation (3.70) on page 96, MR 3rd print 2019
    Eigen::VectorXd vb(6, 1);
    vb << w, v;
    return vb;
}

Eigen::VectorXd screw_axis(const Eigen::Vector3d& q, const Eigen::Vector3d& s, double h) {
    // Equation on page 101, MR 3rd print 2019
    Eigen::VectorXd v(6, 1);
    v << s, -skew_symmetric(s) * q + h * s;
    return v;
}

Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d& tf) {
    // Definition 3.20 on page 98, MR 3rd print 2019
    Eigen::Matrix3d r{ {tf(0,0) , tf(0,1), tf(0,2)},
                            {tf(1,0), tf(1,1), tf(1,2)},
                            {tf(2,0), tf(2,1), tf(2,2)} };
    Eigen::Vector3d p{ tf(0,3), tf(1,3), tf(2,3) };
    Eigen::Matrix3d zeros{ {0,0,0},{0,0,0},{0,0,0} };
    Eigen::MatrixXd adj_matrix(6, 6);
    adj_matrix << r, zeros,
        skew_symmetric(p)* r, r;
    return adj_matrix;
}

double cot(double x) {
    return 1 / (sin(x) / cos(x));
}

Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d& w, double theta) {
    double radians = deg_to_rad(theta);
    // Equation (3.51) on page 82, MR 3rd print 2019
    return Eigen::Matrix3d::Identity() + sin(radians) * skew_symmetric(w)
        + (1 - cos(radians)) * (skew_symmetric(w) * skew_symmetric(w));
}

std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d& r) {
    // Algorithm from page 85 (Equations (3.58)-(3.61) on pages 85-86), MR 3rd print 2019
    Eigen::Vector3d w;
    double theta;
    if (r == Eigen::Matrix3d::Identity()) {
        theta = 0;
    }
    else {
        // Built in function in Eigen for Equation (3.54) on page 84, MR 3rd print 2019
        double trr = r.trace();
        if (trr == -1) {
            // Equal to pi
            theta = acos(0.0) * 2;
            // Equation (3.58) on page 85, MR 3rd print 2019
            w = (1 / sqrt(2 * (1 + r(2, 2))))
                * Eigen::Vector3d{ r(0, 2), r(1, 2), 1 + r(2, 2) };
        }
        else {
            theta = acos((1 / 2) * (trr - 1));
            // Equations directly above Equation (3.53) on page 84, MR 3rd print 2019
            double w_1 = ((1 / (2 * sin(theta))) * (r(2, 1) - r(1, 2)));
            double w_2 = ((1 / (2 * sin(theta))) * (r(0, 2) - r(2, 0)));
            double w_3 = ((1 / (2 * sin(theta))) * (r(1, 0) - r(0, 1)));
            w = Eigen::Vector3d{ w_1, w_2, w_3 };
        }
    }
    return std::pair<Eigen::Vector3d, double>{ w, theta };
}

Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d& w, const Eigen::Vector3d& v, double theta) {
    Eigen::Matrix3d skew_w = skew_symmetric(w);
    Eigen::Matrix4d m_e;
    double radians = deg_to_rad(theta);

    // Proposition 3.25 on page 103, MR 3rd print 2019
    if (w.norm() == 1) {
        Eigen::Matrix3d r = matrix_exponential(w, theta);
        Eigen::Vector3d p = (Eigen::Matrix3d::Identity() * radians + (1 - cos(radians)) * skew_w + (radians - sin(radians)) * skew_w * skew_w) * v;
        m_e = transformation_matrix(r, p);
    }
    else if ((w.norm() == 0) && (v.norm() == 1)) {
        Eigen::Matrix3d r = Eigen::Matrix3d::Identity();
        Eigen::Vector3d p = v * theta;
        m_e = transformation_matrix(r, p);
    }

    return m_e;
}

std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d& t) {
    Eigen::Matrix3d r = t.block<3, 3>(0, 0);
    Eigen::Vector3d p = t.block<3, 1>(0, 3);

    // Algorithm in section 3.3.3.2 on page 104, MR 3rd print 2019
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double h = 0;
    double theta;
    if (r == Eigen::Matrix3d::Identity()) {
        w = Eigen::Vector3d::Zero();
        v = p / p.norm();
        theta = p.norm();
    }
    else {
        std::pair<Eigen::Vector3d, double> m_l = matrix_logarithm(r);
        w = m_l.first;
        theta = m_l.second;
        Eigen::Matrix3d skew_w = skew_symmetric(w);
        v = ((1 / theta) * Eigen::Matrix3d::Identity() - (1 / 2) * skew_w + ((1 / theta) - (1 / 2) * cot(theta / 2)) * skew_w * skew_w) * p;
    }
    return std::pair<Eigen::VectorXd, double>{ screw_axis(w, v, h), theta };
}

Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double>& joint_positions) {
    double l1 = 10.0;
    double l2 = l1;
    double l3 = l2;
    // Equation (4.5) on page 136, MR 3rd print 2019
    Eigen::Matrix4d t_01 = transformation_matrix(rotate_z(joint_positions[0]), Eigen::Vector3d::Zero());
    Eigen::Matrix4d t_12 = transformation_matrix(rotate_z(joint_positions[1]), Eigen::Vector3d{ l1,0,0 });
    Eigen::Matrix4d t_23 = transformation_matrix(rotate_z(joint_positions[2]), Eigen::Vector3d{ l2,0,0 });
    Eigen::Matrix4d t_34 = transformation_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d{ l3,0,0 });
    // Equation (4.4) on page 135, MR 3rd print 2019
    return t_01 * t_12 * t_23 * t_34;
}

Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double>& joint_positions) {
    double l1 = 10.0;
    double l2 = l1;
    double l3 = l2;

    // Example 4.2 from page 142, MR 3rd print 2019
    Eigen::Vector3d w1{ 0,0,1 };
    Eigen::Vector3d v1{ 0, 0, 0 };
    Eigen::Vector3d w2{ 0,0,1 };
    Eigen::Vector3d v2{ 0, -l1, 0 };
    Eigen::Vector3d w3{ 0,0,1 };
    Eigen::Vector3d v3{ 0, -(l1 + l2), 0 };

    Eigen::Matrix4d m = transformation_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d{ l1 + l2 + l3,0,0 });
    Eigen::Matrix4d e_1 = matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Matrix4d e_2 = matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Matrix4d e_3 = matrix_exponential(w3, v3, joint_positions[2]);
    // Equation (4.12) on page 136, MR 3rd print 2019
    return e_1 * e_2 * e_3 * m;
}

Eigen::Matrix4d ur3e_fk_screw(const std::vector<double>& joint_positions) {
    double l1 = -243.5;
    double l2 = -213.2;
    double wi1 = 131.05;
    double wi2 = 92.1;
    double h1 = 151.8;
    double h2 = 85.35;

    // Example 4.5 from page 145, MR 3rd print 2019
    Eigen::Vector3d w1{ 0, 0, 1 };
    Eigen::Vector3d v1{ 0, 0, 0 };
    Eigen::Vector3d w2{ 0, 1, 0 };
    Eigen::Vector3d v2{ -h1, 0, 0 };
    Eigen::Vector3d w3{ 0, 1, 0 };
    Eigen::Vector3d v3{ -h1, 0, l1 };
    Eigen::Vector3d w4{ 0, 1, 0 };
    Eigen::Vector3d v4{ -h1, 0, 0 };
    Eigen::Vector3d w5{ 0, 0, -1 };
    Eigen::Vector3d v5{ -wi1, l1 + l2, 0 };
    Eigen::Vector3d w6{ 0, 1, 0 };
    Eigen::Vector3d v6{ h2 - h1, 0, l1 + l2 };

    Eigen::Matrix4d m{ {-1, 0, 0, l1 + l2},
                   {0, 0, 1, wi1 + wi2},
                   {0, 1, 0, h1 - h2},
                   {0, 0, 0, 1} };
    Eigen::Matrix4d e_1 = matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Matrix4d e_2 = matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Matrix4d e_3 = matrix_exponential(w3, v3, joint_positions[2]);
    Eigen::Matrix4d e_4 = matrix_exponential(w4, v4, joint_positions[3]);
    Eigen::Matrix4d e_5 = matrix_exponential(w5, v5, joint_positions[4]);
    Eigen::Matrix4d e_6 = matrix_exponential(w6, v6, joint_positions[5]);

    return e_1 * e_2 * e_3 * e_4 * e_5 * e_6 * m;
}

Eigen::Matrix4d ur3e_fk_transform(const std::vector<double>& joint_positions) {
    double l1 = -243.5;
    double l2 = -213.2;
    double wi1 = 131.05;
    double wi2 = 92.1;
    double h1 = 151.8;
    double h2 = 85.35;

    // Based on Figure 4.6 (right) on page 146 and Equation (4.5) on page 136, MR 3rd print 2019
    Eigen::Matrix3d r_01 = rotate_z(joint_positions[0]);
    Eigen::Matrix3d r_12 = rotate_y(joint_positions[1]);
    Eigen::Matrix3d r_23 = rotate_y(-joint_positions[2]);
    Eigen::Matrix3d r_34 = rotate_y(joint_positions[3]);
    Eigen::Matrix3d r_45 = rotate_z(-joint_positions[4]);
    Eigen::Matrix3d r_56 = rotate_y(joint_positions[5]);

    // Based on Figure 4.6 (right) on page 146, MR 3rd print 2019
    Eigen::Vector3d p_01{ 0, 0, 0 };
    Eigen::Vector3d p_12{ 0, wi1, h1 };
    Eigen::Vector3d p_23{ l1, 0, 0 };
    Eigen::Vector3d p_34{ l2, 0, 0 };
    Eigen::Vector3d p_45{ 0, 0, -h2 };
    Eigen::Vector3d p_56{ 0, wi2, 0 };
    Eigen::Vector3d p_67{ 0, 0, 0 };

    Eigen::Matrix4d t_01 = transformation_matrix(r_01, p_01);
    Eigen::Matrix4d t_12 = transformation_matrix(r_12, p_12);
    Eigen::Matrix4d t_23 = transformation_matrix(r_23, p_23);
    Eigen::Matrix4d t_34 = transformation_matrix(r_34, p_34);
    Eigen::Matrix4d t_45 = transformation_matrix(r_45, p_45);
    Eigen::Matrix4d t_56 = transformation_matrix(r_56, p_56);
    // Parity with m matrix from PoE implementation
    Eigen::Matrix4d t_67 = transformation_matrix(Eigen::Matrix3d{ {-1,0,0},{0,0,1},{0,1,0} }, p_67);

    // Equation (4.4) on page 135, MR 3rd print 2019
    return t_01 * t_12 * t_23 * t_34 * t_45 * t_56 * t_67;
}

void print_pose(const std::string& label, const Eigen::Matrix4d& tf) {
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);
    Eigen::Vector3d e = euler_zyx_from_rotation_matrix(r);
    std::cout << "Pose \"" << label << "\":" << std::endl;
    std::cout << "Euler ZYX: " << e.transpose() << std::endl;
    std::cout << "Linear position: " << p.transpose() << std::endl << std::endl;
}

void skew_symmetric_test() {
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{ 0.5, 0.5, 0.707107 });
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl << std::endl;
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
    std::cout << "Transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl << std::endl;
}

void transform_vector() {
    Eigen::Matrix3d e = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 60.0, 45.0, 0.0 });
    Eigen::Vector3d p{ 0.0, 0.0, 10.0 };
    Eigen::Matrix4d T = transformation_matrix(e, p);
    Eigen::Vector3d v_a{ 2.5, 3.0, -10.0 };

    // Equation (3.65) on page 88, MR 3rd print 2019
    Eigen::Vector4d  x{ v_a[0], v_a[1], v_a[2], 1.0 };
    Eigen::Vector4d Tx = T * x;
    Eigen::Vector3d v_w = { Tx[0], Tx[1], Tx[2] };

    std::cout << "Transform vector: " << std::endl;
    std::cout << v_w.transpose() << std::endl << std::endl;
}

void euler_zyx_from_rotation_matrix_test() {
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 45, -45.0, 90.0 });
    std::cout << "Euler zyx from rotation matrix: " << std::endl;
    std::cout << euler_zyx_from_rotation_matrix(r) << std::endl << std::endl;
}

void wrench_body_and_sensor_frame() {
    Eigen::Vector3d f_w{ -30,0,0 };
    Eigen::Vector3d m_s{ 0,0,2 };
    Eigen::Vector3d e_ws{ 60,-60,0 };
    Eigen::Matrix3d r = rotation_matrix_from_euler_yzx(e_ws);
    // Equation for moment change of coord frame in summary on page 111, MR 3rd print 2019
    Eigen::Vector3d m_w = r * m_s;
    Eigen::Vector3d f_s = r.transpose() * f_w;
    std::cout << "f_w: " << f_w.transpose() << std::endl;
    std::cout << "m_w: " << m_w.transpose() << std::endl << std::endl;
    std::cout << "f_s: " << f_s.transpose() << std::endl;
    std::cout << "m_s: " << m_s.transpose() << std::endl << std::endl;
}

void sum_of_wrenches_in_different_reference_frames_example() {
    // Example 3.28 on pages 108-109, MR 3rd print 2019
    Eigen::VectorXd f_h(6, 1);
    Eigen::VectorXd f_a(6, 1);
    f_h << 0, 0, 0, 0, -5, 0;
    f_a << 0, 0, 0, 0, 0, 1;
    Eigen::Matrix4d t_hf{
        {1, 0, 0,-0.1},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1} };
    Eigen::Matrix4d t_af{
        {1, 0, 0,-0.25},
        {0, 0, 1, 0},
        {0,-1, 0, 0},
        {0, 0, 0, 1} };

    Eigen::VectorXd f_f(6, 1);
    f_f << (adjoint_matrix(t_hf).transpose() * f_h) + (adjoint_matrix(t_af).transpose() * f_a);

    std::cout << "Calculate the sum of wrenches expressed in different reference frames:" << std::endl;
    std::cout << "f_f: " << f_f.transpose() << std::endl << std::endl;
}

void matrix_exponential_test() {
    // Example 3.12 on page 82, MR 3rd print 2019
    Eigen::Vector3d w{ 0, 0.866, 0.5 };
    double theta = 30;
    std::cout << "Matrix exponential:" << std::endl;
    std::cout << matrix_exponential(w, theta) << std::endl << std::endl;
}

void matrix_logarithm_test() {
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{ 45.0, -45.0, 90.0 });

    std::cout << "Matrix logarithm:" << std::endl;
    std::cout << "w: " << matrix_logarithm(r).first.transpose() << " | theta: " << matrix_logarithm(r).second << std::endl;
    std::cout << "Skew-symmetric of w and theta:" << std::endl;
    std::cout << skew_symmetric(matrix_logarithm(r).first) * rad_to_deg(matrix_logarithm(r).second) << std::endl << std::endl;
}

void planar_3r_fk_transform_test() {
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0},
        {90.0, 0.0, 0.0},
        {0.0, 90.0, 0.0},
        {0.0, 0.0, 90.0},
        {10.0, -15.0, 2.75}
    };
    std::cout << "Planar 3R FK Transform:" << std::endl;
    for (size_t i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d tf = planar_3r_fk_transform(joint_configurations[i]);
        print_pose("Joint " + std::to_string(i + 1), tf);
    }
}

void planar_3r_fk_screw_test() {
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0},
        {90.0, 0.0, 0.0},
        {0.0, 90.0, 0.0},
        {0.0, 0.0, 90.0},
        {10.0, -15.0, 2.75}
    };
    std::cout << "Planar 3R FK Transform using PoE:" << std::endl;
    for (size_t i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d tf = planar_3r_fk_screw(joint_configurations[i]);
        print_pose("Joint " + std::to_string(i + 1), tf);
    }
}

void ur3e_fk_screw_test() {
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}
    };
    std::cout << "Planar 3R FK Transform using PoE:" << std::endl;
    for (size_t i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d tf = ur3e_fk_screw(joint_configurations[i]);
        print_pose("Joint configuration " + std::to_string(i + 1), tf);
    }
}

void ur3e_fk_transform_test() {
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}
    };
    std::cout << "Planar 3R FK Transform using PoE:" << std::endl;
    for (size_t i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d tf = ur3e_fk_transform(joint_configurations[i]);
        print_pose("Joint configuration " + std::to_string(i + 1), tf);
    }
}

int main() {
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    euler_zyx_from_rotation_matrix_test();
    wrench_body_and_sensor_frame();
    sum_of_wrenches_in_different_reference_frames_example();
    matrix_exponential_test();
    matrix_logarithm_test();
    planar_3r_fk_transform_test();
    planar_3r_fk_screw_test();
    ur3e_fk_screw_test();
    ur3e_fk_transform_test();
    return 0;
}
