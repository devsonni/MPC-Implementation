#include <iostream>
#include <casadi/casadi.hpp>
#include <vector>
#include <cmath>
#include <ctime>

using namespace casadi;
using namespace std;

// Shift function
tuple<double, DM, DM, DM> shift_timestep(double T, double t0, DM x0, DM u, Function f_u, DM xs, int mpc_iter) {
    DM f_value = f_u(DMVector{x0, u(Slice(), Slice(0, 1))})[0];
    x0 = DM::full(x0 + T * f_value);
    t0 = t0 + T;
    u = horzcat(u(Slice(), Slice(1, u.size2())), reshape(u(Slice(), Slice(u.size2() - 1, u.size2())), -1, 1));

    double v = 14;
    vector<double> con_t = {v, 0};

    if (mpc_iter >= 300) con_t = {v, -(M_PI / 2) / 24};
    if (mpc_iter >= 360) con_t = {v, 0};
    if (mpc_iter >= 410) con_t = {v, (M_PI / 2) / 24};
    if (mpc_iter >= 470) con_t = {v, 0};
    if (mpc_iter >= 570) con_t = {v, ((11 * M_PI) / 18) / 12};
    if (mpc_iter >= 630) con_t = {v, 0};
    if (mpc_iter >= 780) con_t = {v, ((7 * M_PI) / 18) / 12};
    if (mpc_iter >= 840) con_t = {v, 0};
    if (mpc_iter >= 940) con_t = {v, -(3 * M_PI / 18) / 12};
    if (mpc_iter >= 1000) con_t = {v, 0};
    if (mpc_iter >= 1100) con_t = {v, (3 * M_PI / 18) / 12};
    if (mpc_iter >= 1160) con_t = {v, 0};
    if (mpc_iter >= 1335) con_t = {v, (M_PI / 2) / 12};
    if (mpc_iter >= 1395) con_t = {v, 0};
    if (mpc_iter >= 1535) con_t = {v, (M_PI / 2) / 12};

    DM f_t_value = vertcat(con_t[0] * cos(xs(2)), con_t[0] * sin(xs(2)), con_t[1]);
    xs = DM::full(xs + T * f_t_value);

    return make_tuple(t0, x0, u, xs);
}

// Helper functions
vector<double> DM2Arr(DM dm) {
    vector<double> arr(dm.size1() * dm.size2());
    for (int i = 0; i < dm.size1(); ++i) {
        for (int j = 0; j < dm.size2(); ++j) {
            arr[i * dm.size2() + j] = dm(i, j).scalar();
        }
    }
    return arr;
}

// For plotting a cylinder
tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> data_for_cylinder_along_z(double center_x, double center_y, double radius, double height_z) {
    vector<double> z(50), theta(50);
    for (int i = 0; i < 50; ++i) {
        z[i] = height_z * i / 49.0;
        theta[i] = 2 * M_PI * i / 49.0;
    }

    vector<vector<double>> theta_grid(50, vector<double>(50)), z_grid(50, vector<double>(50));
    for (int i = 0; i < 50; ++i) {
        for (int j = 0; j < 50; ++j) {
            theta_grid[i][j] = theta[j];
            z_grid[i][j] = z[i];
        }
    }

    vector<vector<double>> x_grid(50, vector<double>(50)), y_grid(50, vector<double>(50));
    for (int i = 0; i < 50; ++i) {
        for (int j = 0; j < 50; ++j) {
            x_grid[i][j] = radius * cos(theta_grid[i][j]) + center_x;
            y_grid[i][j] = radius * sin(theta_grid[i][j]) + center_y;
        }
    }

    return make_tuple(x_grid, y_grid, z_grid);
}

// For plotting ellipse
tuple<DM, DM> ellipse(double a_p, double b_p, double x_e_1, double y_e_1) {
    DM x = DM::zeros(101);
    DM y = DM::zeros(101);
    for (int i = 0; i < 100; ++i) {
        double th = 2 * M_PI * i / 99.0;
        x(i) = a_p * sin(th) + x_e_1;
        y(i) = b_p * cos(th) + y_e_1;
    }
    return make_tuple(x, y);
}

int main() {
    // mpc parameters
    double T = 0.2; // discrete step
    int N = 15;     // number of look ahead steps

    // Constraints of UAV with gimbal
    double v_u_min = 14;
    double v_u_max = 30;
    double omega_2_u_min = -M_PI / 30;
    double omega_2_u_max = M_PI / 30;
    double omega_3_u_min = -M_PI / 21;
    double omega_3_u_max = M_PI / 21;

    double omega_1_g_min = -M_PI / 30;
    double omega_1_g_max = M_PI / 30;
    double omega_2_g_min = -M_PI / 30;
    double omega_2_g_max = M_PI / 30;
    double omega_3_g_min = -M_PI / 30;
    double omega_3_g_max = M_PI / 30;

    double theta_u_min = -0.2618;
    double theta_u_max = 0.2618;
    double z_u_min = 75;
    double z_u_max = 150;

    double phi_g_min = -M_PI / 6;
    double phi_g_max = M_PI / 6;
    double theta_g_min = -M_PI / 6;
    double theta_g_max = M_PI / 6;
    double shi_g_min = -M_PI / 2;
    double shi_g_max = M_PI / 2;

    // Symbolic states of UAV with gimbal camera
    SX x_u = SX::sym("x_u");
    SX y_u = SX::sym("y_u");
    SX z_u = SX::sym("z_u");
    SX theta_u = SX::sym("theta_u");
    SX psi_u = SX::sym("psi_u");
    SX phi_g = SX::sym("phi_g");
    SX shi_g = SX::sym("shi_g");
    SX theta_g = SX::sym("theta_g");

    SX states_u = vertcat({x_u, y_u, z_u, theta_u, psi_u, phi_g, shi_g, theta_g});
    int n_states_u = states_u.size1();

    SX v_u = SX::sym("v_u");
    SX omega_2_u = SX::sym("omega_2_u");
    SX omega_3_u = SX::sym("omega_3_u");
    SX omega_1_g = SX::sym("omega_1_g");
    SX omega_2_g = SX::sym("omega_2_g");
    SX omega_3_g = SX::sym("omega_3_g");

    SX controls_u = vertcat({v_u, omega_2_u, omega_3_u, omega_1_g, omega_2_g, omega_3_g});
    int n_controls = controls_u.size1();

    SX rhs_u = vertcat({v_u * cos(psi_u) * cos(theta_u),
                        v_u * sin(psi_u) * cos(theta_u),
                        v_u * sin(theta_u),
                        omega_2_u,
                        omega_3_u,
                        omega_1_g,
                        omega_2_g,
                        omega_3_g});

    Function f_u = Function("f_u", {states_u, controls_u}, {rhs_u});

    SX U = SX::sym("U", n_controls, N);
    SX P = SX::sym("P", n_states_u + 3);
    SX X = SX::sym("X", n_states_u, N + 1);

    X(Slice(), Slice(0, 1)) = P(Slice(0, 8));

    for (int k = 0; k < N; ++k) {
        SX st = X(Slice(), k);
        SX con = U(Slice(), k);
        SX f_value = f_u(st, con);
        SX st_next = st + T * f_value;
        X(Slice(), k + 1) = st_next;
    }

    Function ff = Function("ff", {U, P}, {X});

    // Objective and constraints
    SX obj = 0;
    SX g = SX::sym("g", n_states_u * (N + 1) + 3 * N, 1);

    vector<double> Q = {10, 10, 1, 1, 1, 1, 1, 1};
    vector<double> R = {0.5, 0.05, 0.05, 0.05, 0.05, 0.05};

    for (int k = 0; k < N; ++k) {
        SX st = X(Slice(), k);
        SX con = U(Slice(), k);
        obj = obj + mtimes(mtimes((st - P(Slice(8, 16))).T(), SX::diag(SX(Q))), (st - P(Slice(8, 16)))) +
              mtimes(mtimes(con.T(), SX::diag(SX(R))), con);
    }

    SX OPT_variables = reshape(U, 6 * N, 1);

    g(Slice(0, n_states_u * (N + 1))) = reshape(X, n_states_u * (N + 1), 1);
    int g_index = n_states_u * (N + 1);

    for (int k = 0; k < N; ++k) {
        SX stt = X(Slice(0, 2), k + 1);
        SX p_c = P(Slice(16, 19));
        SX f_p = stt - p_c;
        g(g_index) = SX::norm_2(f_p);
        g_index += 1;
    }

    SXDict nlp_prob = {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}};

    Dict opts;
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = 0;
    opts["ipopt.max_iter"] = 2000;
    opts["ipopt.acceptable_tol"] = 1e-8;
    opts["ipopt.acceptable_obj_change_tol"] = 1e-6;

    Function solver = nlpsol("solver", "ipopt", nlp_prob, opts);

    vector<double> lbx_vec(6 * N, 0), ubx_vec(6 * N, 0);
    vector<double> lbg_vec(n_states_u * (N + 1) + N, 0), ubg_vec(n_states_u * (N + 1) + N, 0);

    for (int k = 0; k < N; ++k) {
        lbx_vec[k * 6 + 0] = v_u_min;
        ubx_vec[k * 6 + 0] = v_u_max;

        lbx_vec[k * 6 + 1] = omega_2_u_min;
        ubx_vec[k * 6 + 1] = omega_2_u_max;

        lbx_vec[k * 6 + 2] = omega_3_u_min;
        ubx_vec[k * 6 + 2] = omega_3_u_max;

        lbx_vec[k * 6 + 3] = omega_1_g_min;
        ubx_vec[k * 6 + 3] = omega_1_g_max;

        lbx_vec[k * 6 + 4] = omega_2_g_min;
        ubx_vec[k * 6 + 4] = omega_2_g_max;

        lbx_vec[k * 6 + 5] = omega_3_g_min;
        ubx_vec[k * 6 + 5] = omega_3_g_max;
    }

    for (int k = 0; k < n_states_u * (N + 1); ++k) {
        lbg_vec[k] = -numeric_limits<double>::infinity();
        ubg_vec[k] = numeric_limits<double>::infinity();
    }

    for (int k = n_states_u * (N + 1); k < n_states_u * (N + 1) + N; ++k) {
        lbg_vec[k] = 0;
        ubg_vec[k] = 30;
    }

    // Initial and desired states
    vector<double> t0(1, 0);
    vector<double> x0 = {0, 0, 100, 0, 0, 0, 0, 0};
    vector<double> xs = {50, 50, 90, 0, 0, 0, 0, 0};
    vector<double> u0(6 * N, 0);
    vector<double> xx = x0;
    vector<double> t(1, 0);

    double mpciter = 0;
    vector<double> xs_p = {0, 0, 100};
    double a_p = 40;
    double b_p = 20;

    while (mpciter < 200) {
        vector<double> x_e = DM2Arr(DM({xs_p[0] + a_p * sin(M_PI * t[0] / 12),
                                        xs_p[1] + b_p * cos(M_PI * t[0] / 12),
                                        90}));
        vector<double> p = x0;
        p.insert(p.end(), xs.begin(), xs.end());
        p.insert(p.end(), x_e.begin(), x_e.end());

        DMDict arg;
        arg["lbx"] = lbx_vec;
        arg["ubx"] = ubx_vec;
        arg["lbg"] = lbg_vec;
        arg["ubg"] = ubg_vec;
        arg["p"] = p;
        arg["x0"] = u0;

        DMDict res = solver(arg);
        DM u = res["x"];
        vector<double> u_m = DM2Arr(reshape(u, 6, N));

        for (int i = 0; i < u_m.size(); ++i) {
            u0[i] = u_m[i];
        }

        vector<double> shift_res = DM2Arr(std::get<1>(shift_timestep(T, t[0], DM(x0), DM(u0), f_u, DM(xs), mpciter)));

        for (int i = 0; i < shift_res.size(); ++i) {
            x0[i] = shift_res[i];
        }

        xx.push_back(x0[0]);
        t[0] = t[0] + T;
        mpciter += 1;

        cout << "MPC iteration: " << mpciter << "   X: " << x0[0] << " Y: " << x0[1] << " Z: " << x0[2] << endl;
    }

    return 0;
}

