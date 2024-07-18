#include <casadi/casadi.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

using namespace casadi;

// Shift function
std::tuple<double, DM, DM, DM> shift_timestep(double T, double t0, DM x0, DM u, Function f_u, DM xs, int mpc_iter) {
    DM f_value = f_u(DMVector{x0, u(Slice(), 0)})[0];
    x0 = x0 + T * f_value;

    t0 = t0 + T;
    DM u0 = horzcat(u(Slice(), Slice(1, u.size2())), reshape(u(Slice(), -1), -1, 1));

    double v = 12;

    std::vector<double> con_t = {v, 0};
    if (mpc_iter >= 500) {
        con_t = {v, M_PI / 100};
    }
    if (mpc_iter >= 1000) {
        con_t = {v, 0};
    }
    if (mpc_iter >= 1500) {
        con_t = {v, M_PI / 100};
    }

    DM f_t_value = vertcat(con_t[0] * cos(xs(2)),
                           con_t[0] * sin(xs(2)),
                           con_t[1]);

    xs = xs + T * f_t_value;
    return std::make_tuple(t0, x0, u0, xs);
}

// For plotting an ellipse
std::pair<DM, DM> ellipse(double a_p, double b_p, double x_e_1, double y_e_1) {
    DM x = DM::zeros(101);
    DM y = DM::zeros(101);
    std::vector<double> th(101);
    for (int i = 0; i < 101; ++i) {
        th[i] = 2 * M_PI * i / 100;
    }
    for (int i = 0; i < 101; ++i) {
        x(i) = a_p * sin(th[i]) + x_e_1;
        y(i) = b_p * cos(th[i]) + y_e_1;
    }
    return {x, y};
}

// Main function
int main() {
    double T = 0.2;
    int N = 15;

    // Constrains of UAV with gimbal
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

    // States constrains of UAV
    double theta_u_min = -0.2618;
    double theta_u_max = 0.2618;
    double z_u_min = 75;
    double z_u_max = 150;

    // States constrains of gimbal
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

    SX states_u = vertcat(x_u, y_u, z_u, theta_u, psi_u, phi_g, shi_g, theta_g);
    int n_states_u = states_u.size1();

    // Controls of UAV that will find by NMPC
    SX v_u = SX::sym("v_u");
    SX omega_2_u = SX::sym("omega_2_u");
    SX omega_3_u = SX::sym("omega_3_u");
    SX omega_1_g = SX::sym("omega_1_g");
    SX omega_2_g = SX::sym("omega_2_g");
    SX omega_3_g = SX::sym("omega_3_g");

    SX controls_u = vertcat(v_u, omega_2_u, omega_3_u, omega_1_g, omega_2_g, omega_3_g);
    int n_controls = controls_u.size1();

    SX rhs_u = vertcat(v_u * cos(psi_u) * cos(theta_u),
                       v_u * sin(psi_u) * cos(theta_u),
                       v_u * sin(theta_u),
                       omega_2_u,
                       omega_3_u,
                       omega_1_g,
                       omega_2_g,
                       omega_3_g);

    Function f_u = Function("f", {states_u, controls_u}, {rhs_u});

    SX U = SX::sym("U", n_controls, N);
    SX P = SX::sym("P", n_states_u + 3);

    SX X = SX::sym("X", n_states_u, N + 1);
    X(Slice(), 0) = P(Slice(0, 8));

    for (int k = 0; k < N; ++k) {
        SX st = X(Slice(), k);
        SX con = U(Slice(), k);
        SX f_value = f_u(SXVector{st, con})[0];
        SX st_next = st + T * f_value;
        X(Slice(), k + 1) = st_next;
    }

    Function ff = Function("ff", {U, P}, {X});

    SX obj = 0;
    SX g;

    for (int k = 0; k < N; ++k) {
        SX stt = X(Slice(0, 8), Slice(0, N));
        SX A = SX::sym("A", N);
        SX B = SX::sym("B", N);
        SX C = SX::sym("C", N);
        SX X_E = SX::sym("X_E", N);
        SX Y_E = SX::sym("Y_E", N);

        double VFOV = 1;
        double HFOV = 1;

        double w1 = 1;
        double w2 = 2;

        SX a = (stt(2, k) * tan(stt(6, k) + VFOV / 2) - stt(2, k) * tan(stt(6, k) - VFOV / 2)) / 2;
        SX b = (stt(2, k) * tan(stt(5, k) + HFOV / 2) - stt(2, k) * tan(stt(5, k) - HFOV / 2)) / 2;

        A(k) = pow(cos(stt(7, k)), 2) / pow(a, 2) + pow(sin(stt(7, k)), 2) / pow(b, 2);
        B(k) = 2 * cos(stt(7, k)) * sin(stt(7, k)) * (1 / pow(a, 2) - 1 / pow(b, 2));
        C(k) = pow(sin(stt(7, k)), 2) / pow(a, 2) + pow(cos(stt(7, k)), 2) / pow(b, 2);

        X_E(k) = stt(0, k) + a + stt(2, k) * tan(stt(6, k) - VFOV / 2);
        Y_E(k) = stt(1, k) + b + stt(2, k) * tan(stt(5, k) - HFOV / 2);

        obj = obj + w1 * sqrt(pow(stt(0, k) - P(8), 2) + pow(stt(1, k) - P(9), 2)) +
              w2 * ((A(k) * pow(P(8) - X_E(k), 2) + B(k) * (P(9) - Y_E(k)) * (P(8) - X_E(k)) + C(k) * pow(P(9) - Y_E(k), 2)) - 1);
    }

    std::vector<double> obs_x = {0, 500, 1000};
    std::vector<double> obs_y = {300, 800, 300};
    double r = 100;

    for (int k = 0; k < N + 1; ++k) {
        for (int j = 0; j < obs_x.size(); ++j) {
            SX X_obs = SX::sym("X_obs", N);
            SX Y_obs = SX::sym("Y_obs", N);
            SX d_obs = SX::sym("d_obs", N);

            X_obs(j) = obs_x[j];
            Y_obs(j) = obs_y[j];
            d_obs(j) = sqrt(pow(X(0, k) - X_obs(j), 2) + pow(X(1, k) - Y_obs(j), 2));

            obj = obj + 0.5 * (1 / (d_obs(j) - r));
        }
    }

    for (int k = 0; k < N + 1; ++k) {
        g = vertcat(g, X(2, k));
        g = vertcat(g, X(3, k));
        g = vertcat(g, X(5, k));
        g = vertcat(g, X(6, k));
        g = vertcat(g, X(7, k));
    }

    for (int k = 0; k < N; ++k) {
        g = vertcat(g, U(0, k));
        g = vertcat(g, U(1, k));
        g = vertcat(g, U(2, k));
        g = vertcat(g, U(3, k));
        g = vertcat(g, U(4, k));
        g = vertcat(g, U(5, k));
    }

    SX OPT_variables = vertcat(reshape(U, -1, 1));
    SX nlp_prob = SX::sym("nlp_prob");
    nlp_prob = vertcat(OPT_variables, P);
    SX lbx = SX::sym("lbx", n_controls * N, 1);
    SX ubx = SX::sym("ubx", n_controls * N, 1);

    for (int i = 0; i < N; ++i) {
        lbx(i * 6 + 0) = v_u_min;
        lbx(i * 6 + 1) = omega_2_u_min;
        lbx(i * 6 + 2) = omega_3_u_min;
        lbx(i * 6 + 3) = omega_1_g_min;
        lbx(i * 6 + 4) = omega_2_g_min;
        lbx(i * 6 + 5) = omega_3_g_min;
        ubx(i * 6 + 0) = v_u_max;
        ubx(i * 6 + 1) = omega_2_u_max;
        ubx(i * 6 + 2) = omega_3_u_max;
        ubx(i * 6 + 3) = omega_1_g_max;
        ubx(i * 6 + 4) = omega_2_g_max;
        ubx(i * 6 + 5) = omega_3_g_max;
    }

    SX lbg = SX::zeros(g.size1(), 1);
    SX ubg = SX::zeros(g.size1(), 1);
    for (int k = 0; k < N + 1; ++k) {
        lbg(k * 5 + 0) = z_u_min;
        lbg(k * 5 + 1) = theta_u_min;
        lbg(k * 5 + 2) = phi_g_min;
        lbg(k * 5 + 3) = shi_g_min;
        lbg(k * 5 + 4) = theta_g_min;
        ubg(k * 5 + 0) = z_u_max;
        ubg(k * 5 + 1) = theta_u_max;
        ubg(k * 5 + 2) = phi_g_max;
        ubg(k * 5 + 3) = shi_g_max;
        ubg(k * 5 + 4) = theta_g_max;
    }

    Dict nlp_solver_opts;
    nlp_solver_opts["ipopt.max_iter"] = 5000;
    nlp_solver_opts["ipopt.tol"] = 1e-8;
    nlp_solver_opts["ipopt.acceptable_tol"] = 1e-6;
    nlp_solver_opts["ipopt.acceptable_obj_change_tol"] = 1e-4;

    Function nlp_solver = nlpsol("nlp_solver", "ipopt", nlp_prob, nlp_solver_opts);
    DMDict arg;
    arg["x0"] = DM::zeros(n_controls * N, 1);
    arg["lbx"] = lbx;
    arg["ubx"] = ubx;
    arg["lbg"] = lbg;
    arg["ubg"] = ubg;
    arg["p"] = DM::zeros(n_states_u + 3, 1);

    double t0 = 0;
    DM x0 = DM::zeros(n_states_u, 1);
    DM xs = DM::zeros(3, 1);
    xs(0) = 500;
    xs(1) = 1500;

    double v = 0;
    double omega = 0;
    DM u0 = DM::zeros(n_controls, N);

    int mpc_iter = 0;

    while (mpc_iter < 500) {
        auto [t0_new, x0_new, u0_new, xs_new] = shift_timestep(T, t0, x0, u0, f_u, xs, mpc_iter);
        t0 = t0_new;
        x0 = x0_new;
        u0 = u0_new;
        xs = xs_new;

        arg["p"] = vertcat(x0, xs);
        arg["x0"] = reshape(u0, n_controls * N, 1);
        DMDict res = nlp_solver(arg);
        DM u = reshape(res["x"], n_controls, N);

        std::cout << "Iteration: " << mpc_iter << std::endl;
        std::cout << "x0: " << x0 << std::endl;
        std::cout << "u: " << u << std::endl;

        mpc_iter++;
    }

    return 0;
}
