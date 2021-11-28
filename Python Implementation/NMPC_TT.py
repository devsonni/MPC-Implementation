import casadi as ca
from casadi import sin, cos, pi, tan
import math
import numpy as np
from time import time

# Shift function
def shift_timestep(T, t0, x0, u, f_u, xs):
    f_value = f_u(x0, u[:, 0])
    print(f_value)
    x0 = ca.DM.full(x0 + (T * f_value))

    t0 = t0 + T
    u0 = ca.horzcat(
        u[:, 1:],
        ca.reshape(u[:, -1], -1, 1)
    )

    con_t = [12, 0.001]   # Linear and angular velocity of target
    f_t_value = ca.vertcat(con_t[0] * cos(xs[2]),
                           con_t[0] * sin(xs[2]),
                           con_t[1])
    xs = ca.DM.full(xs + (T * f_t_value))
    return t0, x0, u0, xs

def DM2Arr(dm):
    return np.array(dm.full())

def SX2Arr(sx):
    return np.array(sx.full())

# mpc parameters
# sim_time = 200    # simulation time
T = 0.2             # discrete step
N = 15              # number of look ahead steps

# Constrains of UAV with gimbal
# input constrains of UAV
v_u_min = 14
v_u_max = 30
omega_2_u_min = -pi / 30
omega_2_u_max = pi / 30
omega_3_u_min = -pi / 21
omega_3_u_max = pi / 21

# input constrains of gimbal
omega_1_g_min = -pi / 30
omega_1_g_max = pi / 30
omega_2_g_min = -pi / 30
omega_2_g_max = pi / 30
omega_3_g_min = -pi / 30
omega_3_g_max = pi / 30

# states constrains of UAV
theta_u_min = -0.2618
theta_u_max = 0.2618
z_u_min = 75
z_u_max = 150

# states constrains of gimbal
phi_g_min = -pi / 6
phi_g_max = pi / 6
theta_g_min = -pi / 6
theta_g_max = pi / 6
shi_g_min = -pi / 2
shi_g_max = pi / 2

# Symbolic states of UAV with gimbal camera

# states of the UAV
x_u = ca.SX.sym('x_u')
y_u = ca.SX.sym('y_u')
z_u = ca.SX.sym('z_u')
theta_u = ca.SX.sym('theta_u')
psi_u = ca.SX.sym('psi_u')
# states of the gimbal
phi_g = ca.SX.sym('phi_g')
shi_g = ca.SX.sym('shi_g')
theta_g = ca.SX.sym('theta_g')

# append all UAV states in one vector
states_u = ca.vertcat(
    x_u,
    y_u,
    z_u,
    theta_u,
    psi_u,
    phi_g,
    shi_g,
    theta_g,
)
n_states_u = states_u.numel()

# Controls of UAV that will find by NMPC
# UAV controls parameters
v_u = ca.SX.sym('v_u')
omega_2_u = ca.SX.sym('omega_2_u')
omega_3_u = ca.SX.sym('omega_3_u')
# Gimbal control parameters
omega_1_g = ca.SX.sym('omega_1_g')
omega_2_g = ca.SX.sym('omega_2_g')
omega_3_g = ca.SX.sym('omega_3_g')

# Appending controls in one vector
controls_u = ca.vertcat(
 v_u,
 omega_2_u,
 omega_3_u,
 omega_1_g,
 omega_2_g,
 omega_3_g,
)
n_controls = controls_u.numel()

# Calculates RHS using control vector and current initial states of UAV and gimbal
rhs_u = ca.vertcat(
    v_u*cos(psi_u)*cos(theta_u),
    v_u*sin(psi_u)*cos(theta_u),
    v_u*sin(theta_u),
    omega_2_u,
    omega_3_u,
    omega_1_g,
    omega_2_g,
    omega_3_g
)

# Non-linear mapping function which is f(x,y)
f_u = ca.Function('f', [states_u, controls_u], [rhs_u])

U = ca.SX.sym('U', n_controls, N)       # Decision Variables
P = ca.SX.sym('P', n_states_u + 3)      # This consists of initial states of UAV with gimbal 1-8 and
# reference states 9-11 (reference states is target's states)

X = ca.SX.sym('X', n_states_u, (N+1))    # Has prediction of states over prediction horizon

# Filling the defined system parameters of UAV
X[:, 0] = P[0:8]       # initial state

for k in range(N):
    st = X[:, k]
    con = U[:, k]
    f_value = f_u(st, con)
    st_next = st + T*f_value
    X[:, k+1] = st_next

ff = ca.Function('ff', [U, P], [X])

# Objective function

obj = 0  # objective function that need to minimize with optimal variables
g = []
# g = ca.SX.sym('g', 8*(N+1))  # constrains of pitch angle theta

for k in range(N):
    stt = X[0:8, 0:N]
    A = ca.SX.sym('A', N)
    B = ca.SX.sym('B', N)
    C = ca.SX.sym('C', N)
    X_E = ca.SX.sym('X_E', N)
    Y_E = ca.SX.sym('Y_E', N)

    VFOV = 1    # Making FOV
    HFOV = 1

    w1 = 1     # MPC weights
    w2 = 1

    a = (stt[2, k] * (tan(stt[6, k] + VFOV / 2)) - stt[2, k] * tan(stt[6, k] - VFOV / 2)) / 2   # FOV Stuff
    b = (stt[2, k] * (tan(stt[5, k] + HFOV / 2)) - stt[2, k] * tan(stt[5, k] - HFOV / 2)) / 2

    A[k] = ((cos(stt[7, k])) ** 2) / a ** 2 + ((sin(stt[7, k])) ** 2) / b ** 2
    B[k] = 2 * cos(stt[7, k]) * sin(stt[7, k]) * ((1 / a ** 2) - (1 / b ** 2))
    C[k] = ((sin(stt[7, k])) ** 2) / a ** 2 + ((cos(stt[7, k])) ** 2) / b ** 2

    X_E[k] = stt[0, k] + a + stt[2, k] * (tan(stt[6, k] - VFOV / 2))                            # Centre of FOV
    Y_E[k] = stt[1, k] + b + stt[2, k] * (tan(stt[5, k] - HFOV / 2))

    obj = obj + w1 * math.sqrt((stt[0, k] - P[8]) ** 2 + (stt[1, k] - P[9]) ** 2) + w2 * ((A[k] * (P[8] - X_E[k]) ** 2 + B[k] * (P[9] - Y_E[k]) * (P[8] - X_E[k]) + C[k] * (P[9] - Y_E[k]) ** 2) - 1)

# Obstacle parameters and virtual radius of uav
x_o_1 = 175
y_o_1 = 779
x_o_2 = -134
y_o_2 = 155
x_o_3 = 441
y_o_3 = 343
obs_r = 30
UAV_r = 5

# compute the constrains, states or inequality constrains
for k in range(N+1):
    g = ca.vertcat(g,
    X[2, k],  # limit on hight of UAV
    X[3, k],  # limit on pitch angle theta
    X[5, k],  # limit on gimbal angle phi
    X[6, k],  # limit on gimbal angle theta
    X[7, k],  # limit on gimbal angle shi
    -math.sqrt((X[0, k]-x_o_1)**2 + (X[1, k]-y_o_1)**2) + (UAV_r+obs_r),  # limit of obstacle-1
    -math.sqrt((X[0, k]-x_o_2)**2 + (X[1, k]-y_o_2)**2) + (UAV_r+obs_r),  # limit of obstacle-2
    -math.sqrt((X[0, k]-x_o_3)**2 + (X[1, k]-y_o_3)**2) + (UAV_r+obs_r)   # limit of obstacle-3
)

# make the decision variables one column vector
OPT_variables = \
    U.reshape((-1, 1))    # Example: 6x15 ---> 90x1 where 6=controls, 16=n+1

nlp_prob = {
    'f': obj,
    'x': OPT_variables,
    'g': g,
    'p': P
}

opts = {
    'ipopt': {
        'max_iter': 100,
        'print_level': 0,
        'acceptable_tol': 1e-8,
        'acceptable_obj_change_tol': 1e-6
    },
    'print_time': 0
}

solver = ca.nlpsol('solver', 'ipopt', nlp_prob, opts)

lbx = ca.DM.zeros((n_controls*N, 1))
ubx = ca.DM.zeros((n_controls*N, 1))
lbg = ca.DM.zeros((n_states_u*(N+1)))
ubg = ca.DM.zeros((n_states_u*(N+1)))

# Constrains on states (Inequality constrains)
lbg[0:128:8] = z_u_min         # z lower bound
lbg[1:128:8] = theta_u_min     # theta lower bound
lbg[2:128:8] = phi_g_min       # phi lower bound
lbg[3:128:8] = theta_g_min     # theta lower bound
lbg[4:128:8] = shi_g_min       # shi lower bound
lbg[6:128:8] = -ca.inf         # Obstacle - 1
lbg[5:128:8] = -ca.inf         # Obstacle - 2
lbg[7:128:8] = -ca.inf         # Obstacle - 3

ubg[0:128:8] = z_u_max         # z lower bound
ubg[1:128:8] = theta_u_max     # theta lower bound
ubg[2:128:8] = phi_g_max       # phi lower bound
ubg[3:128:8] = theta_g_max     # theta lower bound
ubg[4:128:8] = shi_g_max       # shi lower bound
ubg[6:128:8] = 0               # Obstacle - 1
ubg[5:128:8] = 0               # Obstacle - 2
ubg[7:128:8] = 0               # Obstacle - 3

# Constrains on controls, constrains on optimization variable
lbx[0: n_controls*N: n_controls, 0] = v_u_min           # velocity lower bound
lbx[1: n_controls*N: n_controls, 0] = omega_2_u_min     # theta 1 lower bound
lbx[2: n_controls*N: n_controls, 0] = omega_3_u_min     # theta 2 lower bound
lbx[3: n_controls*N: n_controls, 0] = omega_1_g_min     # omega 1 lower bound
lbx[4: n_controls*N: n_controls, 0] = omega_2_g_min     # omega 2 lower bound
lbx[5: n_controls*N: n_controls, 0] = omega_3_g_min     # omega 3 lower bound

ubx[0: n_controls*N: n_controls, 0] = v_u_max           # velocity upper bound
ubx[1: n_controls*N: n_controls, 0] = omega_2_u_max     # theta 1 upper bound
ubx[2: n_controls*N: n_controls, 0] = omega_3_u_max     # theta 2 upper bound
ubx[3: n_controls*N: n_controls, 0] = omega_1_g_max     # omega 1 upper bound
ubx[4: n_controls*N: n_controls, 0] = omega_2_g_max     # omega 2 upper bound
ubx[5: n_controls*N: n_controls, 0] = omega_3_g_max     # omega 3 upper bound

args = {
    'lbg': lbg,  # lower bound for state
    'ubg': ubg,  # upper bound for state
    'lbx': lbx,  # lower bound for controls
    'ubx': ubx   # upper bound for controls
}

# Initial state of UAV, gimbal and target
x_init = 90
y_init = 100
z_init = 80
theta_init = 0
psi_init = 0
phi_init = 0
shi_init = 0
theta_2_init = 0

# target
x_target = 100
y_target = 150
theta_target = 0

t0 = 0
x0 = ca.DM(ca.vertcat(x_init,
                      y_init,
                      z_init,
                      theta_init,
                      psi_init,
                      phi_init,
                      shi_init,
                      theta_2_init))  # initial state

xs = ca.DM(ca.vertcat(x_target,
                      y_target,
                      theta_target))  # initial target state

# xx = DM(state_init)
t = ca.DM(t0)

u0 = ca.DM.zeros((n_controls, N))          # initial control
#X0 = ca.repmat(x0, 1, N+1)                # initial state full
xx = ca.DM.zeros((8,20))                   # change size according to main loop run, works as tracker for target and UAV
ss = ca.DM.zeros((3,20))

xx[:,0] = x0
ss[:,0] = xs

mpc_iter = 0
#cat_states = DM2Arr(X0)
cat_controls = DM2Arr(u0[:, 0])
times = np.array([[0]])

###############################################################################
##############################   Main Loop    #################################

if __name__ == '__main__':
    main_loop = time()  # return time in sec
    while mpc_iter < 1:
        t1 = time()
        args['p'] = ca.vertcat(
            x0,    # current state
            xs   # target state
        )
        # optimization variable current state
        args['x0'] = \
            ca.reshape(u0, n_controls*N, 1)

        sol = solver(
            x0=args['x0'],
            lbx=args['lbx'],
            ubx=args['ubx'],
            lbg=args['lbg'],
            ubg=args['ubg'],
            p=args['p']
        )

        u = ca.reshape(sol['x'], n_controls, N)
        ff_value = ff(u, args['p'])

        #cat_states = np.dstack((
        #    cat_states,
        #    ff_value,
        #))

        cat_controls = np.vstack((
            cat_controls,
            DM2Arr(u[:, 0])
        ))

        #print(cat_controls)

        t = np.vstack((
            t,
            t0
        ))

        # U_x = ca.vertcat(1000,10,3,4,0.1,6)   # For debudding
        # print(U_x)
        t0, x0, u0, xs = shift_timestep(T, t0, x0, u0, f_u, xs)

        # tracking states of target and UAV for plotting
        xx[:, mpc_iter+1] = x0
        ss[:, mpc_iter+1] = xs

        print(xx)
        print(ss)  # Helps in debugging

        # xx ...
        t2 = time()
        print(mpc_iter)
        #print(t2-t1)
        times = np.vstack((
            times,
            t2-t1
        ))

        mpc_iter = mpc_iter + 1



