import casadi as ca
from casadi import sin, cos, pi, tan
import math
import numpy as np
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.cbook import flatten
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
from tvtk.api import tvtk
from mayavi.modules.api import Outline, GridPlane


# Shift function
def shift_timestep(T, t0, x0, u, f_u, xs):
    #print(u[:, 0])
    f_value = f_u(x0, u[:, 0])
    #print(f_value)
    x0 = ca.DM.full(x0 + (T * f_value))

    t0 = t0 + T
    u0 = ca.horzcat(
        u[:, 1:],
        ca.reshape(u[:, -1], -1, 1)
    )

    v = 12

    con_t = [v, 0]   # Linear and angular velocity of target
    if(mpc_iter>=500):
        con_t = [v, (pi)/100] #right turn
    if (mpc_iter >= 1000):
        con_t = [v, 0]  # 50 ahead
    if (mpc_iter >= 1500):
        con_t = [v, (pi)/100]  #left ahead

    f_t_value = ca.vertcat(con_t[0] * cos(xs[2]),
                           con_t[0] * sin(xs[2]),
                           con_t[1])
    #print(xs)
    xs = ca.DM.full(xs + (T * f_t_value))
    #print(xs)
    return t0, x0, u0, xs

def DM2Arr(dm):
    return np.array(dm.full())

def SX2Arr(sx):
    return np.array(sx.full())

# For plotting a cylinder
def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

# For plotting ellipse
def ellipse(a_p, b_p, x_e_1, y_e_1):
    x = ca.DM.zeros((101))
    y = ca.DM.zeros((101))
    th = np.linspace(0,  2*np.pi, 100)
    x = a_p * sin(th) + x_e_1
    y = b_p * cos(th) + y_e_1
    return x, y

# mpc parameters
T = 0.2 # discrete step
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

count = 0
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
    w2 = 2
    #wi[count]
    #wii[count]

    a = (stt[2, k] * (tan(stt[6, k] + VFOV / 2)) - stt[2, k] * tan(stt[6, k] - VFOV / 2)) / 2   # FOV Stuff
    b = (stt[2, k] * (tan(stt[5, k] + HFOV / 2)) - stt[2, k] * tan(stt[5, k] - HFOV / 2)) / 2

    A[k] = ((cos(stt[7, k])) ** 2) / a ** 2 + ((sin(stt[7, k])) ** 2) / b ** 2
    B[k] = 2 * cos(stt[7, k]) * sin(stt[7, k]) * ((1 / a ** 2) - (1 / b ** 2))
    C[k] = ((sin(stt[7, k])) ** 2) / a ** 2 + ((cos(stt[7, k])) ** 2) / b ** 2

    X_E[k] = stt[0, k] + a + stt[2, k] * (tan(stt[6, k] - VFOV / 2))                            # Centre of FOV
    Y_E[k] = stt[1, k] + b + stt[2, k] * (tan(stt[5, k] - HFOV / 2))

    obj = obj + w1 * ca.sqrt((stt[0, k] - P[8]) ** 2 + (stt[1, k] - P[9]) ** 2) + \
          w2 * ((A[k] * (P[8] - X_E[k]) ** 2 + B[k] * (P[9] - Y_E[k]) * (P[8] - X_E[k]) + C[k] * (P[9] - Y_E[k]) ** 2) - 1)
    count += 1

# Obstacle parameters and virtual radius of uav
x_o_1 = 0
y_o_1 = 80
x_o_2 = 500
y_o_2 = 245
x_o_3 = 1000
y_o_3 = 70
x_o_4 = 1500
y_o_4 = 295
x_o_5 = 1765
y_o_5 = 550
x_o_6 = 1500
y_o_6 = 750
x_o_7 = 1000
y_o_7 = 1005
x_o_8 = 500
y_o_8 = 800
x_o_9 = -100
y_o_9 = 950
x_o_10 = -200
y_o_10 = 550
obs_r = 50
UAV_r = 5

# compute the constrains, states or inequality constrains
for k in range(N+1):
    g = ca.vertcat(g,
    X[2, k],  # limit on hight of UAV
    X[3, k],  # limit on pitch angle theta
    X[5, k],  # limit on gimbal angle phi
    X[6, k],  # limit on gimbal angle theta
    X[7, k],  # limit on gimbal angle shi
    -ca.sqrt((X[0, k]-x_o_1)**2 + (X[1, k]-y_o_1)**2) + (UAV_r+obs_r),  # limit of obstacle-1
    -ca.sqrt((X[0, k]-x_o_2)**2 + (X[1, k]-y_o_2)**2) + (UAV_r+obs_r),  # limit of obstacle-2
    -ca.sqrt((X[0, k]-x_o_3)**2 + (X[1, k]-y_o_3)**2) + (UAV_r+obs_r),  # limit of obstacle-3
    -ca.sqrt((X[0, k] - x_o_4) ** 2 + (X[1, k] - y_o_4) ** 2) + (UAV_r + obs_r),    # limit of obstacle-4
    -ca.sqrt((X[0, k] - x_o_5) ** 2 + (X[1, k] - y_o_5) ** 2) + (UAV_r + obs_r),    # limit of obstacle-5
    -ca.sqrt((X[0, k] - x_o_6) ** 2 + (X[1, k] - y_o_6) ** 2) + (UAV_r + obs_r),   # limit of obstacle-6
    -ca.sqrt((X[0, k] - x_o_7) ** 2 + (X[1, k] - y_o_7) ** 2) + (UAV_r + obs_r),    # limit of obstacle-7
    -ca.sqrt((X[0, k] - x_o_8) ** 2 + (X[1, k] - y_o_8) ** 2) + (UAV_r + obs_r),    # limit of obstacle-8
    -ca.sqrt((X[0, k] - x_o_9) ** 2 + (X[1, k] - y_o_9) ** 2) + (UAV_r + obs_r),    # limit of obstacle-9
    -ca.sqrt((X[0, k] - x_o_10) ** 2 + (X[1, k] - y_o_10) ** 2) + (UAV_r + obs_r)   # limit of obstacle-10
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
lbg = ca.DM.zeros(((n_states_u+7)*(N+1)))
ubg = ca.DM.zeros(((n_states_u+7)*(N+1)))

# Constrains on states (Inequality constrains)
lbg[0:240:15] = z_u_min         # z lower bound
lbg[1:240:15] = theta_u_min     # theta lower bound
lbg[2:240:15] = phi_g_min       # phi lower bound
lbg[3:240:15] = theta_g_min     # theta lower bound
lbg[4:240:15] = shi_g_min       # shi lower bound
lbg[6:240:15] = -ca.inf         # Obstacle - 1
lbg[5:240:15] = -ca.inf         # Obstacle - 2
lbg[7:240:15] = -ca.inf         # Obstacle - 3
lbg[8:240:15] = -ca.inf         # Obstacle - 4
lbg[9:240:15] = -ca.inf         # Obstacle - 5
lbg[10:240:15] = -ca.inf         # Obstacle - 6
lbg[11:240:15] = -ca.inf         # Obstacle - 7
lbg[12:240:15] = -ca.inf         # Obstacle - 8
lbg[13:240:15] = -ca.inf         # Obstacle - 9
lbg[14:240:15] = -ca.inf         # Obstacle - 10


ubg[0:240:15] = z_u_max         # z lower bound
ubg[1:240:15] = theta_u_max     # theta lower bound
ubg[2:240:15] = phi_g_max       # phi lower bound
ubg[3:240:15] = theta_g_max     # theta lower bound
ubg[4:240:15] = shi_g_max       # shi lower bound
ubg[5:240:15] = 0               # Obstacle - 1
ubg[6:240:15] = 0               # Obstacle - 2
ubg[7:240:15] = 0               # Obstacle - 3
ubg[8:240:15] = 0               # Obstacle - 4
ubg[9:240:15] = 0               # Obstacle - 5
ubg[10:240:15] = 0              # Obstacle - 6
ubg[11:240:15] = 0              # Obstacle - 7
ubg[12:240:15] = 0              # Obstacle - 8
ubg[13:240:15] = 0              # Obstacle - 9
ubg[14:240:15] = 0              # Obstacle - 10

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

# target
x_target = 100
y_target = 150
theta_target = 0

t0 = 0
x0 = [99, 150, 80, 0, 0, 0, 0, 0]
xs = ca.DM(ca.vertcat(x_target,
                      y_target,
                      theta_target))  # initial target state

# xx = DM(state_init)
t = ca.DM(t0)
loop_run = 2000

u0 = ca.DM.zeros((n_controls, N))          # initial control
xx = ca.DM.zeros((8,loop_run+1))                  # change size according to main loop run, works as tracker for target and UAV
ss = ca.DM.zeros((3,loop_run+1))
x_e_1 = ca.DM.zeros((loop_run+1))
y_e_1 = ca.DM.zeros((loop_run+1))

xx[:,0] = x0
ss[:,0] = xs

mpc_iter = 0
#loop_run = 822
cat_controls = DM2Arr(u0[:, 0])
times = np.array([[0]])

###############################################################################
##############################   Main Loop    #################################

if __name__ == '__main__':
    main_loop = time()  # return time in sec
    while mpc_iter < loop_run:
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

        cat_controls = np.vstack((
            cat_controls,
            DM2Arr(u[:, 0])
        ))

        #print(cat_controls)

        t = np.vstack((
            t,
            t0
        ))

        t0, x0, u0, xs = shift_timestep(T, t0, x0, u, f_u, xs)

        # tracking states of target and UAV for plotting
        xx[:, mpc_iter+1] = x0
        ss[:, mpc_iter+1] = xs

        # xx ...
        t2 = time()
        print(mpc_iter)
        #print(t2-t1)
        times = np.vstack((
            times,
            t2-t1
        ))

        mpc_iter = mpc_iter + 1

        a_p = (x0[2] * (tan(x0[6] + VFOV / 2)) - x0[2] * tan(x0[6] - VFOV / 2)) / 2 # For plotting FOV
        b_p = (x0[2] * (tan(x0[5] + HFOV / 2)) - x0[2] * tan(x0[5] - HFOV / 2)) / 2
        x_e_1[mpc_iter] = x0[0] + a_p + x0[2] * (tan(x0[6] - VFOV / 2))
        y_e_1[mpc_iter] = x0[1] + b_p + x0[2] * (tan(x0[5] - HFOV / 2))


################################################################################
################################# Plotting Stuff ###############################

# Collecting data for obstacles
Xc_1,Yc_1,Zc_1 = data_for_cylinder_along_z(x_o_1,y_o_1,obs_r,80)
Xc_2,Yc_2,Zc_2 = data_for_cylinder_along_z(x_o_2,y_o_2,obs_r,80)
Xc_3,Yc_3,Zc_3 = data_for_cylinder_along_z(x_o_3,y_o_3,obs_r,80)
Xc_4,Yc_4,Zc_4 = data_for_cylinder_along_z(x_o_4,y_o_4,obs_r,80)
Xc_5,Yc_5,Zc_5 = data_for_cylinder_along_z(x_o_5,y_o_5,obs_r,80)
Xc_6,Yc_6,Zc_6 = data_for_cylinder_along_z(x_o_6,y_o_6,obs_r,80)
Xc_7,Yc_7,Zc_7 = data_for_cylinder_along_z(x_o_7,y_o_7,obs_r,80)
Xc_8,Yc_8,Zc_8 = data_for_cylinder_along_z(x_o_8,y_o_8,obs_r,80)
Xc_9,Yc_9,Zc_9 = data_for_cylinder_along_z(x_o_9,y_o_9,obs_r,80)
Xc_10,Yc_10,Zc_10 = data_for_cylinder_along_z(x_o_10,y_o_10,obs_r,80)
x_e, y_e = ellipse(a_p, b_p, x_e_1[mpc_iter], y_e_1[mpc_iter])

# Plotting starts from here
xx1 = np.array(xx)
xs1 = np.array(ss)
ss1 = np.zeros((loop_run+1))
x_e_2 = np.array(x_e)
y_e_2 = np.array(y_e)

fig = plt.figure(figsize=(14, 9))
ax = plt.axes(projection = '3d')
ax.plot3D(xx1[0, 0:loop_run-1], xx1[1, 0:loop_run-1], xx1[2, 0:loop_run-1], linewidth = "3", color = "red")
ax.plot3D(xs1[0,0:loop_run-1], xs1[1,0:loop_run-1], ss1[0:loop_run-1], '--', linewidth = "3", color = "blue")
#ax.plot3D(xx1[0, 0:loop_run-1], xx1[1, 0:loop_run-1], ss1[0:loop_run-1], linewidth = "2", color = "green")
ax.plot3D(x_e_2[0:100, 0], y_e_2[0:100, 0] , ss1[0:100], linewidth = "2", color = "black")
ax.plot_surface(Xc_1, Yc_1, Zc_1, antialiased=False, color="blue")  # 0bstacles
ax.plot_surface(Xc_2, Yc_2, Zc_2, antialiased=False, color="blue")
ax.plot_surface(Xc_3, Yc_3, Zc_3, antialiased=False, color="blue")
ax.set_title('UAV FOLLOWS TARGET')


fig1 = mlab.figure()
mlab.clf()  # Clear the figure
mlab.plot3d(xx1[0, 0:loop_run-1], xx1[1, 0:loop_run-1], xx1[2, 0:loop_run-1], tube_radius=5, color=(0,0,1))
mlab.plot3d(xs1[0,0:loop_run-1], xs1[1,0:loop_run-1], ss1[0:loop_run-1], tube_radius=5, color=(1,0,0))
mlab.plot3d(xx1[0, 0:loop_run-1], xx1[1, 0:loop_run-1], ss1[0:loop_run-1], tube_radius=3, color=(0,.7,0))
mlab.plot3d(x_e_2[0:100, 0], y_e_2[0:100, 0] , ss1[0:100], tube_radius=3, color=(0,0,0))
mlab.mesh(Xc_1, Yc_1, Zc_1, color=(0,1,0.6))
mlab.mesh(Xc_2, Yc_2, Zc_2)
mlab.mesh(Xc_3, Yc_3, Zc_3)
mlab.mesh(Xc_4,Yc_4,Zc_4)
mlab.mesh(Xc_5,Yc_5,Zc_5)
mlab.mesh(Xc_6,Yc_6,Zc_6)
mlab.mesh(Xc_7,Yc_7,Zc_7)
mlab.mesh(Xc_8,Yc_8,Zc_8)
mlab.mesh(Xc_9,Yc_9,Zc_9)
mlab.mesh(Xc_10,Yc_10,Zc_10)
mlab.title('Tracking UAV')
s = [0, 1900, -300, 1000, 0, 170]
mlab.orientation_axes(xlabel='X-axis',ylabel='Y-axis',zlabel='Z-axis')
#mlab.outline(color=(0, 0, 0), extent=s)
mlab.axes(color=(0, 0, 0), extent=s)
#yy = np.arange(-300,1000,100)
#zz = np.arange(0,300,10)
#y0 = np.zeros_like(yy)
#z0 = np.zeros_like(zz)
#mlab.plot3d(y0,yy,y0,line_width=5,tube_radius=5)
#mlab.plot3d(z0,z0,zz,line_width=5,tube_radius=5)
#xx = yy = zz = np.arange(0,1900,100)
#yx = np.zeros_like(xx)
#mlab.plot3d(xx, yx, yx,line_width=5,tube_radius=5)
mlab.show()


# Calculating error between UAV and FOV
error = ca.DM.zeros(loop_run+1)
for i in range(loop_run):
    error[i] = ca.sqrt((x_e_1[i+1]-ss[0,i])**2 + (y_e_1[i+1]-ss[1,i])**2)

itr = loop_run
error1 = np.array(error)
sum_err = sum(error1[0:itr])
print(sum(error1[0:itr]))

fig = plt.figure()
plt.subplot(121)
plt.plot(error[0:itr], linewidth= "2", color = "red")
plt.xlabel("Iteration")
plt.ylabel("Error")
plt.subplot(122)
plt.bar('Error' , sum_err, color="red")

plt.show()
