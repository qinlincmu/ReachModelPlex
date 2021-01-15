"""

simulator for AA project

author:Qin Lin (qinlin@andrew.cmu.edu)

"""
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import os
import gif
from mpmath import *
import copy
import time as tm
import datetime
from cffi import FFI
from mpl_toolkits.mplot3d import Axes3D

sys.path.append(os.path.dirname(os.path.abspath(__file__))
               + "/aa_planner/")
sys.path.append(os.path.dirname(os.path.abspath(__file__))
               + "/utility/")
try:
    from policy import Policy
    print ('CPO policy imported!')
    import vehicle_model
    import verification
    import visualize
    import flowstar
    import constant
    import modelplex
    import tubempc
except:
    raise

MAX_TIME = 5.0  # max simulation time

DT = 0.1  # [s] time tick



save_gif = True


vehicle_type = 'RC_Car'
flowstar_stepsize = 0.1


class State:
    """
    vehicle state class
    """

    def __init__(self, x=0, x_dot=0.0, y=0.0, y_dot=0.0, vx=0.0, vy=0.0, yaw=0.0, yaw_dot=0.0, v=0.0):
        self.x = x
        self.x_dot = x_dot
        self.y = y
        self.y_dot = y_dot
        self.yaw = yaw
        self.yaw_dot = yaw_dot
        self.v = v
        self.vx = vx
        self.vy = vy




def tire_dyn_r(v_x, wheel_vx, alpha):
    """
        function: tire force calculate
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()

    ##find longitudinal wheel slip K (kappa)
    if np.abs(wheel_vx-v_x) < 0.01 or (np.abs(wheel_vx) < 0.01 and np.abs(v_x) < 0.01):
        K = 0
    elif np.abs(v_x) < 0.01:   #infinite slip, longitudinal saturation
        K = math.inf
        Fx = np.sign(wheel_vx)*vehicle.mu*vehicle.load_r
        Fy = 0
        return Fx, Fy
    else:
        K = (wheel_vx-v_x)/np.abs(v_x)
    
    ###instead of avoiding -1, now look for positive equivalent
    if K < 0:
        spin_dir = -1
        K = np.abs(K)
    else:
        spin_dir = 1
    
    #alpha > pi/2 cannot be adapted to this formula
    #because of the use of tan(). Use the equivalent angle instead.
    #alpha > pi/2 means vehicle moving backwards
    #Fy sign has to be reversed, but the *sign(alpha) will take care of it
    if np.abs(alpha) > np.pi/2:
        alpha = (np.pi-np.abs(alpha))*np.sign(alpha)
 
    gamma = np.sqrt(vehicle.C_x**2 * (K/(1+K))**2 + vehicle.C_alpha**2 * (np.tan(alpha)/(1+K))**2 )
    
    if gamma <= 3*vehicle.mu*vehicle.load_r:
        F = gamma - 1/(3*vehicle.mu*vehicle.load_r)*gamma**2 + 1/(27*vehicle.mu**2*vehicle.load_r**2)*gamma**3
    else:
        #more accurate modeling with peak friction value
        F = vehicle.mu_s*vehicle.load_r
    
    if gamma == 0:
        Fx = 0
        Fy = 0
    else:
        Fx = vehicle.C_x/gamma * (K/(1+K)) * F * spin_dir
        Fy = -vehicle.C_alpha/gamma * (np.tan(alpha)/(1+K)) * F

    return Fx, Fy


def tire_dyn_f(alpha):
    """
        function: tire force calculate
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    #alpha > pi/2 cannot be adapted to this formula
    #because of the use of tan(). Use the equivalent angle instead.
    #alpha > pi/2 means vehicle moving backwards
    #Fy sign has to be reversed, but the *sign(alpha) will take care of it
    if np.abs(alpha) > pi/2:
        alpha = (np.pi-np.abs(alpha))*np.sign(alpha)
    
    alpha_sl = math.atan(3*vehicle.mu*vehicle.load_f/vehicle.C_alpha)
    if np.abs(alpha) <= alpha_sl:
        Fy = -vehicle.C_alpha*np.tan(alpha) + vehicle.C_alpha**2/(3*vehicle.mu*vehicle.load_f)*np.abs(np.tan(alpha))*np.tan(alpha) - vehicle.C_alpha**3/(27*vehicle.mu**2*vehicle.load_f**2)*np.tan(alpha)**3
    else:
        Fy = -vehicle.mu*vehicle.load_f*np.sign(alpha)
    return Fy

def update_state_ffast(state, vc, delta):
    """
        Function: state update based on dynamic bicycle model used in ffast
        Input: state and control action
        Output: return updated state 
    """

    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    if np.abs(state.vx) < 0.01 and np.abs(state.vy) <0.01:
        alpha_f = 0
        alpha_r = 0
    else:
        alpha_f = math.atan2((state.vy+vehicle.L_f*state.yaw_dot), state.vx)-delta
        alpha_r = math.atan2((state.vy-vehicle.L_r*state.yaw_dot), state.vx)
    F_yf = tire_dyn_f(alpha_f)
    F_xr, F_yr = tire_dyn_r(state.vx, vc, alpha_r)

    T_z = vehicle.L_f*F_yf*np.cos(delta)-vehicle.L_r*F_yr
    ma_x = F_xr-F_yf*np.sin(delta)
    ma_y = F_yf*np.cos(delta)+F_yr

    ####without damping
    yawdot_dot = T_z/vehicle.I_z
    vx_dot = ma_x/vehicle.m + state.yaw_dot*state.vy
    vy_dot = ma_y/vehicle.m - state.yaw_dot*state.vx


    ####with damping
    # yawdot_dot = T_z/vehicle.I_z -0.02*state.yaw_dot
    # vx_dot = ma_x/vehicle.m + state.yaw_dot*state.vy -0.025*state.vx
    # vy_dot = ma_y/vehicle.m - state.yaw_dot*state.vx -0.025*state.vy

    ###translate to inertial frame
    state.v = math.sqrt(state.vx**2+state.vy**2)
    beta = math.atan2(state.vy, state.vx)

    state.x_dot = state.v*math.cos(beta+state.yaw)
    state.y_dot = state.v*math.sin(beta+state.yaw)

    state.x = state.x + state.x_dot*DT
    state.y = state.y + state.y_dot*DT
    state.yaw = state.yaw + state.yaw_dot * DT
    state.vx = state.vx+vx_dot * DT
    state.vy = state.vy+vy_dot * DT
    state.yaw_dot = state.yaw_dot + yawdot_dot*DT

    return state
def update_state_ge(state, v_cmd, delta):
    """
        Function: state update based on numerically stable dynamic bicycle model (Ref: https://arxiv.org/abs/2011.09612)
        Input: state and control action
        Output: return updated state 
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()

    state.x = state.x + DT*(state.vx*np.cos(state.yaw)-state.vy*np.sin(state.yaw))
    state.y = state.y + DT*(state.vy*np.cos(state.yaw)+state.vx*np.sin(state.yaw))
    state.yaw = state.yaw + state.yaw_dot * DT
    state.vx = v_cmd
    state.vy = (vehicle.m*state.vx*state.vy+DT*(vehicle.L_f*vehicle.C_alpha-vehicle.L_r*vehicle.C_alpha)*state.yaw_dot-DT*vehicle.C_alpha*delta*state.vx-DT*vehicle.m*state.vx**2*state.yaw_dot)/(vehicle.m*state.vx-DT*2*vehicle.C_alpha)
    state.yaw_dot = (vehicle.I_z*state.vx*state.yaw_dot+DT*(vehicle.L_f*vehicle.C_alpha-vehicle.L_r*vehicle.C_alpha)*state.vy-DT*vehicle.L_f*vehicle.C_alpha*delta*state.vx)/(vehicle.I_z*state.vx-DT*(vehicle.L_f**2*vehicle.C_alpha+vehicle.L_r**2*vehicle.C_alpha))

    state.v = math.sqrt(state.vx**2+state.vy**2)
    beta = math.atan2(state.vy, state.vx)

    state.x_dot = state.v*math.cos(beta+state.yaw)
    state.y_dot = state.v*math.sin(beta+state.yaw)
    return state



def update_state_linear_tire(state, v_cmd, delta):
    """
        Function: state update based on dynamic bicycle model with linear tire model
        Input: state and control action
        Output: return updated state 
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    if np.abs(state.vy) <= 0.01 or np.abs(state.vx) <= 0.01:
        state = update_state_slip(state, v_cmd, delta)
        return state
    else:
        theta_f = (state.vy+vehicle.L_f*state.yaw_dot)/state.vx
        theta_r = (state.vy-vehicle.L_r*state.yaw_dot)/state.vx
        
        a22 = -4.0*vehicle.C_alpha/(vehicle.m*state.vx)
        a24 = -state.vx-(2.0*vehicle.C_alpha*vehicle.L_f-2.0*vehicle.C_alpha*vehicle.L_r)/(vehicle.m*state.vx)
        a42 = (-2.0*vehicle.L_f*vehicle.C_alpha+2.0*vehicle.L_r*vehicle.C_alpha)/(vehicle.I_z*state.vx)
        a44 = (-2.0*vehicle.L_f*vehicle.L_f*vehicle.C_alpha-2.0*vehicle.L_r*vehicle.L_r*vehicle.C_alpha)/(vehicle.I_z*state.vx)

        b2 = (2.0*vehicle.C_alpha)/vehicle.m
        b4 = (2.0*vehicle.L_f*vehicle.C_alpha)/(vehicle.I_z)

        vy_dot = state.vy*a22+state.yaw_dot*a24+delta*b2
        yaw_dotdot = state.vy*a42+state.yaw_dot*a44+delta*b4

    state.vx = v_cmd

    beta = math.atan2(state.vy, state.vx)
    state.v = math.sqrt(state.vx**2+state.vy**2)
    state.x_dot = state.v*math.cos(beta+state.yaw)
    state.y_dot = state.v*math.sin(beta+state.yaw)

    state.x = state.x + state.x_dot*DT
    state.y = state.y + state.y_dot*DT
    state.yaw = state.yaw + state.yaw_dot * DT
    
    state.vy = state.vy+vy_dot * DT
    state.yaw_dot = state.yaw_dot + yaw_dotdot*DT

    return state

def update_state(state, a, delta):
    """
        Function: state update based on kinematic bicycle model without slip
        Input: state and control action
        Output: return updated state 
    """
    
    state.x_dot = state.v * math.cos(state.yaw) 
    state.x = state.x + state.x_dot* DT
    state.y_dot = state.v * math.sin(state.yaw)
    state.y = state.y + state.y_dot* DT

    rot = np.array(
               [
                 [np.cos(state.yaw), -np.sin(state.yaw)],
                 [np.sin(state.yaw), np.cos(state.yaw)]
               ]

              )
    rot_inv = np.array(
               [
                 [np.cos(state.yaw), np.sin(state.yaw)],
                 [-np.sin(state.yaw), np.cos(state.yaw)]
               ]

              )
    
    state.vx = rot_inv.dot(np.array([[state.x_dot], [state.y_dot]]))[0, 0]
    state.vy = rot_inv.dot(np.array([[state.x_dot], [state.y_dot]]))[1, 0]

    state.yaw_dot = state.v / WB * math.tan(delta)
    state.yaw = state.yaw + state.yaw_dot * DT
    state.v = state.v + a * DT

    return state

def update_state_slip(state, v, delta):
    """
        Function: state update based on kinematic bicycle model considering slip
        Input: state and control action
        Output: return updated state 
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    beta = math.atan(math.tan(delta)*vehicle.L_r/vehicle.L)

    
    state.x_dot = v * math.cos(beta+state.yaw) 
    state.x = state.x + state.x_dot* DT
    state.y_dot = v * math.sin(beta+state.yaw)
    state.y = state.y + state.y_dot* DT


    
    state.vx = float(v*cos(beta))
    state.vy = float(v*sin(beta))
    state.v = np.sqrt(state.vx**2+state.vy**2)

    state.yaw_dot = v*math.tan(delta)*math.cos(beta) / vehicle.L
    state.yaw = state.yaw + state.yaw_dot * DT

    return state

def predict_state(state, v, delta, state_type):
    """
        Function: predict the next state based on control input
        Input: state: current state (class); v, delta: control input; state_type: "class" or "list"
        Output: return state represented in list 
    """

    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()

    beta = math.atan(math.tan(delta)*vehicle.L_r/vehicle.L)

    if state_type == "class":
        x_dot = v * math.cos(beta+state.yaw) 
        x = state.x + x_dot* DT
        y_dot = v * math.sin(beta+state.yaw)
        y = state.y + y_dot* DT
        vx = float(v*cos(beta))
        vy = float(v*sin(beta))
        yaw_dot = v*math.tan(delta)*math.cos(beta) / vehicle.L
        yaw = state.yaw + yaw_dot * DT
        return [x, y, 0, yaw, vx, vy, yaw_dot] #self.state
    else:
        x_dot = float(v * math.cos(beta+state[3]))
        x = state[0] + x_dot* DT
        y_dot = float(v * math.sin(beta+state[3]))
        y = state[1] + y_dot* DT
        vx = float(v*cos(beta))
        vy = float(v*sin(beta))
        yaw_dot = float(v*math.tan(delta)*math.cos(beta) / vehicle.L)
        yaw = state[3] + yaw_dot * DT
        return [x, y, 0, yaw, vx, vy, yaw_dot]
    
def compute_curvature(delta):
    """
        Function: compute curvature based on proposed steering angle
        Input: steering angle
        Output: curvature
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    beta = math.atan(math.tan(delta)*vehicle.L_r/vehicle.L)
    curvature = math.tan(delta)*math.cos(beta)/vehicle.L
    return curvature

@gif.frame
def plot(t, wp_x, wp_y, wp_k, state, clock, di, current_state, veri_para, x_segs_list, y_segs_list, vr_reach_segs, reach_verdict_value,
         mode, vr_monitor_state_list, vr_monitor_control_list, xg_list, yg_list, x_rec, y_rec, yaw_rec):


    fig = plt.figure(num=None, figsize=(10,15))
    ax1 = fig.add_subplot(411)
    plt.cla()
    # for stopping simulation with the esc key.
    plt.gcf().canvas.mpl_connect('key_release_event',
        lambda event: [exit(0) if event.key == 'escape' else None])
    #plt.plot(x, y, "ob", label="trajectory")
    ax1.plot(wp_x, wp_y, "xr", label="waypoints")
    visualize.plot_car(ax1, state.x, state.y, state.yaw, steer=di)
    #ax1.axis("equal")
    ax1.grid(True)
    ax1.set_title("Time[s]:" + str(round(clock, 2))
          + ", speed[km/h]:" + str(round(state.v * 3.6, 2)) + ", waypoint: (" + str(round(wp_x[-1], 2))+ str(round(wp_y[-1], 2))+ str(round(wp_k[-1], 2))+") " + "yaw: "+str(round(state.yaw, 2)))
    visualize.plot_reachflow(ax1, current_state, veri_para, x_segs_list, y_segs_list, vr_reach_segs, mode)

    ax2 = fig.add_subplot(423)
    ax2.set_title("reachability verdict value")
    ax2.plot(t, reach_verdict_value, "xr")
    ax3 = fig.add_subplot(424)
    ax3.set_title("reachability verdict value")
    ax3.plot(t, reach_verdict_value, "xr")

    ax4 = fig.add_subplot(425)
    ax4.set_title("Control verdict value")
    ax4.plot(t, vr_monitor_control_list[0], "xr")
    ax5 = fig.add_subplot(426)
    ax5.set_title("Control verdict ID")
    ax5.plot(t, vr_monitor_control_list[1], "xr")

    
    ax6 = fig.add_subplot(427, projection='3d')
    ax6.scatter(x_rec, y_rec, yaw_rec, c = 'b', marker='o')
    ax6.set_xlabel('X-axis')
    ax6.set_ylabel('Y-axis')
    ax6.set_zlabel('YAW-axis')

    # ax6 = fig.add_subplot(427)
    # ax6.set_title("xg")
    # ax6.plot(t, xg_list, "xr")
    # ax7 = fig.add_subplot(428)
    # ax7.set_title("yg")
    # ax7.plot(t, yg_list, "xr")
def sample_control(monitor, state, proposed_action, waypoint, curr, next, time, plantParamsNext, params, plantnext, post_pre, clock):
    
    delta_s = 0.3
    v_s = 0.0
    v_i = proposed_action[0] #not sample on velocity
    pos_action = []
    pos_state = []
    pos_x = []
    pos_y = []
    pos_yaw = []
    #print("clock: ", clock, arange(proposed_action[1]-delta_s, proposed_action[1]+delta_s, 0.01))
    #print("action-waypoint", proposed_action[1], waypoint)
    if v_s == 0:
        v_container = [proposed_action[0]]
    else:
        v_container = arange(proposed_action[0]-v_s, proposed_action[0]+v_s, 0.1)
    start = datetime.datetime.now()
    for delta_i in arange(proposed_action[1]-delta_s, proposed_action[1]+delta_s, 0.01):
        for v_i in v_container:
            next_state = predict_state(state, v_i, delta_i, "list")
            #next_state = state
            #next_waypoint = waypoint
            #print(next_waypoint)
            next_waypoint = [waypoint[0], waypoint[1], compute_curvature(delta_i)]
            
            #next = modelplex.currConvert(next, curr, next_state, next_waypoint, time)
            next = modelplex.currConvert(next, curr, state, next_waypoint, time)
            plantParamsNext = modelplex.plantParamsConvert(plantParamsNext, params, next)
            plantnext = modelplex.plantcurrConvert(plantnext, curr, next)
            post_pre = modelplex.postConvert(post_pre, next, [v_i, delta_i])
            ctrl_monitor_verdict = monitor.boundaryDist(next[0], post_pre[0], params)
            
            ctrl_verdict_value = float(ctrl_monitor_verdict.val)
            ctrl_verdict_id = float(ctrl_monitor_verdict.id)
            #print("next_wp", next_waypoint, ctrl_verdict_value, ctrl_verdict_id, next_state[1])
            if float(ctrl_verdict_value) > 0:
                ctrl_verdict_value = 1.0 #for better visualization
            else:
                ctrl_verdict_value, ctrl_verdict_id = check_verdict(state[0], state[1], waypoint[0], waypoint[1], ctrl_verdict_value, ctrl_verdict_id, clock)
            #print(float(ctrl_monitor_verdict.val))
            if ctrl_verdict_value >0:
                pos_action.append([v_i, delta_i])
                pos_state.append(next_state)
                pos_x.append(next_state[0])
                pos_y.append(next_state[1])
                pos_yaw.append(next_state[3])
    end = datetime.datetime.now()
    dt = (end-start).total_seconds()*1000
    dis = 1000
    print("sampling runtime: %5.2f ms\n" %dt)
    if len(pos_action) == 0:
        print('clock: %5.2f warning!!!!!!!!! No possible fallback exists\n' %(clock))
        return proposed_action, predict_state(state, proposed_action[0], proposed_action[1], "list"), [], [], []
    else:
        print("clock: %5.2f found %d recovery states\n" %(clock, len(pos_action)))
    for i in range(len(pos_action)):
        if abs(proposed_action[1]-pos_action[i][1]) <= dis:
            dis = abs(proposed_action[1]-pos_action[i][1])
            ind = i
    #print (pos_x)
    #print (pos_y)
    #print (pos_yaw)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos_x, pos_y, pos_yaw, c = 'b', marker='o')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('yaw')
    ax.set_title("Recovery states sampled from ModelPlex monitor")

    fig.savefig(str(clock)+'.png')
    return pos_action[ind], pos_state[ind], pos_x, pos_y, pos_yaw
def check_verdict(x, y, waypoint_x, waypoint_y, value, id, clock):
    """
        Function: regulate verdict value, when the vehicle's position is very close to the waypoint, 
                  monitor will have "false alarm" because the predicted position will pass the waypoint
        Input: vehicle's current position, waypoint, verdict value & id from monitor, timestamp
        Output: regulated verdict value & id
    """
    if id == -15 or id == -16 or clock == 0:
        if np.sqrt((x-waypoint_x)**2+(y-waypoint_y)**2) < 0.2:
            return 1, 1
        if clock == 0: #remove the false alarm in the beginning
            return 1, 1
    return value, id
def do_simulation(initial_state):
    """
    Simulation
    """
    Model_class = vehicle_model.load_model("RC_Car")
    vehicle = Model_class()
    state = initial_state
    ODE = "nonlinear_without_beta" #"nonlinear_without_beta", "linear", "nonlinear_with_beta", "linear_tire"
    uncertainty_estimation_method  = 'gp' #option['gp', 'sampling']SUE-GP and sampling distribution approach
    mode = "rounded_square" #["straight", "circle", "rounded_square"],
    clock = 0.0
    x_list = [state.x]
    x_dot_list = [state.x_dot]
    y_list = [state.y]
    y_dot_list = [state.y_dot]
    z_list = [0.0]
    yaw_list = [state.yaw]
    yaw_dot_list = [state.yaw_dot]
    v_list = [state.v]
    t_list = []
    d_list = []
    a_list = []
    wp_x_list = []
    wp_y_list = []
    wp_cur_list = []
    timestamp_list = []
    vc_list = []
    xg_list = []
    yg_list = []
    pl = Policy(planner_mode="rounded_square")#straight, rounded_square, circle

    ffi = FFI()
    monitor = ffi.dlopen(os.getcwd()+"/utility/libmonitor.so")
    ffi.cdef("""
        typedef struct plantparameters {
            long double A;
            long double B;
            long double T;
            long double a;
            long double eps;
            long double k;
            long double vh;
            long double vl;
        } plantparameters;

        typedef struct plantstate {
            long double t;
            long double t_0;
            long double v;
            long double v_0;
            long double xg;
            long double xg_0;
            long double yg;
            long double yg_0;
        } plantstate;

        typedef struct parameters {
            long double A;
            long double B;
            long double T;
            long double eps;
        } parameters;

        typedef struct state {
            long double a;
            long double k;
            long double t;
            long double v;
            long double vh;
            long double vl;
            long double xg;
            long double yg;
        } state;

        typedef struct input input;
        typedef struct verdict { int id; long double val; } verdict;
        verdict plantBoundaryDist(plantstate pre, plantstate curr, const plantparameters* const params);
        verdict boundaryDist(state pre, state curr, const parameters* const params);
    """)

    plantParams = ffi.new("struct plantparameters*")
    plantParamsNext = ffi.new("struct plantparameters*")
    params = ffi.new("struct parameters*")
    plantpre = ffi.new("struct plantstate*")
    plantcurr = ffi.new("struct plantstate*")
    plantnext = ffi.new("struct plantstate*")
    pre = ffi.new("struct state*")
    curr = ffi.new("struct state*")
    next = ffi.new("struct state*") #future state from prediction
    post = ffi.new("struct state*") #future state for keymera approach
    post_pre = ffi.new("struct state*") #future state for keymera approach adding model-based prediction of state
    #post_sim = ffi.new("struct state*") #future state for keymera approach


    last_state_timestamp = 1
    ifplot = False #False
    total_t = 0 #total running time, what unit?
    counter = 0 #interation number
    horizon = 1 #feature step in integer
    max_t = 0   #max runtime in each iteration
    min_t = 100 #min runtime in each iteration
    vr_reach = 0
    vr_monitor = 0
    vr_reach_list = []
    vr_monitor_plant_list = [[], []] #value&ID
    vr_monitor_control_list = [[], []] #value&ID

    ctrl_verdict_value = 0

    ###############initialization############# 
    params.A = constant.MAX_MOTOR_ACCEL
    params.B = constant.MAX_MOTOR_ACCEL
    params.T = constant.T_CYCLE_TIME
    params.eps = 0.5

    curr.a = 0.0
    curr.k = 0.0
    curr.t = 0.0
    curr.v = 0.0
    curr.vh = 10.0
    curr.vl = 0.0
    curr.xg = 0.0
    curr.yg = 0.0

    plantcurr.t = 0.0
    plantcurr.t_0 = 0.0
    plantcurr.v = 0.0
    plantcurr.v_0 = 0.0
    plantcurr.xg = 0.0
    plantcurr.xg_0 = 0.0
    plantcurr.yg = 0.0
    plantcurr.yg_0 = 0.0
        


    counter_timeout = 0
    bad_position_x = []
    bad_position_y = []
    counter_inside = 0 #how many future state is actually included by reachability 
    safety_violated = False #set true if unsafe
    uncertainty_estimation_method = 'gp' #option['gp', 'sampling']SUE-GP and sampling distribution approach
    x_list_buffer, y_list_buffer = [], []
      
    counter_eg_x, counter_eg_y = [], []
    counter_flow_x, counter_flow_y = [], []
    x_all, y_all, ind_all, flow_x_all, flow_y_all = [], [], [], [], []
    rec_x_all, rec_y_all, rec_yaw_all = [], [], []
    delta_t_all = []

    reach_verdict_value_all = []

    verify_result = [[], []]


    frames = []
    start_sec = tm.time()
    bad_cnt = 0
    data_f = open ("datalog.txt", "w")


    pre_ctrl_verdict_value = 1
    prev_action = [0, 0]
    state_buffer = []
    control_buffer = []
    cnt = 0
    while MAX_TIME >= clock:

        current_state_timestamp = clock
        x0 = [state.x, state.y, state.yaw, state.vx, state.vy, state.yaw_dot]  # current state

        action, waypoint, curvature = pl.get_action(x0) #get from aa_planner, [velocity, angle]

        current_state = [state.x, state.y, 0, state.yaw, state.vx, state.vy, state.yaw_dot] #self.state
        

        proposed_action = action[:]
        current_waypoint = [waypoint[0], waypoint[1], curvature]
        data_f.write(str(state.x)+" "+str(state.y)+" "+str(0)+" "+str(state.yaw)+" "+str(state.vx)+" "+str(state.vy)+" "+str(state.yaw_dot)+" ")
        data_f.write(str(current_waypoint[0])+" "+str(current_waypoint[1])+" "+str(current_waypoint[2])+" ")
        data_f.write(str(proposed_action[0])+" "+str(proposed_action[1]))
        data_f.write("\n")
        #####################modelplex verification#################
        #pre = curr

        pre.a = float(curr.a)
        pre.k = float(curr.k)
        pre.t = float(curr.t)
        pre.v = float(curr.v)
        pre.vh = float(curr.vh)
        pre.vl = float(curr.vl)
        pre.xg = float(curr.xg)
        pre.yg = float(curr.yg)
        
        plantpre.t = float(plantcurr.t)
        plantpre.t_0 = float(plantcurr.t_0)
        plantpre.v = float(plantcurr.v)
        plantpre.v_0 = float(plantcurr.v_0)
        plantpre.xg = float(plantcurr.xg)
        plantpre.xg_0 = float(plantcurr.xg_0)
        plantpre.yg = float(plantcurr.yg)
        plantpre.yg_0 = float(plantcurr.yg_0)

        time = DT
        curr = modelplex.currConvert(curr, pre, current_state, current_waypoint, time)
        
        plantParams = modelplex.plantParamsConvert(plantParams, params, curr)
        plantcurr = modelplex.plantcurrConvert(plantcurr, pre, curr)
        

        ########################simulation of bad control
        inject_bad_control = True
        if inject_bad_control == True:
            bad_ctr_sequence = [[0.7, -0.7], [0.7, -0.7], [0.7, -0.7]]
            if(state.x > 1.0 and cnt < len(bad_ctr_sequence)):
                proposed_action[0] = bad_ctr_sequence[cnt][0]
                proposed_action[1] = bad_ctr_sequence[cnt][1]
                print("clock %5.2f:, inject bad control (%5.2f, %5.2f)" %(clock, bad_ctr_sequence[cnt][0], bad_ctr_sequence[cnt][1]))
                print("action got from CPO: (%5.2f, %5.2f)" %(action[0], action[1]))
                print("x-y position: (%5.2f, %5.2f); waypoint: (%5.2f, %5.2f)\n" %(current_state[0], current_state[1], current_waypoint[0], current_waypoint[1]))
                cnt += 1
        ################################
        #next_state = predict_state(state, proposed_action[0], proposed_action[1], "class")
        next_state = current_state
        #next_waypoint = current_waypoint #less sensitive to bad steering
        next_waypoint = [waypoint[0], waypoint[1], compute_curvature(proposed_action[1])] #modify the waypoint's curvature, more sensitive to bad steering
        next = modelplex.currConvert(next, curr, next_state, next_waypoint, time)
        plantParamsNext = modelplex.plantParamsConvert(plantParamsNext, params, next)
        plantnext = modelplex.plantcurrConvert(plantnext, curr, next)
        
        ##################state verification using monitor##################
        plant_verdict = monitor.plantBoundaryDist(plantpre[0], plantcurr[0], plantParams)
        plant_verdict_value = float(plant_verdict.val)
        plant_verdict_id = int(plant_verdict.id)
        if plant_verdict_value > 0:
            plant_verdict_value = 1.0
        else:
            plant_verdict_value, plant_verdict_id = check_verdict(state.x, state.y, current_waypoint[0], current_waypoint[1], plant_verdict_value, plant_verdict_id, clock)
        #####################################################################

        ##################control verification using monitor#################
        post = modelplex.postConvert(post, curr, proposed_action)
        post_pre = modelplex.postConvert(post_pre, next, proposed_action)
        
        #ctrl_monitor_verdict = monitor.boundaryDist(curr[0], post[0], params) #original ModelPlex control checking
        ctrl_monitor_verdict = monitor.boundaryDist(next[0], post_pre[0], params)
        ctrl_verdict_value = float(ctrl_monitor_verdict.val)
        ctrl_verdict_id = float(ctrl_monitor_verdict.id)

        if float(ctrl_verdict_value) > 0:
            ctrl_verdict_value = 1.0 #for better visualization
        else:
            ctrl_verdict_value, ctrl_verdict_id = check_verdict(state.x, state.y, current_waypoint[0], current_waypoint[1], ctrl_verdict_value, ctrl_verdict_id, clock)
        ###################################################################

        last_state_timestamp = current_state_timestamp

        ##################control verification using reachability##############
        speed_bound = 0.01
        delta_bound = 0.01
        uncertainty_control = [speed_bound, delta_bound]
        uncertainty_state = [0.01, 0.01, 0, 0.01, 0.01, 0.01, 0.01] #(x, y, z, yaw, vx, vy, yaw_dot)
        start = datetime.datetime.now()
        
        x_bound, y_bound, yaw_bound, vx_bound, vy_bound, yawdot_bound = flowstar.execute_flowstar(vehicle_type, ODE, horizon, current_state,
                          proposed_action, waypoint[0], waypoint[1], uncertainty_state, uncertainty_control, flowstar_stepsize, vehicle.L, vehicle.width)

        #####apply full brake
        last2_state_x = 0.5*(x_bound[-1][0]+x_bound[-1][1])
        last2_state_y = 0.5*(y_bound[-1][0]+y_bound[-1][1])
        last2_state_yaw = 0.5*(yaw_bound[-1][0]+yaw_bound[-1][1])
        last2_state_vx = 0.5*(vx_bound[-1][0]+vx_bound[-1][1])
        last2_state_vy = 0.5*(vy_bound[-1][0]+vy_bound[-1][1])
        last2_state_yawdot = 0.5*(yawdot_bound[-1][0]+yawdot_bound[-1][1])
        last2_state = [last2_state_x, last2_state_y, 0, last2_state_yaw, last2_state_vx, last2_state_vy, last2_state_yawdot]
        
        last2_state_x_uc = 0.5*np.abs(x_bound[-1][1]-x_bound[-1][0]-vehicle.L)
        last2_state_y_uc = 0.5*np.abs(y_bound[-1][1]-y_bound[-1][0]-vehicle.width)
        last2_state_yaw_uc = 0.5*np.abs(yaw_bound[-1][1]-yaw_bound[-1][0])
        last2_state_vx_uc = 0.5*np.abs(vx_bound[-1][1]-vx_bound[-1][0])
        last2_state_vy_uc = 0.5*np.abs(vy_bound[-1][1]-vy_bound[-1][0])
        last2_state_yawdot_uc = 0.5*np.abs(yawdot_bound[-1][1]-yawdot_bound[-1][0])
        last2_state_uc = [last2_state_x_uc, last2_state_y_uc, [], last2_state_yaw_uc, last2_state_vx_uc, last2_state_vy_uc, last2_state_yawdot_uc]
        last_proposed_action = [proposed_action[0]/2, 0]#[proposed_action[0]/2, 0] #TODO, minimum instaneously deceleration?
        uncertainty_control = [0.01, 0.01]

        x_bound_last, y_bound_last, yaw_bound_last, vx_bound_last, vy_bound_last, yawdot_bound_last = flowstar.execute_flowstar(vehicle_type, ODE, horizon, last2_state,
                          last_proposed_action, waypoint[0], waypoint[1], last2_state_uc, uncertainty_control, flowstar_stepsize, vehicle.L, vehicle.width)
        end = datetime.datetime.now()
        dt = (end-start).total_seconds()*1000
        #print("Reachability runtime %5.2f" %dt)
        x_segs = []
        y_segs = []
        for i in range(len(x_bound)):
            x_segs.append([x_bound[i][0], x_bound[i][1], x_bound[i][1], x_bound[i][0]])
            y_segs.append([y_bound[i][0], y_bound[i][0], y_bound[i][1], y_bound[i][1]])
        
        #append last reachable set of full brake
        for i in range(len(x_bound_last)):
            x_segs.append([x_bound_last[i][0], x_bound_last[i][1], x_bound_last[i][1], x_bound_last[i][0]])
            y_segs.append([y_bound_last[i][0], y_bound_last[i][0], y_bound_last[i][1], y_bound_last[i][1]])

        veri_para = verification.Verify_para()
        vr_reach_segs = list() #reachability result for each flow segment


        vr_reach = verification.verify_reach(x_segs, y_segs, current_state, current_waypoint) #liveness-like verification using reachability
        for i in range(len(x_segs)):
            vr_collision = verification.verify_dubin(veri_para, mode, x_segs[i], y_segs[i], current_state, current_waypoint)
            #vr_reach_segs.append(vr_collision or vr_reach) #reachability verifies "reach-and-avoid"
            vr_reach_segs.append(vr_collision) #reachability verifies only "avoid"
        if True in vr_reach_segs:
            reach_verdict_value = -1
        else:
            reach_verdict_value = 1
        
        pc = 0.01
        rec_x, rec_y, rec_yaw = [], [], []
        if ctrl_verdict_value < 0 or plant_verdict_value <0:
            print("DETECT: x-y (%5.2f, %5.2f), waypoint (%5.2f, %5.2f)" %(state.x, state.y, current_waypoint[0], current_waypoint[1]))
            print ("clock: %5.2f, reach verdict %d, control (%5.2f, %d),  plant (%5.2f, %d)\n" %(clock, reach_verdict_value, ctrl_verdict_value, ctrl_verdict_id, plant_verdict_value, plant_verdict_id))
        fallback = True
        if fallback == True:
            if pre_ctrl_verdict_value > 0 and ctrl_verdict_value < 0: #modelPlex violation just starts
                op_control, op_state, rec_x, rec_y, rec_yaw = sample_control(monitor, current_state, prev_action, current_waypoint, curr, next, time, plantParamsNext, params, plantnext, post_pre, clock)
                control_buffer.append(op_control)
                state_buffer.append(op_state)
            if pre_ctrl_verdict_value<0 and ctrl_verdict_value < 0 and reach_verdict_value>0: #small modelPlex violation continues
                op_control, op_state, rec_x, rec_y, rec_yaw = sample_control(monitor, state_buffer[-1], control_buffer[-1], current_waypoint, curr, next, time, plantParamsNext, params, plantnext, post_pre, clock)
                control_buffer.append(op_control)
                state_buffer.append(op_state)
            if ctrl_verdict_value < 0 and reach_verdict_value<0: #large modelPlex violation
                op_control, op_state, rec_x, rec_y, rec_yaw = sample_control(monitor, state_buffer[-1], control_buffer[-1], current_waypoint, curr, next, time, plantParamsNext, params, plantnext, post_pre, clock)
                control_buffer.append(op_control)
                state_buffer.append(op_state)

                bk_state = state_buffer[-1]
                tgt_state = predict_state(bk_state, state_buffer[-1][0], state_buffer[-1][1], "list")
                tgt_x = tgt_state[0]
                tgt_y = tgt_state[1]
                tgt_v = np.sqrt(tgt_state[4]**2+tgt_state[5]**2)
                tgt_yaw = tgt_state[3]

                X_int = np.array(
                        [
                            [state.x],
                            [state.y],
                            [state.v],
                            [state.yaw]
                          ]
                     )
                X_ref = np.zeros((2, 4, 1))
                X_ref[0,0,0] = tgt_x
                X_ref[0,1,0] = tgt_y
                X_ref[0,2,0] = tgt_v
                X_ref[0,3,0] = tgt_yaw

                X_ref[1,0,0] = tgt_x
                X_ref[1,1,0] = tgt_y
                X_ref[1,2,0] = tgt_v
                X_ref[1,3,0] = tgt_yaw
                start = datetime.datetime.now()
                mpc_a, mpc_d = tubempc.mpc(X_int, X_ref)
                proposed_action[0] = state.v+mpc_a
                proposed_action[1] = mpc_d
                print("clock: %5.2f,  mpc control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (acc: %5.2f, angle: %5.2f)" %(clock, mpc_a, mpc_d))
                end = datetime.datetime.now()
                dt = (end-start).total_seconds()*1000
                print("Tube MPC runtime %5.2f" %dt)

        prev_action = [proposed_action[0], proposed_action[1]]
        
        pre_ctrl_verdict_value = ctrl_verdict_value
        
        di = proposed_action[1]
        ai = proposed_action[0]-state.v


        reach_verdict_value_all.append(reach_verdict_value)
        rec_x_all.append(rec_x)
        rec_y_all.append(rec_y)
        rec_yaw_all.append(rec_yaw)

        x_list.append(state.x)
        x_dot_list.append(state.x_dot)
        y_list.append(state.y)
        y_dot_list.append(state.y_dot)
        z_list.append(0.0)
        yaw_list.append(state.yaw)
        yaw_dot_list.append(state.yaw_dot)
        v_list.append(state.v)
        t_list.append(clock)
        d_list.append(di)
        a_list.append(ai)
        vc_list.append(action[0])
        wp_x_list.append(waypoint[0])
        wp_y_list.append(waypoint[1])
        wp_cur_list.append(curvature)
        timestamp_list.append(start_sec+clock)
        xg_list.append(curr.xg)
        yg_list.append(curr.yg)

        #vr_monitor_plant_list[0].append(1 if plant_verdict_value>0 else plant_verdict_value)
        vr_monitor_plant_list[0].append(plant_verdict_value)
        vr_monitor_plant_list[1].append(plant_verdict_id)
        #vr_monitor_control_list[0].append(1 if ctrl_verdict_value>0 else ctrl_verdict_value)
        vr_monitor_control_list[0].append(ctrl_verdict_value)
        vr_monitor_control_list[1].append(ctrl_verdict_id)
        if save_gif:
            frame = plot(t_list, wp_x_list, wp_y_list, wp_cur_list, state, clock, di, current_state, veri_para, x_segs, y_segs, vr_reach_segs, reach_verdict_value_all,
                        mode, vr_monitor_plant_list, vr_monitor_control_list, xg_list, yg_list, rec_x, rec_y, rec_yaw)
            frames.append(frame)
        #state = update_state(state, ai, di)
        state = update_state_slip(state, proposed_action[0], proposed_action[1])
        #state = update_state_ffast(state, proposed_action[0], proposed_action[1])
        #state = update_state_linear_tire(state, proposed_action[0], proposed_action[1])
        #state = update_state_ge(state, proposed_action[0], proposed_action[1])
        clock = clock + DT

    

    if save_gif:
        gif.save(frames, "simulation.gif", duration=100)
    return t_list, x_list, y_list, yaw_list, v_list, d_list, a_list, reach_verdict_value_all, vr_monitor_control_list

def main():
    print(__file__ + " start!!")

    initial_state = State(x=0.0, x_dot=0.0, y=0.0, y_dot=0.0, vx=0.0, vy=0.0, yaw=0.0, yaw_dot=0.0, v=0.0)

    
    t, x, y, yaw, v, d, a, reach, monitor_ctl = do_simulation(initial_state)
    
    st = int(1.5/DT)
    et = int(2.5/DT)
    plt.subplots()
    plt.plot(x[st:et], y[st:et], "-r")
    plt.grid(True)
    plt.xlabel("x position [m]")
    plt.ylabel("y position [m]")
    plt.title("x y position")
    plt.legend()

    # plt.subplots()
    # plt.plot(t[st:et], y[st:et], "-r", label="y")
    # plt.grid(True)
    # plt.xlabel("Time [s]")
    # plt.ylabel("y position [m]")
    # plt.legend()

    plt.subplots()
    plt.plot(t[st:et], np.array(yaw[st:et])*1, "-r", label="yaw")
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.ylabel('yaw rad $10^{-3}$')
    plt.legend()

    plt.subplots()
    plt.plot(t[st:et], np.array(d[st:et])*1, "-r", label="steering angle")
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.ylabel('Steering angle [rad]')
    plt.title("Steering angle proposed by NN controller")
    plt.legend()


    plt.subplots()
    plt.plot(t[st:et], monitor_ctl[0][st:et], "-b", label="ModelPlex control verdict")
    plt.plot(t[st:et], reach[st:et], "-r", label="Reachability verdict")
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.ylabel('Verdict value')
    #plt.title("Steering angle proposed by NN controller")
    plt.legend()

    #print(yaw)
    #print(y)
    plt.show()



if __name__ == '__main__':
    main()