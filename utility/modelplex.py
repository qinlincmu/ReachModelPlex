#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu
import math
import numpy as np

def normalize_angle(angle):
    """
    Normalize angle from [-pi, pi].
    """
    angle = angle % (2*np.pi)
    if (angle >= np.pi):
        angle -= 2*np.pi
    return angle


def rotate_state(state, angle):
    """
    Rotate a state [x, y, yaw, x_dot, y_dot, yaw_dot] by a specified angle
    in radians.
    """
    if len(state) != 6:
        raise ValueError('Invalid state; dimension mismatch')

    x = state[0]
    y = state[1]
    yaw = state[2]

    new_state = np.array(state)

    new_state[0] = x*np.cos(angle) - y*np.sin(angle)
    new_state[1] = y*np.cos(angle) + x*np.sin(angle)
    new_state[2] = normalize_angle(yaw + angle)

    return new_state


def translate_state(state, translation):
    """
    Translate a state [x, y, yaw, x_dot, y_dot, yaw_dot] by a translation
    [dx, dy].
    """
    if len(state) != 6:
        raise ValueError('Invalid state; dimension mismatch')
    if len(translation) != 2:
        raise ValueError('Invalid translation amount; dimension mismatch')

    state = np.array(state)
    state[0] += translation[0]
    state[1] += translation[1]
    return state
class state:
    def __init__(self, a=0, k=0, t=0, v=0, vh=0, vl=0, xg=0, yg=0):
        self.a = a
        self.k = k
        self.t = t
        self.v = v
        self.vh = vh
        self.vl =vl
        self.xg = xg
        self.yg =yg

class parameters:
    def __init__(self, A=0, B=0, T=0, eps=0):
        self.A = A
        self.B = B
        self.T  = T
        self.eps = eps

class plantparameters:
    def __init__(self, A=0, B=0, T=0, a=0, eps=0, k=0, vh=0, vl=0): 
        self.A = A
        self.B = B
        self.T = T
        self.a = a
        self.eps = eps
        self.k = k
        self.vh = vh
        self.vl = vl
        
class plantstate:
    def __init__(self,t=0, t_0=0, v=0, v_0=0, xg=0, xg_0=0, yg=0, yg_0=0):
        self.t = t
        self.t_0 = t_0
        self.v = v
        self.v_0 = v_0
        self.xg = xg
        self.xg_0 = xg_0
        self.yg = yg
        self.yg_0 = yg_0

def convertState(current_state):
    '''
        Function: converting current robot's state, the resulting state will be used for modelplex's state
        Output: converted state will be used for modelplex's state
    '''
    xr = current_state[0]
    yr = current_state[1]
    zt = current_state[3]
    xg = current_state[7]
    yg = current_state[8]
    k =  current_state[9]

    coord_trafo_state = list()
    coord_trafo_state.append(0.0)
    coord_trafo_state.append(0.0)
    coord_trafo_state.append(0.0)
    coord_trafo_state.append(0.0)
    coord_trafo_state.append(current_state[4])
    coord_trafo_state.append(current_state[5])
    coord_trafo_state.append(current_state[6])
    convxg = (xg-xr)*math.cos(zt) + (yg-yr)*math.sin(zt)
    convyg = -(xg-xr)*math.sin(zt) + (yg-yr)*math.cos(zt)
    #print("x-y: ", xr, yr)
    #print("convxg-convyg: ", convxg, convyg)
    #/*     Car           Monitor
    # *       x           y
    # *       ^           ^
    # *       |           |
    # * y <---+           +---> x
    # */
    coord_trafo_state.append(-convyg)
    coord_trafo_state.append(convxg)
    coord_trafo_state.append(-k)               

    converted_state = convertUnits(coord_trafo_state)

    return converted_state

def convertUnits(coord_trafo_state):
    '''
        Function: converting units
    '''
    converted_state = list()
    converted_state.append(coord_trafo_state[0]*10.0) #xr [dm] @note = 0
    converted_state.append(coord_trafo_state[1]*10.0) #yr [dm] @note = 0
    converted_state.append(coord_trafo_state[2]*10.0) #zr [dm] @note = 0
    converted_state.append(coord_trafo_state[3])     #orientation radians (unused in monitor)
    converted_state.append(coord_trafo_state[4]*10.0) #vx [dm/s]
    converted_state.append(coord_trafo_state[5]*10.0) #vy [dm/s] @note = 0
    converted_state.append(coord_trafo_state[6]*10.0) #vz [dm/s] @note = 0
    converted_state.append(coord_trafo_state[7]*10.0) #xg [dm]
    converted_state.append(coord_trafo_state[8]*10.0) #yg [dm]
    converted_state.append(coord_trafo_state[9]*100.0) #k  [centi-(meters^-1)]

    return converted_state

def convertAction(proposed_action):
    '''
        Function: converting proposed action to future state
        Input: robot's standard control action from ROS
        Output: converted future state
    '''
    vset = proposed_action[0]
    steer = proposed_action[1]
    
    wheelbase = 0.257
    d_com_rear = 0.1294
    
    # sanity check: is steering angle roughly planned curvature
    # convert steering angle into curve radius (bicycle model)
    # simpler alternative: double r = wheelbase/tan(steer);
    cot = 1.0/math.tan(steer) if steer != 0 else 0.0
    r = math.sqrt(d_com_rear*d_com_rear + wheelbase*wheelbase*cot*cot) if steer != 0 else 0.0
    if steer > 0:
        ksteer = 100.0/r
    elif steer < 0:
        ksteer = -100.0/r
    else:
        ksteer = 0
    
    ksteer2 = 100.0/(wheelbase*cot) if steer != 0 else 0.0
    
    # actual acceleration profile to compute acceleration from vset unavailable, so we pretend to reach vset from current v within T_CYCLE_TIME 
    a = proposed_action[2] #[dm/s^2]
    t = 0.0      #timer expected to be reset [ds]
    vh = proposed_action[3] #[dm/s]    
    vl = proposed_action[4] #[dm/s]
    xg = proposed_action[5] #[dm]
    yg = proposed_action[6] #[dm]
    k  = proposed_action[7] #curvature [centi-(meters^-1)]
    
            
    converted_action = list()
    converted_action.append(a)
    converted_action.append(k)
    converted_action.append(t)
    converted_action.append(vh)
    converted_action.append(vl)
    converted_action.append(xg)
    converted_action.append(yg)
    return converted_action

def extCtrl(pre, proposed_action):
    '''
        Function: converting proposed action and current state to future state
        Input:
            pre: current modelplex state
            proposed_action: robot's standard control action from ROS
        Output: converted future state
    '''
    ctrl_action = proposed_action[:]
    ctrl_action.append(pre.a)
    ctrl_action.append(pre.vh)
    ctrl_action.append(pre.vl)
    ctrl_action.append(pre.xg)
    ctrl_action.append(pre.yg)
    ctrl_action.append(pre.k)
    converted_action = convertAction(ctrl_action)
    
    result = state()
    result.a = converted_action[0]
    result.k = converted_action[1]
    result.t = converted_action[2]
    result.v = pre.v
    result.vh = converted_action[3]
    result.vl = converted_action[4]
    result.xg = converted_action[5]
    result.yg = converted_action[6]
    return result

def currConvert(curr, pre, current_state, current_waypoint, time):
    '''
        Function: convert robot's standard current state and waypoint to modelplex's state
        Input: 
             curr (current modelplex state without values)
             pre (previous modelplex state)
             current_state (robot's standard current state from ROS)
             current_waypoint (robot's standard current waypoints from ROS)
             time: time elapsed from last control loop to the current one
        Output: curr with asigned values
    '''
    ctrl_state = current_state[:]
    ctrl_state.extend(current_waypoint)
    converted_state = convertState(ctrl_state)
    if time > 0:
        curr.a = (converted_state[4] - float(pre.v))/time
    else:
        curr.a = 0
    curr.k = converted_state[9]
    curr.t = time*10.0
    curr.v = converted_state[4]
    curr.vh = float(pre.vh)
    curr.vl = float(pre.vl)
    curr.xg = converted_state[7]
    curr.yg = converted_state[8]
    return curr

def plantParamsConvert(plantParams, params, curr):
    '''
        Function: getting modelplex's plant parameters
        Input: 
             plantParams (plant parameters without values)
             params (modelplex's parameters)
             curr (current modelplex state)
        Output: plantParams with asigned values
    '''
    plantParams.A =float(params.A)
    plantParams.B = float(params.B)
    plantParams.T = float(params.T)
    plantParams.a = float(curr.a)
    plantParams.eps = float(params.eps)
    plantParams.k = float(curr.k)
    plantParams.vh = float(curr.vh)
    plantParams.vl = float(curr.vl)
    return plantParams

def plantcurrConvert(plantcurr, pre, curr):
    '''
        Function: getting modelplex's plant state
        Input: 
             plantcurr (current plant state without values)
             pre (previous modelplex's state)
             curr (current modelplex's state)
        Output: plantcurr with asigned values
    '''
    plantcurr.t = float(curr.t)
    plantcurr.t_0 = float(pre.t)
    plantcurr.v = float(curr.v)
    plantcurr.v_0 = float(pre.v)
    plantcurr.xg = float(curr.xg)
    plantcurr.xg_0 = float(pre.xg)
    plantcurr.yg = float(curr.yg)
    plantcurr.yg_0 = float(pre.yg)
    return plantcurr

def postConvert(post, curr, proposed_action):
    '''
        Function: getting modelplex's future state
        Input: 
             post (future state)
             curr (current modelplex's state)
             proposed_action (robot's standard control action from ROS)
        Output: post with asigned values
    '''
    ext_ctrl = extCtrl(curr, proposed_action)
    post.a = float(ext_ctrl.a)
    post.k = float(ext_ctrl.k)
    post.t = float(ext_ctrl.t)
    post.v = float(ext_ctrl.v)
    post.vh = float(ext_ctrl.vh)
    post.vl = float(ext_ctrl.vl)
    post.xg = float(ext_ctrl.xg)
    post.yg = float(ext_ctrl.yg)
    return post