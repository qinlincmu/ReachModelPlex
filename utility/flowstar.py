#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu
import math_tool
from interval import interval, inf, imath#
import vehicle_model
import interval_estimation
import os
import subprocess
from threading import Timer
import numpy as np

FLOWSTAR_TIMEOUT = 1  # 50ms, 0.01


def execute_flowstar(vehicle_type, ODE, horizon, state, command, waypoint_x, waypoint_y,
                    uncertainty_state, uncertainty_control,
                    stepsize, length, width):
    '''
        Function: External execution of reachability computation by calling ./RC_bicycle
        Input: type of dynamic model, state (x, y, z, yaw, vx, vy, yaw_dot),
               command, waypoint coordinate, control bounds
        Output: flowpipe points for visulization, safety indicator
    '''
    uncertainty_x = uncertainty_state[0]
    uncertainty_y = uncertainty_state[1]
    uncertainty_yaw = uncertainty_state[3]
    uncertainty_vx = uncertainty_state[4]
    uncertainty_vy = uncertainty_state[5]
    uncertainty_yawdot = uncertainty_state[6]

    uncertainty_v = uncertainty_control[0]
    uncertainty_delta = uncertainty_control[1]

    horizon = stepsize*horizon

    #uncertainty_angle = delta_bound #control steering angle uncertainty, SUE-GP
    #uncertainty_speed = speed_bound #control velocity uncertainty, SUE-GP

    rot = np.array(
        [
        [np.cos(state[3]), -np.sin(state[3])],
        [np.sin(state[3]), np.cos(state[3])]
        ]
        )
    x_dot = rot.dot(np.array([[state[3]], [state[4]]]))[0, 0]
    y_dot = rot.dot(np.array([[state[3]], [state[4]]]))[1, 0]

    pos_x = [-uncertainty_x, uncertainty_x]
    pos_x_I = interval(pos_x[0], pos_x[1])
    pos_y = [-uncertainty_y, uncertainty_y]
    pos_y_I = interval(pos_y[0], pos_y[1])
    yaw = [state[3]-uncertainty_yaw, state[3]+uncertainty_yaw]
    yaw_I = interval(yaw[0], yaw[1])
    #print(yaw)
    # x_dot = [math_tool.positive(x_dot-uncertainty_speed_state), math_tool.positive(x_dot+uncertainty_speed)]
    # x_dot_I = interval([x_dot[0], x_dot[1]])

    vx = [math_tool.positive(state[4]-uncertainty_vx), math_tool.positive(state[4]+uncertainty_vx)]
    vx_I = interval(vx[0], vx[1])

    vy = [math_tool.positive(state[5]-uncertainty_vy), math_tool.positive(state[5]+uncertainty_vy)]
    vy_I = interval(vy[0], vy[1])

    # y_dot = [math_tool.positive(y_dot-uncertainty_speed), math_tool.positive(y_dot+uncertainty_speed)]
    # y_dot_I = interval([y_dot[0], y_dot[1]])
    yaw_dot = [math_tool.positive(state[6]-uncertainty_yawdot), math_tool.positive(state[6]+uncertainty_yawdot)]
    yaw_dot_I = interval(yaw_dot[0], yaw_dot[1])
    v = [command[0]-uncertainty_v, command[0]+uncertainty_v]
    v_I = interval(v[0], v[1])
    delta = [command[1]-uncertainty_delta, command[1]+uncertainty_delta]
    delta_I = interval(delta[0], delta[1])
    
    model_class = vehicle_model.load_model(vehicle_type)
    model = model_class()
    # alpha_f = interval_estimation.initial_alpha_f_interval(model, y_dot_I, yaw_dot_I, x_dot_I, delta_I)
    # alpha_r = interval_estimation.initial_alpha_r_interval(model, y_dot_I, yaw_dot_I, x_dot_I)

    ####################for model without beta#########################
    exec_file = os.path.join(ODE, './RC_bicycle')
    if ODE == 'nonlinear_without_beta':
        args = [exec_file, str(pos_x_I[0][0]), str(pos_x_I[1][0]), str(pos_y_I[0][0]), str(pos_y_I[1][0]),
                str(yaw_I[0][0]), str(yaw_I[1][0]), str(delta_I[0][0]), str(delta_I[1][0]),
                str(v_I[0][0]), str(v_I[1][0]), str(horizon), str(stepsize), str(model.L)]
        #print (args)
        #print ('x_I', pos_x_I, pos_x_I[0][0], pos_x_I[0][1])
        # print ('y_I', pos_y_I)
        # print ('yaw_I', yaw_I)
        # print ('vx_I', vx_I)
        # print ('vy_I', vy_I)
        # print ('yaw_dot_I', yaw_dot_I)
        # print ('v_I', v_I)
        # print ('delta_I', delta_I)
    ####################for model with beta############################
    elif ODE == 'nonlinear_with_beta':
        beta_I = interval_estimation.initial_beta_interval(model, delta_I)
        co = './RC_bicycle '+str(pos_x_I[0][0]) + ' '+ str(pos_x_I[0][1]) +' ' + str(pos_y_I[0][0]) + ' '+ str(pos_y_I[0][1]) +\
             ' ' + str(yaw_I[0][0]) +' '+ str(yaw_I[0][1]) + ' ' + str(beta_I[0][0]) + ' ' +str(beta_I[0][1])+\
             ' '+str(v_I[0][0]) + ' '+str(v_I[0][1])
        os.system('cd ' + ODE +';'+co)
    elif ODE == 'linear':
        a13 = -1*command[0]*np.sin(state[3])
        a23 = command[0]*np.cos(state[3])
        b11 = np.cos(state[3])
        b21 = np.sin(state[3])
        b31 = np.tan(command[1])/0.12
        b32 = command[0]/(0.12*np.cos(command[1])*np.cos(command[1])) 
        co = './RC_bicycle '+ str(a13) +' '+str(a23) + ' ' + str(b11) + ' '+str(b21)+' '+str(b31)+' '+ str(b32) +\
            ' '+str(pos_x_I[0][0]) + ' '+ str(pos_x_I[0][1]) +' ' + str(pos_y_I[0][0]) + ' '+ str(pos_y_I[0][1]) + ' ' + str(yaw_I[0][0]) +\
            ' '+ str(yaw_I[0][1]) + ' ' +str(v_I[0][0]) + ' '+str(v_I[0][1]) + ' '+ str(delta_I[0][0]) +\
            ' '+ str(delta_I[0][1]) + ' '+str(-1*speed_bound)+' '+str(speed_bound) + ' '+str(-1*delta_bound) +' '+str(delta_bound)
        os.system('cd ' + ODE +';'+co)
    else:#'linear_tire'
        beta = interval_estimation.initial_beta_interval(model, delta_I)
        alpha_f = interval_estimation.initial_alpha_f_interval(model, y_dot_I, yaw_dot_I, x_dot_I, delta_I)
        alpha_r = interval_estimation.initial_alpha_r_interval(model, y_dot_I, yaw_dot_I, x_dot_I)
        Fxf = [0, 0] #fxi = C*kappa, 0 if no skid
        Fxr = [0, 0]
        Fyf = [min(-1*model.C_alpha*interval(alpha_f)[0][0], -1*model.C_alpha*interval(alpha_f)[0][1]), max(-1*model.C_alpha*interval(alpha_f)[0][0], -1*model.C_alpha*interval(alpha_f)[0][1])]
        Fyr = [min(-1*model.C_alpha*interval(alpha_r)[0][0], -1*model.C_alpha*interval(alpha_r)[0][1]), max(-1*model.C_alpha*interval(alpha_r)[0][0], -1*model.C_alpha*interval(alpha_r)[0][1])]
        args = [exec_file, str(pos_x[0]), str(pos_x[1]), str(pos_y[0]), str(pos_y[1]), str(yaw[0]), str(yaw[1]), str(yaw_dot[0]), str(yaw_dot[1]),
               str(Fxf[0]), str(Fxf[1]), str(Fxr[0]), str(Fxr[1]), str(Fyf[0]), str(Fyf[1]), str(Fyr[0]), str(Fyr[1]),
               str(beta[0]), str(beta[1]), str(command_speed[0]), str(command_speed[1]), str(command_delta[0]), str(command_delta[1]), 
               str(vx[0]), str(vx[1]), str(vy[0]), str(vy[1]), str(y_dot[0]), str(y_dot[1]),
               str(model.I_z), str(model.L_f), str(model.L_r), str(model.m), str(horizon), str(stepsize)]
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    timer = Timer(FLOWSTAR_TIMEOUT, process.kill)
    try:
        timer.start()
        (output, err) = process.communicate()
    finally:
        timer.cancel()
        if not timer.is_alive():
            print("Flowstar timed out at %.3f! Now running backup reachflow prediction")
            return
    output = output.decode(encoding='utf8')
    line = output.strip().split(';')
    x_bound = list()
    y_bound = list()
    yaw_bound = list()
    vx_bound = list()
    vy_bound = list()
    yawdot_bound = list()
    seg_num = len(line)-1 if line[-1]=='' else len(line)
    for i in range(seg_num):
        seg = line[0].split(',')
        x_bound.append((float(seg[0])+state[0]-0.5*length, float(seg[1])+state[0]+0.5*length)) #add up x global position by shifting
        y_bound.append((float(seg[2])+state[1]-0.5*length, float(seg[3])+state[1]+0.5*length)) #add up y global position by shifting
        yaw_bound.append((float(seg[4]), float(seg[5])))
        if ODE == 'nonlinear_without_beta':
            vx_bound.append((vx_I[0][0], vx_I[0][1]))
            vy_bound.append((vy_I[0][0], vy_I[0][1]))
            yawdot_bound.append((yaw_dot_I[0][0], yaw_dot_I[0][1]))
    return x_bound, y_bound, yaw_bound, vx_bound, vy_bound, yawdot_bound
def parse_flow(flow_x, flow_y):
    x_list, y_list = [], []
    x_list_i, y_list_i = [], []
    for i in range(len(flow_x.split(';'))):
        x_temp = flow_x.split(';')[i] #x position in one segment
        y_temp = flow_y.split(';')[i] #y position in one segment
        x_list.append(x_temp.split(',')) #list, one element is a list of x position in one segment
        y_list.append(y_temp.split(',')) #list, one element is a list of y position in one segment
            
    for i in range(len(x_list)):
        x_list_i.append([float(item) for item in x_list[i]]) #x position in one segment
        y_list_i.append([float(jtem) for jtem in y_list[i]]) #y position in one segment
    return x_list_i, y_list_i