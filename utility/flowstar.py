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

FLOWSTAR_TIMEOUT = 10  # 50ms, 0.01

def add_sign(number):
    if number >= 0:
        return str(number)
    else:
        return "(-1)*"+str(-1*number)

def execute_flowstar(DT, vehicle_type, ODE, horizon, state, command, waypoint_x, waypoint_y,
                    uncertainty_state, uncertainty_control,
                    stepsize, length, width, full_brake):
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

    vx = [state[4]-uncertainty_vx, state[4]+uncertainty_vx]
    vx_I = interval(vx[0], vx[1])

    vy = [state[5]-uncertainty_vy, state[5]+uncertainty_vy]
    vy_I = interval(vy[0], vy[1])

    yaw_dot = [state[6]-uncertainty_yawdot, state[6]+uncertainty_yawdot]
    yaw_dot_I = interval(yaw_dot[0], yaw_dot[1])
    v = [command[0]-uncertainty_v, command[0]+uncertainty_v]
    v_I = interval(v[0], v[1])
    delta = [command[1]-uncertainty_delta, command[1]+uncertainty_delta]
    delta_I = interval(delta[0], delta[1])
    if full_brake == False:
        min_a = min((v[0]-state[4])/DT, (v[1]-state[4])/DT)
        max_a = max((v[0]-state[4])/DT, (v[1]-state[4])/DT)
        a_I = interval(min_a, max_a)
    else:
        a_I = interval(-8.01, -7.99)#interval(-10.01, -9.99) #TODO: maximum acc/dec set to 2 in rc-car
    model_class = vehicle_model.load_model(vehicle_type)
    vehicle = model_class()

    ####################for model without beta#########################
    exec_file = os.path.join(ODE, './RC_bicycle')
    if ODE == 'nonlinear_without_beta':
        args = [exec_file, str(pos_x_I[0][0]), str(pos_x_I[1][0]), str(pos_y_I[0][0]), str(pos_y_I[1][0]),
                str(yaw_I[0][0]), str(yaw_I[1][0]), str(delta_I[0][0]), str(delta_I[1][0]),
                str(v_I[0][0]), str(v_I[1][0]), str(a_I[0][0]), str(a_I[1][0]), str(horizon), str(stepsize), str(vehicle.L)]
    ####################for model with beta############################
    elif ODE == 'nonlinear_with_beta':
        beta_I = interval_estimation.initial_beta_interval(vehicle, delta_I)
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
        a22 = -4.0*vehicle.C_alpha/(vehicle.m*command[0])
        a24 = -command[0]-(2.0*vehicle.C_alpha*vehicle.L_f-2.0*vehicle.C_alpha*vehicle.L_r)/(vehicle.m*command[0])
        a42 = (-2.0*vehicle.L_f*vehicle.C_alpha+2.0*vehicle.L_r*vehicle.C_alpha)/(vehicle.I_z*command[0])
        a44 = (-2.0*vehicle.L_f*vehicle.L_f*vehicle.C_alpha-2.0*vehicle.L_r*vehicle.L_r*vehicle.C_alpha)/(vehicle.I_z*command[0])

        b2 = (2.0*vehicle.C_alpha)/vehicle.m
        b4 = (2.0*vehicle.L_f*vehicle.C_alpha)/(vehicle.I_z)

        a24_1 =  add_sign(a24)
        a44_1 = add_sign(a44)
        b2_1 = add_sign(b2)
        b4_1 = add_sign(b4)
        # print('x', str(pos_x_I[0][0]), str(pos_x_I[1][0]))
        # print('y', str(pos_y_I[0][0]), str(pos_y_I[1][0]))
        # print('vx', str(vx_I[0][0]), str(vx_I[1][0]))
        # print('vy', str(vy_I[0][0]), str(vy_I[1][0]))
        # print('yaw', str(yaw_I[0][0]), str(yaw_I[1][0]))
        # print(state[6]-uncertainty_yawdot, state[6]+uncertainty_yawdot)
        # print(yaw_dot_I)
        # print(delta_I)
        # print('dyaw', str(yaw_dot_I[0][0]), str(yaw_dot_I[1][0]))
        # print('delta', str(delta_I[0][0]), str(delta_I[1][0]))

        args = [exec_file, str(pos_x_I[0][0]), str(pos_x_I[1][0]), str(pos_y_I[0][0]), str(pos_y_I[1][0]),
               str(vx_I[0][0]), str(vx_I[1][0]), str(vy_I[0][0]), str(vy_I[1][0]), str(yaw_I[0][0]), str(yaw_I[1][0]),
               str(yaw_dot_I[0][0]), str(yaw_dot_I[1][0]), str(a_I[0][0]), str(a_I[1][0]), str(delta_I[0][0]), str(delta_I[1][0]),
               str(a22), add_sign(a24), str(a42), add_sign(a44), add_sign(b2), add_sign(b4), str(horizon), str(stepsize)]
        #print(args[1:])
        #print(str(a22)+"*vy+"+a24_1+"*dpsi+"+b2_1+"*delta")
        #print(str(a42)+"*vy+"+a44_1+"*dpsi+"+b4_1+"*delta")
        # beta = interval_estimation.initial_beta_interval(model, delta_I)
        # alpha_f = interval_estimation.initial_alpha_f_interval(model, y_dot_I, yaw_dot_I, x_dot_I, delta_I)
        # alpha_r = interval_estimation.initial_alpha_r_interval(model, y_dot_I, yaw_dot_I, x_dot_I)
        # Fxf = [0, 0] #fxi = C*kappa, 0 if no skid
        # Fxr = [0, 0]
        # Fyf = [min(-1*model.C_alpha*interval(alpha_f)[0][0], -1*model.C_alpha*interval(alpha_f)[0][1]), max(-1*model.C_alpha*interval(alpha_f)[0][0], -1*model.C_alpha*interval(alpha_f)[0][1])]
        # Fyr = [min(-1*model.C_alpha*interval(alpha_r)[0][0], -1*model.C_alpha*interval(alpha_r)[0][1]), max(-1*model.C_alpha*interval(alpha_r)[0][0], -1*model.C_alpha*interval(alpha_r)[0][1])]
        # args = [exec_file, str(pos_x[0]), str(pos_x[1]), str(pos_y[0]), str(pos_y[1]), str(yaw[0]), str(yaw[1]), str(yaw_dot[0]), str(yaw_dot[1]),
        #        str(Fxf[0]), str(Fxf[1]), str(Fxr[0]), str(Fxr[1]), str(Fyf[0]), str(Fyf[1]), str(Fyr[0]), str(Fyr[1]),
        #        str(beta[0]), str(beta[1]), str(command_speed[0]), str(command_speed[1]), str(command_delta[0]), str(command_delta[1]), 
        #        str(vx[0]), str(vx[1]), str(vy[0]), str(vy[1]), str(y_dot[0]), str(y_dot[1]),
        #        str(model.I_z), str(model.L_f), str(model.L_r), str(model.m), str(horizon), str(stepsize)]
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
    #print("line", line)
    x_bound = list()
    y_bound = list()
    yaw_bound = list()
    vx_bound = list()
    vy_bound = list()
    yawdot_bound = list()
    seg_num = len(line)-1 if line[-1]=='' else len(line)
    #print("seg", seg_num)
    if ODE == 'nonlinear_without_beta':
        for i in range(seg_num):
            seg = line[0].split(',')
            x_bound.append((float(seg[0])+state[0]-0.5*length, float(seg[1])+state[0]+0.5*length)) #add up x global position by shifting
            y_bound.append((float(seg[2])+state[1]-0.5*length, float(seg[3])+state[1]+0.5*length)) #add up y global position by shifting
            yaw_bound.append((float(seg[4]), float(seg[5])))
            # kinematic bicycle model doesn't compute vx, vy, yawdot 
            vx_bound.append((vx_I[0][0], vx_I[0][1]))
            vy_bound.append((vy_I[0][0], vy_I[0][1]))
            yawdot_bound.append((yaw_dot_I[0][0], yaw_dot_I[0][1]))
    if ODE == 'linear_tire':
        for i in range(seg_num):
            seg = line[i].split(',')
            #print(len(seg))
            x_bound.append((float(seg[0])+state[0]-0.5*length, float(seg[1])+state[0]+0.5*length)) #add up x global position by shifting
            y_bound.append((float(seg[2])+state[1]-0.5*length, float(seg[3])+state[1]+0.5*length)) #add up y global position by shifting
            vx_bound.append((float(seg[4]), float(seg[5])))
            vy_bound.append((float(seg[6]), float(seg[7])))
            yaw_bound.append((float(seg[8]), float(seg[9])))
            yawdot_bound.append((float(seg[10]), float(seg[11])))
    #print('xxbb', x_bound)
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