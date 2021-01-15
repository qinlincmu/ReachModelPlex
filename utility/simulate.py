#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu
import math
import vehicle_model
import copy
import numpy as np

def simulate(vehicle_type, old_state, command, ODE, DT):
    '''
        Function: model-based state update
        vehicle_type: "RC_Car"/"MRZR"
        old_state: [pos_x, pos_y, pos_z, yaw, vx, vy, yaw_dot]
        command: [delta, v]
        ODE: options "nonlinear_without_beta", "linear", "nonlinear_with_beta", "linear_tire"
        DT: sampling time
    '''
    state_predict = list()
    delta = command[1]
    v = command[0]
    Model_class = vehicle_model.load_model(vehicle_type)
    model = Model_class()
    if ODE == "nonlinear_without_beta":
        x_dot = v * math.cos(old_state[3])
        x = old_state[0] + x_dot*DT
        y_dot = v * math.sin(old_state[3])
        y = old_state[1] + y_dot* DT
        rot = np.array(
               [
                 [np.cos(old_state[3]), -np.sin(old_state[3])],
                 [np.sin(old_state[3]), np.cos(old_state[3])]
               ]
              )
        rot_inv = np.linalg.inv(rot)
        vx = rot_inv.dot(np.array([[x_dot], [y_dot]]))[0, 0]
        vy = rot_inv.dot(np.array([[x_dot], [y_dot]]))[1, 0]
        yaw_dot = v / model.L * math.tan(delta)
        yaw = old_state[3] + yaw_dot * DT
        state_predict = [x, y, 0, yaw, vx, vy, yaw_dot]
    if ODE == "linear_tire":
        vx = old_state[4]
        vy = old_state[5]
        phi = old_state[3]
        alpha_f = math.atan2((vy+model.L_f*old_state[6])/vx,1)-delta if vx!= 0 else 0
        alpha_r = math.atan2((vy-model.L_f*old_state[6])/vx,1) if vx!= 0 else 0
        k = 0
        F_xf = model.C_x*k
        F_xr = model.C_x*k
        F_yf = -1*model.C_alpha*alpha_f
        F_yr = -1*model.C_alpha*alpha_r
        vx_dot = old_state[6]*vy + (F_xr-F_yf*math.sin(delta))/model.m
        vy_dot = -1*old_state[6]*vx + (F_xr*math.cos(delta)+F_yf)/model.m
        yaw_dotdot = (model.L_f*F_yf-model.L_r*F_yr)/model.I_z
        
        # v = math.sqrt(v_x*v_x+v_y*v_y)
        # beta = math.atan2(v_y/v_x, 1)
        # x_dot = v*math.cos(beta+old_state[3])# or x_dot*cos(phi)-y_dot*sin(phi)
        # y_dot = v*math.sin(beta+old_state[3])# or x_dot*sin(phi)+y_dot*cos(phi)

        #yaw_dot = v / model.L * math.tan(delta) #should update it?

        vx = old_state[4] + vx_dot * DT
        vy = old_state[5] + vy_dot * DT
        yaw_dot = old_state[6] + yaw_dotdot * DT

        x_dot =  vx*np.cos(phi)-vy*np.sin(phi)
        y_dot =  vx*np.sin(phi)+vy*np.cos(phi)

        x = old_state[0] + x_dot*DT
        y = old_state[1] + y_dot*DT
        yaw = old_state[3] + yaw_dot * DT


        state_predict = [x, y, 0, yaw, vx, vy, yaw_dot]
    return state_predict

