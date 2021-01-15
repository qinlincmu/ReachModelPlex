#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu

from interval import interval, inf, imath#

def initial_beta_interval(model, delta_I):
    '''
       Function: compute the interval of beta for the kinematic bicycle model
       Input: steering angle interval
    '''
    L = model.L
    L_r = model.L_r
    beta = (lambda delta: imath.atan(L_r * imath.tan(delta)/L))(delta_I)
    beta_min = min(beta[0][0], beta[0][1])
    beta_max = max(beta[0][0], beta[0][1])
    return [beta_min, beta_max]
def initial_alpha_f_interval(model, v_y_I, yaw_dot_I, v_x_I, delta_I):
    '''
       Function: compute the interval of beta for the kinematic bicycle model
       Input: intervals of vy, yaw_dot, vx, delta
    '''
    L_f = model.L_f
    fun = lambda v_y,yaw_dot,v_x,delta: imath.atan((v_y+L_f*yaw_dot)/v_x)-delta
    alpha_f = fun(v_y_I, yaw_dot_I, v_x_I, delta_I)
    return [alpha_f[0][0], alpha_f[0][1]]
def initial_alpha_r_interval(model, v_y_I, yaw_dot_I, v_x_I):
    '''
       Function: compute the interval of beta for the kinematic bicycle model
       Input: intervals of vy, yaw_dot, vx
    '''
    L_f = model.L_f
    fun = lambda v_y,yaw_dot,v_x: imath.atan((v_y-L_f*yaw_dot)/v_x)
    alpha_r = fun(v_y_I, yaw_dot_I, v_x_I)
    return [alpha_r[0][0], alpha_r[0][1]]