#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu

import numpy as np
import math
import matplotlib.pyplot as plt
import math_tool

class Verify_para:
    def __init__(self, INNER_LINE=0.35, OUTER_LINE = 0.35, INNER_R = 0.65, OUTER_R = 1.35):
    #def __init__(self, INNER_LINE=0.35, OUTER_LINE = 0.35, INNER_R = 0.85, OUTER_R = 1.15):
        self.INNER_LINE = INNER_LINE #max offset for straightline (left)
        self.OUTER_LINE = OUTER_LINE #max offset for straightline (right)
        self.INNER_R = INNER_R #min R for circle (inner circle)
        self.OUTER_R = OUTER_R #max R for circle (outer circle)

def verify_reach(x_segs, y_segs, current_state, current_waypoint):
    '''
        Funtion: verify flowpipes approach waypoint
        Input: flowpipes x_segs, y_segs
               current_state
               current_waypoint as goal
        Output: diagnostic result, False safe, True unsafe       
    '''
    result = False
    d = list()  
    dis = np.sqrt((current_state[0]-current_waypoint[0])**2+(current_state[1]-current_waypoint[1])**2)
    d.append(dis)
    if dis <= 0.05:
        return False
    
    for i in range(len(x_segs)):
        x_mid = 0.5*(min(x_segs[i])+max(x_segs[i]))
        y_mid = 0.5*(min(y_segs[i])+max(y_segs[i]))
        d.append(np.sqrt((x_mid-current_waypoint[0])**2+(y_mid-current_waypoint[1])**2))
    for i in range(1, len(d)):
        if d[i] <= d[i-1]:
            continue
        else:
            result = True
            return result
    return result

def verify_dubin(veri_para, mode, x_list, y_list, current_state, current_waypoint):
    '''
        Funtion: verify intersection between one segment and corridor
        Input: veri_para: allowed offset from path
               mode: options: "straight", "circle", "rounded_square"
               current_state: [state.x, state.y, 0, state.yaw, state.vx, state.vy, state.yaw_dot]
               current_waypoint: [waypoint[0], waypoint[1], curvature]
        Output: diagnostic result, False safe, True unsafe       
    '''        
 
    INNER_LINE = veri_para.INNER_LINE
    OUTER_LINE = veri_para.OUTER_LINE
    INNER_R = veri_para.INNER_R
    OUTER_R = veri_para.OUTER_R
    dt = 0.1
    result = False 
    left = False
    right = False

    x = current_state[0]
    y = current_state[1]
    
    #key: current waypoint, value: center of the circle
    dict_waypont_center = {
        (2, 1, 1): (1, 1),
        (1, 3, 1): (1, 2),
        (-1, 2, 1): (0, 2),
        (0, 0, 1): (0, 1) 
    }
    #key: current waypoint, value: previous waypoint
    dict_pre_waypont = {
        (0, 0, 1): (-1, 1, 0),
        (1, 0, 0): (0, 0, 1),
        (2, 1, 1): (1, 0, 0),
        (2, 2, 0): (2, 1, 1),
        (1, 3, 1): (2, 2, 0),
        (0, 3, 0): (1, 3, 1),
        (-1, 2, 1): (0, 3, 0),
        (-1, 1, 0): (-1, 2, 1) 
    }

    if mode == "straight":
        sampling_x = np.arange(x, x+10, dt)
        sampling_y = 0*sampling_x
        sampling_y_l = sampling_y-INNER_LINE
        sampling_y_r = sampling_y+OUTER_LINE
        for i in range(len(x_list)):
            index = find_closet_index(x_list[i], sampling_x)
            if y_list[i] > sampling_y_l[index]:
                left = True
                break
            elif y_list[i] < sampling_y_r[index]:
                right = True
                break
            else:
                continue
        result = left | right
        return result
    elif mode == "circle":
        for i in range(len(x_list)):
            left =  ((x_list[i]-0)^2 + (x_list[i]-1)^2 < INNER_R^2)
            right = ((x_list[i]-0)^2 + (x_list[i]-1)^2 > OUTER_R^2)
            if left | right:
                break
            else:
                continue
        result = left | right
        return result
    elif mode == "rounded_square":
        if current_waypoint[2] == 0:#straght line
            if current_waypoint[0] == dict_pre_waypont[tuple(current_waypoint)][0]: #up or down
                for i in range(len(x_list)):
                    if x_list[i] > current_waypoint[0]+INNER_LINE or x_list[i] < current_waypoint[0]-INNER_LINE:
                        result = True
                        return result
            else: #left or right
                for i in range(len(y_list)):
                    if y_list[i] > current_waypoint[1]+INNER_LINE or y_list[i] < current_waypoint[1]-INNER_LINE:
                        result = True
                        return result
        else: #curve
            center = dict_waypont_center[tuple(current_waypoint)]
            for i in range(len(x_list)):
                if (x_list[i]-center[0])**2+(y_list[i]-center[1])**2>OUTER_R**2:
                    right = True
                    break
            for i in range(len(x_list)):
                if (x_list[i]-center[0])**2+(y_list[i]-center[1])**2<INNER_R**2:
                    left = True
                    break
            result = left | right
            return result        

def sampling_dubin_path(initial, target, r, dt):
    yaw = []
    if initial[0] == target[0]: #x doesn't change
        y = math_tool.range_by_step(initial[1], target[1]+dt, dt)
        x = [initial[0] for item in y]
    elif initial[1] == target[1]: #y doesn't change
        x = math_tool.range_by_step(initial[0], target[0]+dt, dt)
        y = [initial[1] for item in x]
    else:
        x = []
        y = []
        for item in math_tool.range_by_step(initial[2], target[2], dt):
            #print(item)
            xdot = r * math.cos(item)
            ydot = r * math.sin(item)
            if item == initial[2]:
                x.append(initial[0])
                y.append(initial[1])
            else:
                x.append(x[-1]+xdot*dt)
                y.append(y[-1]+ydot*dt)
    for i in range(len(x)):
        if i == 0:
            yaw.append(initial[2])
        else:
            dx = x[i]-x[i-1]
            dy = y[i]-y[i-1]
            yaw.append(math.atan2(dy, dx))
    return x, y, yaw
def get_closet_ref(x, y, yaw, current_waypoint):
    dict_waypont = {
        (1, 0, 0): (0, 0, 0, 0),
        (2, 1, 1): (1, 0, 0, np.pi/2),
        (2, 2, 0): (2, 1, np.pi/2, np.pi/2),
        (1, 3, 1): (2, 2, np.pi/2, np.pi),
        (0, 3, 0): (1, 3, np.pi, np.pi),
        (-1, 2, 1): (0, 3, -np.pi, -np.pi/2),
        (-1, 1, 0): (-1, 2, -np.pi/2, -np.pi/2),
        (0, 0, 1): (-1, 1, -np.pi/2, 0)}
    yaw0 =  dict_waypont[current_waypoint][2]
    yaw1 =  dict_waypont[current_waypoint][3]
    target = ((current_waypoint[0], current_waypoint[1], yaw1))
    initial = (dict_waypont[current_waypoint][0], dict_waypont[current_waypoint][1],yaw0)
    x_r, y_r, yaw_r = sampling_dubin_path(initial, target, 1, 0.1)
    d_min = 100
    ind = 0
    for i in range(len(x_r)):
        dis = np.sqrt((x-x_r[i])**2+(y-y_r[i])**2)
        if  dis <= d_min:
            d_min = dis
            ind = i
    return x_r[ind], y_r[ind], yaw_r[ind]
if __name__ == '__main__':

    dict_waypont = {
        (1, 0, 0): (0, 0, 0, 0),
        (2, 1, 1): (1, 0, 0, np.pi/2),
        (2, 2, 0): (2, 1, np.pi/2, np.pi/2),
        (1, 3, 1): (2, 2, np.pi/2, np.pi),
        (0, 3, 0): (1, 3, np.pi, np.pi),
        (-1, 2, 1): (0, 3, -np.pi, -np.pi/2),
        (-1, 1, 0): (-1, 2, -np.pi/2, -np.pi/2),
        (0, 0, 1): (-1, 1, -np.pi/2, 0) 
    }# key: current waypoint (end point); value previous waypoint x,y,previous yaw, current (end point) yaw
    current_waypoint = (0, 0, 1)
    yaw0 =  dict_waypont[current_waypoint][2]
    yaw1 =  dict_waypont[current_waypoint][3]
    target = ((current_waypoint[0], current_waypoint[1], yaw1))
    initial = (dict_waypont[current_waypoint][0], dict_waypont[current_waypoint][1],yaw0)
    x, y, yaw = sampling_dubin_path(initial, target, 1, 0.1)

    plt.plot(x, y, 'r')
    plt.show()
