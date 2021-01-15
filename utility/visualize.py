#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu
import numpy as np
import matplotlib.pyplot as plt
import math

def plot_reachflow_modelplex(fig, state, veri_para, x_list_i, y_list_i, vr_reach_seg, mode, verify_result):
    INNER_LINE = veri_para.INNER_LINE
    OUTER_LINE = veri_para.OUTER_LINE
    INNER_R = veri_para.INNER_R
    OUTER_R = veri_para.OUTER_R

    ax1 = fig.add_subplot(211)

    for i in xrange(len(x_list_i)):
        if vr_reach_seg[i]:
            ax1.plot(x_list_i[i], y_list_i[i], 'r')
        else:
            ax1.plot(x_list_i[i], y_list_i[i], 'g')
    if mode == "circle":
        plt.plot(state[0], state[1], 'r*', markersize=12)
        th = np.arange(0*np.pi, 2*np.pi+np.pi/10, np.pi/10)
        plt.plot(1 * np.cos(th), 1 * np.sin(th)+1, 'b', linewidth = 3.0)
        plt.xlim(-2, 2)
        plt.ylim(-1, 3)
    if mode == "rounded_square":
        ax1.plot(state[0], state[1], 'r*', markersize=12)#current position
        ax1.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1), 'b', linewidth = 3.0 )#1st
        ax1.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1)-INNER_LINE, 'r', linewidth = 3.0 )#1st left
        ax1.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1)+INNER_LINE, 'r', linewidth = 3.0 )#1st right

        th = np.arange(1.5*np.pi, 2*np.pi+np.pi/10, np.pi/10)
        ax1.plot(1 * np.cos(th)+1, 1 * np.sin(th)+1, 'b', linewidth = 3.0)#2nd
        ax1.plot(INNER_R * np.cos(th)+1, INNER_R * np.sin(th)+1, 'r', linewidth = 3.0)#2nd left
        ax1.plot(OUTER_R * np.cos(th)+1, OUTER_R * np.sin(th)+1, 'r', linewidth = 3.0)#2nd right

        ax1.plot(2+0*np.arange(1, 2.1, 0.1), np.arange(1, 2.1, 0.1), 'b', linewidth = 3.0)#3rd
        ax1.plot(2+0*np.arange(1, 2.1, 0.1)-INNER_LINE, np.arange(1, 2.1, 0.1), 'r', linewidth = 3.0)#3rd left
        ax1.plot(2+0*np.arange(1, 2.1, 0.1)+OUTER_LINE, np.arange(1, 2.1, 0.1), 'r', linewidth = 3.0)#3rd right

        th = np.arange(0*np.pi, 0.5*np.pi+np.pi/10, np.pi/10)
        ax1.plot(1 * np.cos(th)+1, 1 * np.sin(th)+2, 'b', linewidth = 3.0)#4th
        ax1.plot(INNER_R * np.cos(th)+1, INNER_R * np.sin(th)+2, 'r', linewidth = 3.0)#4th left
        ax1.plot(OUTER_R * np.cos(th)+1, OUTER_R * np.sin(th)+2, 'r', linewidth = 3.0)#4th right

        ax1.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1), 'b', linewidth = 3.0)#5th
        ax1.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1)-INNER_LINE, 'r', linewidth = 3.0)#5th left      
        ax1.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1)+OUTER_LINE, 'r', linewidth = 3.0)#5th 

        th = np.arange(0.5*np.pi, 1.0*np.pi+np.pi/10, np.pi/10)
        ax1.plot(1 * np.cos(th)+0, 1 * np.sin(th)+2, 'b', linewidth = 3.0)#6th
        ax1.plot(INNER_R * np.cos(th)+0, INNER_R * np.sin(th)+2, 'r', linewidth = 3.0)#6th left
        ax1.plot(OUTER_R * np.cos(th)+0, OUTER_R * np.sin(th)+2, 'r', linewidth = 3.0)#6th right

        ax1.plot(-1+0*np.arange(0, 1.1, 0.1), np.arange(1, 2.1, 0.1), 'b', linewidth = 3.0)#7th
        ax1.plot(-1+0*np.arange(0, 1.1, 0.1)+INNER_LINE, np.arange(1, 2.1, 0.1), 'r', linewidth = 3.0)#7th left
        ax1.plot(-1+0*np.arange(0, 1.1, 0.1)-OUTER_LINE, np.arange(1, 2.1, 0.1), 'r', linewidth = 3.0)#7th right

        th = np.arange(np.pi, 1.5*np.pi+np.pi/10, np.pi/10)
        ax1.plot(1 * np.cos(th)+0, 1 * np.sin(th)+1, 'b', linewidth = 3.0)#8th
        ax1.plot(INNER_R * np.cos(th)+0, INNER_R * np.sin(th)+1, 'r', linewidth = 3.0)#8th left
        ax1.plot(OUTER_R * np.cos(th)+0, OUTER_R * np.sin(th)+1, 'r', linewidth = 3.0)#8th right
        ax1.set_xlim(-2, 3)
        ax1.set_ylim(-1, 4)

        ax2 = fig.add_subplot(212)
        if (len(verify_result[0]) > 5) & (len(verify_result[1]) > 5):
            ax2.plot(verify_result[0][-5:], 'b-*', label = 'Reach_flow')
            ax2.plot(verify_result[1][-5:], 'r', label = 'Keymera Monitor')
            ax2.legend()
            ax2.set_xlim(0, 4)
            ax2.set_ylim(-0.5, 1.5)
        fig.canvas.draw()
        plt.clf()

def plot_reachflow(ax, state, veri_para, x_list_i, y_list_i, vr_reach_seg, mode):
    
    INNER_LINE = veri_para.INNER_LINE
    OUTER_LINE = veri_para.OUTER_LINE
    INNER_R = veri_para.INNER_R
    OUTER_R = veri_para.OUTER_R
    line_width_reach = 1

    #ax1 = fig.add_subplot(211)
    if len(x_list_i) == 0:
        print('did not get any reachable segs')
        return
    if len(x_list_i) != len(vr_reach_seg):
        print('reachable segs size error')
        return

    for i in range(len(x_list_i)):
        if vr_reach_seg[i]:
            ax.plot(x_list_i[i]+[x_list_i[i][0]], y_list_i[i]+[y_list_i[i][0]], 'r', line_width_reach)
        else:
            ax.plot(x_list_i[i]+[x_list_i[i][0]], y_list_i[i]+[y_list_i[i][0]], 'green', line_width_reach)
    if mode == "circle":
        plt.plot(state[0], state[1], 'r*', markersize=12)
        th = np.arange(0*np.pi, 2*np.pi+np.pi/10, np.pi/10)
        plt.plot(1 * np.cos(th), 1 * np.sin(th)+1, 'b', line_width_reach)
        plt.xlim(-2, 2)
        plt.ylim(-1, 3)
    if mode == "rounded_square":
        ax.plot(state[0], state[1], 'b*', markersize=12)#current position
        #ax.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1), 'b', linewidth = 3.0 )#1st
        ax.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1)-INNER_LINE, 'r', line_width_reach)#1st left
        ax.plot(np.arange(0, 1.1, 0.1), 0*np.arange(0, 1.1, 0.1)+INNER_LINE, 'r', line_width_reach)#1st right

        th = np.arange(1.5*np.pi, 2*np.pi+np.pi/10, np.pi/10)
        #ax.plot(1 * np.cos(th)+1, 1 * np.sin(th)+1, 'b', linewidth = 3.0)#2nd
        ax.plot(INNER_R * np.cos(th)+1, INNER_R * np.sin(th)+1, 'r', line_width_reach)#2nd left
        ax.plot(OUTER_R * np.cos(th)+1, OUTER_R * np.sin(th)+1, 'r', line_width_reach)#2nd right

        #ax.plot(2+0*np.arange(1, 2.1, 0.1), np.arange(1, 2.1, 0.1), 'b', linewidth = 3.0)#3rd
        ax.plot(2+0*np.arange(1, 2.1, 0.1)-INNER_LINE, np.arange(1, 2.1, 0.1), 'r', line_width_reach)#3rd left
        ax.plot(2+0*np.arange(1, 2.1, 0.1)+OUTER_LINE, np.arange(1, 2.1, 0.1), 'r', line_width_reach)#3rd right

        th = np.arange(0*np.pi, 0.5*np.pi+np.pi/10, np.pi/10)
        #ax.plot(1 * np.cos(th)+1, 1 * np.sin(th)+2, 'b', linewidth = 3.0)#4th
        ax.plot(INNER_R * np.cos(th)+1, INNER_R * np.sin(th)+2, 'r', line_width_reach)#4th left
        ax.plot(OUTER_R * np.cos(th)+1, OUTER_R * np.sin(th)+2, 'r', line_width_reach)#4th right

        #ax.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1), 'b', linewidth = 3.0)#5th
        ax.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1)-INNER_LINE, 'r', line_width_reach)#5th left      
        ax.plot(np.arange(0, 1.1, 0.1), 3+0*np.arange(0, 1.1, 0.1)+OUTER_LINE, 'r', line_width_reach)#5th 

        th = np.arange(0.5*np.pi, 1.0*np.pi+np.pi/10, np.pi/10)
        #ax.plot(1 * np.cos(th)+0, 1 * np.sin(th)+2, 'b', linewidth = 3.0)#6th
        ax.plot(INNER_R * np.cos(th)+0, INNER_R * np.sin(th)+2, 'r', line_width_reach)#6th left
        ax.plot(OUTER_R * np.cos(th)+0, OUTER_R * np.sin(th)+2, 'r', line_width_reach)#6th right

        #ax.plot(-1+0*np.arange(0, 1.1, 0.1), np.arange(1, 2.1, 0.1), 'b', linewidth = 3.0)#7th
        ax.plot(-1+0*np.arange(0, 1.1, 0.1)+INNER_LINE, np.arange(1, 2.1, 0.1), 'r', line_width_reach)#7th left
        ax.plot(-1+0*np.arange(0, 1.1, 0.1)-OUTER_LINE, np.arange(1, 2.1, 0.1), 'r', line_width_reach)#7th right

        th = np.arange(np.pi, 1.5*np.pi+np.pi/10, np.pi/10)
        #ax.plot(1 * np.cos(th)+0, 1 * np.sin(th)+1, 'b', linewidth = 3.0)#8th
        ax.plot(INNER_R * np.cos(th)+0, INNER_R * np.sin(th)+1, 'r', line_width_reach)#8th left
        ax.plot(OUTER_R * np.cos(th)+0, OUTER_R * np.sin(th)+1, 'r', line_width_reach)#8th right
        ax.set_xlim(-1.5, 2.5)
        ax.set_ylim(-1, 4)
def plot_car(ax, x, y, yaw, steer=0.0, cabcolor="-r", truckcolor="-k", linewidth = 1.5):
    # Vehicle parameters for visualization only
    LENGTH = 0.257  # [m]
    WIDTH = 0.257/2  # [m]
    BACKTOWHEEL = 0.05  # [m]
    WHEEL_LEN = 0.05  # [m]
    WHEEL_WIDTH = 0.02  # [m]
    TREAD = 0.05  # [m]
    WB = 0.257  # [m]

    outline = np.array([[-BACKTOWHEEL, (LENGTH - BACKTOWHEEL), (LENGTH - BACKTOWHEEL), -BACKTOWHEEL, -BACKTOWHEEL],
                        [WIDTH / 2, WIDTH / 2, - WIDTH / 2, -WIDTH / 2, WIDTH / 2]])

    fr_wheel = np.array([[WHEEL_LEN, -WHEEL_LEN, -WHEEL_LEN, WHEEL_LEN, WHEEL_LEN],
                         [-WHEEL_WIDTH - TREAD, -WHEEL_WIDTH - TREAD, WHEEL_WIDTH - TREAD, WHEEL_WIDTH - TREAD, -WHEEL_WIDTH - TREAD]])

    rr_wheel = np.copy(fr_wheel)

    fl_wheel = np.copy(fr_wheel)
    fl_wheel[1, :] *= -1
    rl_wheel = np.copy(rr_wheel)
    rl_wheel[1, :] *= -1

    Rot1 = np.array([[math.cos(yaw), math.sin(yaw)],
                     [-math.sin(yaw), math.cos(yaw)]])
    Rot2 = np.array([[math.cos(steer), math.sin(steer)],
                     [-math.sin(steer), math.cos(steer)]])

    fr_wheel = (fr_wheel.T.dot(Rot2)).T
    fl_wheel = (fl_wheel.T.dot(Rot2)).T
    fr_wheel[0, :] += WB
    fl_wheel[0, :] += WB

    fr_wheel = (fr_wheel.T.dot(Rot1)).T
    fl_wheel = (fl_wheel.T.dot(Rot1)).T

    outline = (outline.T.dot(Rot1)).T
    rr_wheel = (rr_wheel.T.dot(Rot1)).T
    rl_wheel = (rl_wheel.T.dot(Rot1)).T

    outline[0, :] += x
    outline[1, :] += y
    fr_wheel[0, :] += x
    fr_wheel[1, :] += y
    rr_wheel[0, :] += x
    rr_wheel[1, :] += y
    fl_wheel[0, :] += x
    fl_wheel[1, :] += y
    rl_wheel[0, :] += x
    rl_wheel[1, :] += y


    ax.plot(np.array(outline[0, :]).flatten(),
             np.array(outline[1, :]).flatten(), truckcolor, linewidth=linewidth)
    ax.plot(np.array(fr_wheel[0, :]).flatten(),
             np.array(fr_wheel[1, :]).flatten(), truckcolor, linewidth=linewidth)
    ax.plot(np.array(rr_wheel[0, :]).flatten(),
             np.array(rr_wheel[1, :]).flatten(), truckcolor, linewidth=linewidth)
    ax.plot(np.array(fl_wheel[0, :]).flatten(),
             np.array(fl_wheel[1, :]).flatten(), truckcolor, linewidth=linewidth)
    ax.plot(np.array(rl_wheel[0, :]).flatten(),
             np.array(rl_wheel[1, :]).flatten(), truckcolor, linewidth=linewidth)
    ax.plot(x, y, "*")
