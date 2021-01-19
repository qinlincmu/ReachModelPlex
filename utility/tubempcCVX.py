#!/usr/bin/env python

import math
import cvxpy
from cvxpy import *
import numpy as np
#import gurobipy as gp
#from gurobipy import GRB
from polytope_2_7 import Polytope
import scipy.linalg as sp_linalg

np.set_printoptions(threshold=np.inf, precision=5, suppress=True)

T = 1
MAX_ITER = 1
DU_TH = 0.1  # iteration finish param

NX = 4  # x = [x, y, v, yaw]
NU = 2  # u = [accel, steer]

# R = np.diag([1, 1])  # input cost matrix

# Q = np.diag([100.0, 100.0, 1.0, 100.0]) # state cost matrix

R = np.diag([1, 1])#np.diag([0.01, 0.01])  # input cost matrix
Rd = np.diag([0.01, 0.01])  # input difference cost matrix
Q = np.diag([100.0, 100.0, 0.5, 100])#np.diag([100.0, 100.0, 1.0, 100.0])#np.diag([1.0, 1.0, 0.5, 0.5])  # state cost matrix
Qf = Q  # state final matrix

MAX_STEER = np.deg2rad(90.0)  # maximum steering angle [rad]
MAX_DSTEER = np.deg2rad(30.0)  # maximum steering speed [rad/s]
MAX_SPEED = 2#55.0 / 3.6  # maximum speed [m/s]
MIN_SPEED = 0.1#-20.0 / 3.6  # minimum speed [m/s]
MAX_ACCEL = 2.0  # maximum accel [m/ss]

# Disturbance

class State:
    """
    vehicle state class
    """

    def __init__(self, x=0.0, y=0.0, yaw=0.0, v=0.0):
        self.x = x
        self.y = y
        self.yaw = yaw
        self.v = v
        self.predelta = None

W_vertex =  0.005*np.array([[1, 1, 1, 1],
               		  [1, -1, 1, 1],
               		  [-1, -1, 1, 1],
               		  [-1, 1, 1, 1],
               		  [1, 1, 1, -1],
               		  [1, -1, 1, -1],
               		  [-1, -1, 1, -1],
               		  [-1, 1, 1, -1],
               		  [1, 1, -1, 1],
               		  [1, -1, -1, 1],
               		  [-1, -1, -1, 1],
               		  [-1, 1, -1, 1],
               		  [1, 1, -1, -1],
               		  [1, -1, -1, -1],
               		  [-1, -1, -1, -1],
               		  [-1, 1, -1, -1]], dtype=np.float)

W = Polytope(W_vertex)

Uc = Polytope(np.array([[2, 6.07],
    					[2, -6.07],
    					[-2, -6.07],
    					[-2, 6.07]]))



DT = 0.1

# Vehicle parameters
WB = 0.257 # Wheel base [m]

def get_linear_model_matrix(v, phi, delta):
    # v : absolute velocity
    # phi : orientation in radians
    # delta : steering (0)
    A = np.zeros((NX, NX))
    A[0, 0] = 1.0
    A[1, 1] = 1.0
    A[2, 2] = 1.0
    A[3, 3] = 1.0
    A[0, 2] = DT * math.cos(phi)
    A[0, 3] = - DT * v * math.sin(phi)
    A[1, 2] = DT * math.sin(phi)
    A[1, 3] = DT * v * math.cos(phi)
    A[3, 2] = DT * math.tan(delta) / WB

    B = np.zeros((NX, NU))
    B[2, 0] = DT
    B[3, 1] = DT * v / (WB * math.cos(delta) ** 2)

    C = np.zeros(NX)
    C[0] = DT * v * math.sin(phi) * phi
    C[1] = - DT * v * math.cos(phi) * phi
    C[3] = - DT * v * delta / (WB * math.cos(delta) ** 2)

    return A, B, C
def update_state(state, a, delta):

    # input check
    if delta >= MAX_STEER:
        delta = MAX_STEER
    elif delta <= -MAX_STEER:
        delta = -MAX_STEER

    state.x = state.x + state.v * math.cos(state.yaw) * DT
    state.y = state.y + state.v * math.sin(state.yaw) * DT
    state.yaw = state.yaw + state.v / WB * math.tan(delta) * DT
    state.v = state.v + a * DT

    if state. v > MAX_SPEED:
        state.v = MAX_SPEED
    elif state. v < MIN_SPEED:
        state.v = MIN_SPEED

    return state
def get_nparray_from_matrix(x):
    return np.array(x).flatten()

def predict_motion(x0, oa, od, xref):
    xbar = xref * 0.0
    for i, _ in enumerate(x0):
        xbar[i, 0] = x0[i]

    state = State(x=x0[0], y=x0[1], yaw=x0[3], v=x0[2])
    for (ai, di, i) in zip(oa, od, range(1, T + 1)):
        state = update_state(state, ai, di)
        xbar[0, i] = state.x
        xbar[1, i] = state.y
        xbar[2, i] = state.v
        xbar[3, i] = state.yaw

    return xbar

def linear_mpc_control_robust(xref, xbar, x0, dref):
    """
    linear robust mpc control

    xref: reference point
    xbar: operational point
    x0: initial state
    dref: reference steer angle
    """
    x = cvxpy.Variable((NX, T + 1))
    u = cvxpy.Variable((NU, T))

    cost = 0.0
    constraints = []
    region_V = np.array(
                  [
                    [-10, -10],
                    [-10, 10],
                    [10, 10],
                    [10, -10]
                  ]
                 ) #feasible zone, TODO, this is useful for presence of obstacle, we can use iris to get convex-free area
    Xc_robust = Polytope(region_V)

    for t in range(T):
        cost += cvxpy.quad_form(u[:, t], R)

        if t != 0:
            cost += cvxpy.quad_form(xref[:, t] - x[:, t], Q)

        A, B, C = get_linear_model_matrix(
            xbar[2, t], xbar[3, t], dref[0, t])
        constraints += [x[:, t + 1] == A * x[:, t] + B * u[:, t] + C]

        XX = sp_linalg.solve_discrete_are(A, B, Q, R)
 

        K = np.dot(np.linalg.pinv(R + np.dot(B.T, np.dot(XX, B))), np.dot(B.T, np.dot(XX, A)))

        # X and U bounds polytope calculations
        Ak = np.array(A - B.dot(K))

        orig_Z = W + W.mult(Ak)

        Z = Polytope(orig_Z.V[:,0:2])
        
        Uc_robust = Uc

        try:
            Xc_robust = Xc_robust - Z
            Uc_robust = Uc_robust - orig_Z.mult(K)
            print("robust set obtained")
            print("xc robust: ", Xc_robust.A, Xc_robust.b, Xc_robust.A.shape, Xc_robust.b.shape)
            print("uc robust: ", Uc_robust.A, Uc_robust.b, Uc_robust.A.shape, Uc_robust.b.shape)
        except:
            print("can not get robust set")
            return None, None, None, None, None, None

        if t < (T - 1):
            cost += cvxpy.quad_form(u[:, t + 1] - u[:, t], Rd)
            constraints += [cvxpy.abs(u[1, t + 1] - u[1, t]) <=
                            MAX_DSTEER * DT]
        constraints += [abs(x[:, t]) <= max(Xc_robust.b)]
        constraints += [abs(u[0, t]) <= max(Xc_robust.b)]
        constraints += [abs(u[0, t]) <= min(Xc_robust.b)]
        # constraints += [Uc_robust.A*u[:, t] <= Uc_robust.b]
    cost += cvxpy.quad_form(xref[:, T] - x[:, T], Qf)

    constraints += [x[:, 0] == x0]
    constraints += [x[2, :] <= MAX_SPEED]
    constraints += [x[2, :] >= MIN_SPEED]
    constraints += [cvxpy.abs(u[0, :]) <= MAX_ACCEL]
    constraints += [cvxpy.abs(u[1, :]) <= MAX_STEER]

    prob = cvxpy.Problem(cvxpy.Minimize(cost), constraints)
    prob.solve(solver=cvxpy.ECOS, verbose=False)

    if prob.status == cvxpy.OPTIMAL or prob.status == cvxpy.OPTIMAL_INACCURATE:
        ox = get_nparray_from_matrix(x.value[0, :])
        oy = get_nparray_from_matrix(x.value[1, :])
        ov = get_nparray_from_matrix(x.value[2, :])
        oyaw = get_nparray_from_matrix(x.value[3, :])
        oa = get_nparray_from_matrix(u.value[0, :])
        odelta = get_nparray_from_matrix(u.value[1, :])

    else:
        print("Error: Cannot solve mpc..")
        oa, odelta, ox, oy, oyaw, ov = None, None, None, None, None, None

    return oa, odelta, ox, oy, oyaw, ov


    result1, result2, X = None, None, None



def linear_mpc_control(xref, xbar, x0, dref):
    """
    normal linear mpc control

    xref: reference point
    xbar: operational point
    x0: initial state
    dref: reference steer angle
    """

    x = cvxpy.Variable((NX, T + 1))
    u = cvxpy.Variable((NU, T))

    cost = 0.0
    constraints = []

    for t in range(T):
        cost += cvxpy.quad_form(u[:, t], R)

        if t != 0:
            cost += cvxpy.quad_form(xref[:, t] - x[:, t], Q)

        A, B, C = get_linear_model_matrix(
            xbar[2, t], xbar[3, t], dref[0, t])
        constraints += [x[:, t + 1] == A * x[:, t] + B * u[:, t] + C]

        if t < (T - 1):
            cost += cvxpy.quad_form(u[:, t + 1] - u[:, t], Rd)
            constraints += [cvxpy.abs(u[1, t + 1] - u[1, t]) <=
                            MAX_DSTEER * DT]

    cost += cvxpy.quad_form(xref[:, T] - x[:, T], Qf)

    constraints += [x[:, 0] == x0]
    constraints += [x[2, :] <= MAX_SPEED]
    constraints += [x[2, :] >= MIN_SPEED]
    constraints += [cvxpy.abs(u[0, :]) <= MAX_ACCEL]
    constraints += [cvxpy.abs(u[1, :]) <= MAX_STEER]

    prob = cvxpy.Problem(cvxpy.Minimize(cost), constraints)
    prob.solve(solver=cvxpy.ECOS, verbose=False)

    if prob.status == cvxpy.OPTIMAL or prob.status == cvxpy.OPTIMAL_INACCURATE:
        ox = get_nparray_from_matrix(x.value[0, :])
        oy = get_nparray_from_matrix(x.value[1, :])
        ov = get_nparray_from_matrix(x.value[2, :])
        oyaw = get_nparray_from_matrix(x.value[3, :])
        oa = get_nparray_from_matrix(u.value[0, :])
        odelta = get_nparray_from_matrix(u.value[1, :])

    else:
        print("Error: Cannot solve mpc..")
        oa, odelta, ox, oy, oyaw, ov = None, None, None, None, None, None

    return oa, odelta, ox, oy, oyaw, ov
def mpc(xref, x0, dref, oa, odelta):
    """
    MPC contorl with updating operational point iteraitvely
    """

    if oa is None or od is None:
        oa = [0.0] * T
        od = [0.0] * T

    xbar = predict_motion(x0, oa, od, xref)
    #oa, od, ox, oy, oyaw, ov = linear_mpc_control(xref, xbar, x0, dref)
    oa, od, ox, oy, oyaw, ov = linear_mpc_control_robust(xref, xbar, x0, dref)
    # for i in range(MAX_ITER): #iterative version of mpc, due to real-time requirement, we don't use it
    #     xbar = predict_motion(x0, oa, od, xref)
    #     poa, pod = oa[:], od[:]
    #     oa, od, ox, oy, oyaw, ov = linear_mpc_control(xref, xbar, x0, dref)
    #     print("dref od ..............", dref.shape, od.shape) #1,6 5
    #     #oa, od, ox, oy, oyaw, ov = linear_mpc_control(xref, xbar, x0, pod) #modified by qin 
    #     du = sum(abs(oa - poa)) + sum(abs(od - pod))  # calc u change value
    #     print(du)
    #     # if du <= DU_TH:
    #     #     break
    # else:
    #     print("Iterative is max iter")

    return oa, od, ox, oy, oyaw, ov