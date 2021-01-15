#!/usr/bin/env python

import math
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from polytope_2_7 import Polytope
import scipy.linalg as sp_linalg

np.set_printoptions(threshold=np.inf, precision=5, suppress=True)

T, N = 1, 1 # horizon length for MPC (both T and N should be same)

NX = 4  # x = [x, y, v, yaw]
NU = 2  # u = [accel, steer]

R = np.diag([1, 1])  # input cost matrix

Q = np.diag([100.0, 100.0, 1.0, 100.0]) # state cost matrix

# Disturbance

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

Uc = Polytope(np.array([[1.25, 6.07],
    					[1.25, -6.07],
    					[-1.25, -6.07],
    					[-1.25, 6.07]]))

DT = 0.1

# Vehicle parameters
WB = 0.257 #, L = 2.89, 2.89  # Wheel base [m]
MAX_STEER = np.deg2rad(80.0)  # maximum steering angle [rad]
MAX_ACCEL = 2  # maximum accel [m/ss]

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

# GUROBI variables
m = gp.Model("MPC")
m.setParam('OutputFlag', False)
x = np.array([[m.addVar(name="x"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(T + 1)] for j in range(NX)])
u = np.array([[m.addVar(name="u"+str(i)+"0", vtype=GRB.CONTINUOUS, lb=-MAX_ACCEL, ub=MAX_ACCEL) for i in range(T)],
              [m.addVar(name="u"+str(i)+"1", vtype=GRB.CONTINUOUS, lb=-MAX_STEER, ub=MAX_STEER) for i in range(T)]])
xref_var = np.array([[m.addVar(name="ref"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(T + 1)] for j in range(NX)])
a = np.array([[m.addVar(name="a"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(NX)] for j in range(NX)])
b = np.array([[m.addVar(name="b"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(NU)] for j in range(NX)])
c = np.array([m.addVar(name="c"+str(i), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(NX)])
x0_var = np.array([m.addVar(name="n0"+str(i), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(2)])
result1_var = np.array([[m.addVar(name="result1"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(3)] for j in range(10)])
result2_var = np.array([[m.addVar(name="result2"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(3)] for j in range(4)])
X_var = np.array([[m.addVar(name="X"+str(i)+str(j), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=GRB.INFINITY) for i in range(3)] for j in range(10)])

cost = gp.QuadExpr()
cost += (xref_var[:, 0] - x[:, 0]).dot((Q).dot(xref_var[:, 0] - x[:, 0]))
for t in range(T):
    cost += u[:, t].dot(R.dot(u[:, t]))
    cost += (xref_var[:, t + 1] - x[:, t + 1]).dot((Q).dot(xref_var[:, t + 1] - x[:, t + 1]))

    for n in range(NX):
        m.addConstr(x[n, t + 1] == a.dot(x[:, t])[n] + b.dot(u[:, t])[n] + c[n])

    # if result1 is not None:
    for result in result1_var:
        m.addConstr(result[0] * x[0, t + 1] + result[1] * x[1, t + 1] <= result[2])
    for result in result2_var:
        m.addConstr(result[0] * u[0, t] + result[1] * u[1, t] <= result[2])

for n in range(2, NX):
    m.addConstr(x[n,0] == x0_var[n - 2])

for result in X_var:
    m.addConstr(result[0] * x[0, 0] + result[1] * x[1, 0] <= result[2])
m.update()
m.setObjective(cost, GRB.MINIMIZE)


# Set computations and GUROBI optimization for MPC
def linear_mpc_control(xref, ubar, x0, obstacles):
    oa = ubar[0,0:N]   
    odelta = ubar[1,0:N]  

    A, B, C = get_linear_model_matrix(x0[2], x0[3], ubar[1,0]) 

    result1, result2, X = None, None, None
    XX = sp_linalg.solve_discrete_are(A, B, Q, R)
 

    K = np.dot(np.linalg.pinv(R + np.dot(B.T, np.dot(XX, B))), np.dot(B.T, np.dot(XX, A)))

    # X and U bounds polytope calculations
    Ak = np.array(A - B.dot(K))

    orig_Z = W + W.mult(Ak)

    Z = Polytope(orig_Z.V[:,0:2])
    Xc_robust = Polytope(obstacles)
    Uc_robust = Uc
    try:
        
        Xc_robust = Xc_robust - Z
        Uc_robust = Uc_robust - orig_Z.mult(K)
        print("robust set obtained")
    except:
    	return oa.T, odelta.T

    result1 = np.hstack((Xc_robust.A, Xc_robust.b))
    result1 = np.vstack((result1, np.zeros((10 - result1.shape[0], 3))))
    result2 = np.hstack((Uc_robust.A, Uc_robust.b))
    
    x_init = Z.add(x0[0:2])
        
    X = np.hstack((x_init.A, x_init.b))
    X = np.vstack((X, np.zeros((10 - X.shape[0], 3))))

    for n in range(NX):
        for t in range(T + 1):
            xref_var[n,t].setAttr(GRB.Attr.LB, max(min(xref[n,t], 1e20), -1e20))
            xref_var[n,t].setAttr(GRB.Attr.UB, max(min(xref[n,t], 1e20), -1e20))

        c[n].setAttr(GRB.Attr.LB, max(min(C[n], 1e20), -1e20))
        c[n].setAttr(GRB.Attr.UB, max(min(C[n], 1e20), -1e20))
        for n2 in range(NX):
            a[n, n2].setAttr(GRB.Attr.LB, max(min(A[n, n2], 1e20), -1e20))
            a[n, n2].setAttr(GRB.Attr.UB, max(min(A[n, n2], 1e20), -1e20))
        for n2 in range(NU):
            b[n, n2].setAttr(GRB.Attr.LB, max(min(B[n, n2], 1e20), -1e20))
            b[n, n2].setAttr(GRB.Attr.UB, max(min(B[n, n2], 1e20), -1e20))

    for i in range(10):
        for j in range(3):
            result1_var[i,j].setAttr(GRB.Attr.LB, max(min(result1[i,j], 1e20), -1e20))
            result1_var[i,j].setAttr(GRB.Attr.UB, max(min(result1[i,j], 1e20), -1e20))
            X_var[i, j].setAttr(GRB.Attr.LB, max(min(X[i, j], 1e20), -1e20))
            X_var[i, j].setAttr(GRB.Attr.UB, max(min(X[i, j], 1e20), -1e20))

    for i in range(4):
        for j in range(3):
            result2_var[i,j].setAttr(GRB.Attr.LB, max(min(result2[i,j], 1e20), -1e20))
            result2_var[i,j].setAttr(GRB.Attr.UB, max(min(result2[i,j], 1e20), -1e20))

    for n in range(2, 4):
        x0_var[n - 2].setAttr(GRB.Attr.LB, max(min(x0[n], 1e20), -1e20))
        x0_var[n - 2].setAttr(GRB.Attr.UB, max(min(x0[n], 1e20), -1e20))

    m.update()
    m.optimize()  
    oa = ubar[0,0:N]
    odelta = ubar[1,0:N]
    if m.status == GRB.OPTIMAL:
        for i, uu in enumerate(u[0, :]):
            oa[i] = uu.X
            if (i == 0):
                oa[i] -= K[0,0] * (x0[0] - x[0,0].X) + K[0,1] * (x0[1] - x[1,0].X)
        for i, uu in enumerate(u[1, :]):
            odelta[i] = uu.X
            if (i == 0):
                odelta[i] -= K[1, 0] * (x0[0] - x[0,0].X) + K[1,1] * (x0[1] - x[1,0].X)
    else:
        print("Error: Cannot solve mpc..")
    return oa.T, odelta.T


def mpc(X0, X_ref):
    U = np.zeros((N, 2, 1))  # warm start for U(control). Has no significance with frenet path planning. Was used with CILQR planner to get an initial U from CILQR.
    region_V = np.array(
                  [
                    [-10, -10],
                    [-10, 10],
                    [10, 10],
                    [10, -10]
                  ]
                 )

    U[0:N,0,0], U[0:N,1,0] = linear_mpc_control(X_ref.T[0], U.T[0], X0.T[0], region_V)
    return U[0][0][0], U[0][1][0]
