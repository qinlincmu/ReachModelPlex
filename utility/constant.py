#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu

T_CYCLE_TIME = 1.0     #// [ds]
MAX_MOTOR_ACCEL = 25.0 #// [dm/s^2], see vesc_interface.yaml
DT = 0.1 #simulation sanpling time: 0.02 for dynamic model, 0.1 for kinematic model
ODE = "nonlinear_without_beta" #"nonlinear_without_beta", "linear_tire"
flowstar_stepsize = 0.1 #0.001 for dynamic model, 0.1 for kinematic model
horizon = int(DT/flowstar_stepsize) #steps for reachability analysis, current control
horizon_b = 1*horizon #2 for dynamic experiment, 1 for kinematic experiment. 2 means predicting 2-step control cycles for full braking, 1 means 1-step
inject_bad_control = True #switch for on/off disturbance
fallback = True #switch for on/off fallback control