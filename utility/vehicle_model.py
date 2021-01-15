#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu

#RC-car
#L=0.257, L_f=0.137, L_r=0.12, C_alpha=56.4, I_z=0.0558
#matlab_simulator: 
# vehicle.L = 0.257;          % wheelbase (m)
# vehicle.L_f = 0.115;        % CoG to front axle (m)
# vehicle.L_r = 0.142;        % CoG to rear axle (m)
# vehicle.tw = 0.165;         % trackwidth (m)

# vehicle.wheel_dia = 0.0323; % wheel diameter (m)
# vehicle.wheel_w = 0.025;    % wheel width (m)

# vehicle.m = 2.596;          % mass (kg)
# g = 9.81;                   % gravity const (m/s^2)

# vehicle.load_f = vehicle.m*g*vehicle.L_r/vehicle.L;
# vehicle.load_r = vehicle.m*g*vehicle.L_f/vehicle.L;

# vehicle.C_x = 103.94;       % longitudinal stiffness (N)
# vehicle.C_alpha = 56.4;     % cornering stiffness (N)
# vehicle.I_z = 0.0558;       % rotation inertia (kgm^2)
# vehicle.mu = 1.37;          % friction coefficient
# vehicle.mu_slide = 1.96;    % sliding friction coefficient

#MRZR
#m: 879.0
#L_f = 1.364
#L_r = 1.364
#load_f = 4307.1
#load_r = 4307.1
#C_x = 13782
#C_alpha = 68912
#I_z = 1020.0
#mu_k = 1.37
#mu_s = 1.96
#tw = 1.295
#wheel_dia = 0.66
#wheel_w = 0.229
#r_min = 5.0
#max_wp_dist = 15.0
#max_v = 6.0
def load_model(vehicle_type):
    if vehicle_type == "RC_Car":
        class Model:
            def __init__(self, m=2.645, L=0.257, width = 0.129, L_f=0.137, L_r=0.12, load_f=12.1, load_r=13.8178, C_x=103.94, C_alpha=56.4, I_z=0.0558, mu=1.37, mu_s=1.96):
            #def __init__(self, m=2.672, L=0.257, L_f=0.1343, L_r=0.12, load_f=12.1, load_r=13.8178, C_x=103.94, C_alpha=56.4, I_z=0.0558, mu=1.37, mu_s=1.96):
                self.m = m
                self.L = L
                self.width = width
                self.L_f = L_f
                self.L_r = L_r
                self.load_f = load_f
                self.load_r = load_r
                self.C_x = C_x
                self.C_alpha = C_alpha
                self.I_z = I_z
                self.mu = mu
                self.mu_s = mu_s
    if vehicle_type == "MRZR":
        class Model:
            def __init__(self, m=879, L=2.728, L_f=1.364, L_r=1.364, load_f=4307.1, load_r=4307.1, C_x=13786, C_alpha=68912, I_z=1020, mu_k=1.37, mu_s=1.96, tw=1.295, wheel_dia=0.66, wheel_w=0.229, r_min=5.0, max_wp_dist=15.0, max_v=6.0):
                self.m = m
                self.L = L
                self.L_f = L_f
                self.L_r = L_r
                self.load_f = load_f
                self.load_r = load_r
                self.C_x = C_x
                self.C_alpha = C_alpha
                self.I_z = I_z
                self.mu_k = mu_k
                self.mu_s = mu_s
                self.tw = tw
                self.wheel_dia = wheel_dia
                self.wheel_w = wheel_w
                self.r_min = r_min
                self.max_wp_dist = max_wp_dist
                self.max_v = max_v
    return Model