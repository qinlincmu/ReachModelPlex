# ReachModelPlex
Project of integrating reachability and ModelPlex Monitor

## All python dependencies

    matplotlib==3.0.3

    numpy==1.18.2

    scipy==1.4.1

    cffi==1.14.2

    mpmath==1.1.0

    pyinterval==1.2.1.dev0

    pycddlib==2.1.2

    cvxpy==1.1.5

## set up reachability tool
    - put "flowstar" file and model file (nonlinear_without_beta) file in the same level of directory
    - install gcc-8 and g++-8, so far flow* only works with at least version-8 compilers
    - install latest m4, gmp, mpfr, gsl, glpk, bison, and flex. (one possible issue in Ubuntu 18, please install flex 2.6.3 instaed of 2.6.4)
    - type "make" in the flowstar file, then type "sudo ldconfig"
    - type "make" in the model file, if everything goes well, you'll find an executable RC_bicycle program
    NOTE: so far the prototype is only working in Linux, only Ubuntu16.04 LTS has been tested.

## To run the simulation, we just execute the main python file: simulator.py

## Introduction of all subfolders:
    aa_planner: store policy of CPO learning-based controller
    flowstar: library of flowstar - a reachability tool
    nonlinear_without_beta: kinematic bicycle model of the vehicle
    utility: subfunctions used
            constant.py: parameters used for ModelPlex monitor
            flowstar.py: execution of flowstar in python
            interval_estimation: variable interval estimation using interval arithmetic
            libmonitor.so: shared library compiled from monitor's condition (C code)
            math_tool: some simple math utility functions
            modelplex: conversion functions to transfer robot's state into modelplex's state for further use in safety checking using modelplex
            polytope_2_7: external function used for representing disturbance for tube mpc
            simulate.py: model-based state prediction
            tubempc.py: optimal control using tubempc
            vehicle_model.py: vehicle model's parameters
            verification: verify reachability's safety (intersection between reachable set and unsafe set)
            visualize: visulization functions
            
## demo explaination:
The demo is done in the rounded-square mode of aa_planner (see below). The starting point is the bottom left point. The eight waypoints defined in the global     frame are (0, 0), (1, 0), (2, 1), (2, 2), (1, 3), (0, 3), (-1, 2), (-1, 1). https://youtu.be/VmE23Y745GI

The double-side red track is the safety boudary used by reachability analysis. Any intersection of the reachable set (green boxes) and the red line will be       predicted as dangerous, when fallback controller will be called. Note that in the video, the green box turns into red is not because of real collision, but       just having risk if we still insist apply our current bad control. The reachabilty verdict value is from the result of reachability analysis (positive means       safe, negative means unsafe). Same for the modelplex value. ModelPlex also outputs verdict ID. The detailed explanation of the ID is:
    
    -1: failed to reset time
    
    -2: controller modified speed measurement
    
    -3: control choice will violate lower speed limit
    
    -4: control choice will violate upper speed limit
    
    -5: control choice will make car go in reverse
    
    -6: control choice exceeds upper acceleration limit
    
    -7: control choice exceeds braking limit
    
    -8: brake configuration incompatible with speed limit configuration
    
    -9: acceleration configuration incompatible with speed limit configuration
    
    -10: invalid speed limit configuration (violates vl < vh)
    
    -11: invalid speed limit configuration (violates 0<=vl)
    
    -12: invalid goal (not on tube to origin)
    
    -13: invalid goal (not on tube to origin)
    
    -14: invalid tube (tube width exceeds curve radius)
    
    -15: invalid tube (goal not ahead of car)
    
    -16: invalid steering (steering not in direction of goal)
    
