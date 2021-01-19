# ReachModelPlex
Project of integrating reachability and ModelPlex Monitor

1. All python dependencies

    matplotlib==3.0.3

    numpy==1.18.2

    scipy==1.4.1

    cffi==1.14.2

    mpmath==1.1.0

    pyinterval==1.2.1.dev0

    pycddlib==2.1.2

    cvxpy==1.1.5

2. set up reachability tool
    2.1. put "flowstar" file and model file (nonlinear_without_beta) file in the same level of directory
    2.2. install gcc-8 and g++-8, so far flow* only works with at least version-8 compilers
    2.3. install latest m4, gmp, mpfr, gsl, glpk, bison, and flex. (one possible issue in Ubuntu 18, please install flex 2.6.3 instaed of 2.6.4)
    2.4. type "make" in the flowstar file, then type "sudo ldconfig"
    2.5. type "make" in the model file, if everything goes well, you'll find an executable RC_bicycle program
    NOTE: so far the prototype is only working in Linux, only Ubuntu16.04 LTS has been tested.

3. To run the simulation, we just execute the main python file: simulator.py

4. Introduction of all subfolders:
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
            
5. demo explaination:
    the demo is done in the rounded-square mode of aa_planner (see below). The starting point is the bottom left point. The eight waypoints defined in the global     frame are (0, 0), (1, 0), (2, 1), (2, 2), (1, 3), (0, 3), (-1, 2), (-1, 1). https://youtu.be/VmE23Y745GI
