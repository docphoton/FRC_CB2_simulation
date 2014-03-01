FRC_CB2_simulation
==================

Physics Simulations of robot/ball dynamics

    Included are C-files that perform numerical simulations of various dynamical 
    aspects of the FRC2014 challenge and the FRC team 4601 Robot's solution to those challenges. 
    
    ball_simulation.c  - solves the equations of motion for  a spherical ball, not spinning, projected 
           through air (includes air drag). Uses gsl libraries. 
    
    ball_launch.c - solves the equations for the ball rolling off the catapult as it is thrown.
           Assumes that (1) the ball remains a rigid body and (2) rolls without slipping as it is thrown.
           uses gsl libraries. 
