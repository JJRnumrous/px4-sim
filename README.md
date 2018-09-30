# PX4 Simulator
Simulator for PX4 SITL and HIL, written in C

## Setup
Run the following commands to pull the MAVLink library:
-   `git submodule init`
-   `git submodule update`

## Run
To run the simulation, run the following command:
-   `make clean && make && ./px4_quad_sim`

## Files
-   px4_quad_sim: the entry point and main loop of the simulator
-   quad_model: implementation of the quadcopter physics
-   sensor_model: implementation of the sensors of the quadcopter
-   px_sim_communication: implementation of the communication between the simulator and PX4

-   quad_parameters: the parameters of the quadcopter, sensors and environment
-   utils: utility functions needed by 1 or more of the above-mentioned files

## Customize Model
To customize the model, simply change the parameters in `quad-parameters.h`

## SITL vs HIL
To switch between SITL and HIL, change the constant `HIL` in `px4_quad_sim.h`
-   HIL: 0 - SITL
-   HIL: 1 - HIL

## Loop frequency
The loop frequency is denoted by the constant `SENSOR_FREQ` in `px4_quad_sim.h`.
This can be changed to change the frequency of the loop.