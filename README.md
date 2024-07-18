# NMPC Implementation 🕹️     

Model Predictive Control (MPC) uses a system model to predict future states based on optimal predicted inputs within a prediction horizon. The control applies only one input, then repeats the process to compensate for unmeasured noise or disturbances.

In this implementation, the system is an unmanned aerial vehicle (UAV) tracking a mobile vehicle. The cost function, minimized for predicted inputs, is derived from the distance between the UAV and the moving target. This code leverages the CasADi framework for NMPC.

## Repository Structure

This repository contains three language implementations: Python, C++, and MATLAB. The MATLAB implementation includes additional models such as state prediction of the target and dynamic obstacle avoidance modules.

- **"State predictive model of target"**: Contains the target's model, serving as a moving reference for the UAV and providing the initial cost value.
- **"NMPC_TT"**: Contains the NMPC code.

## UAV Tracking Target Illustration ✈️

### UAV without Gimbal
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/TrackWTG1.jpg">
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/TrackTWG.jpg">

---

### UAV with 3-DoF Gimbal
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/TargetTrack6.jpg">
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/TargetTrack5.jpg">
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/TargetTrack3.jpg">

---

## Without Obstacle Avoidance 🏢🏗️
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/Obstacle2.jpg">

## With Obstacle Avoidance 🏢🏗️
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/Obstacle3.jpg">
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/Imgs/obstacle5.jpg">

### Source Code Reference 🔗
This code is the source code of the paper: [NMPC-based UAV 3D Target Tracking In The Presence Of Obstacles and Visibility Constraints](https://ieeexplore.ieee.org/document/9476710).
