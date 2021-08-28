# NMPC Implementation 🕹️     

👉 For simple intiution Model predictive control uses the system model to do prediction of future states of system for some predicted optimal inputs in prediction horizon, and applys only one input and does the same process again to compansate the unmeasured noise or disturbance in system.      
👉 In this instance our system is unmaned arieal vehicle and it's following the target which is mobile vehicle.       
👉 The cost function that sould be minimize for predicted inputs is made by the distance between UAV and target which is moving.   
👉 This code uses the casadi framework to code NMPC.

✨ "State predictive model of target" folder of repository has the targets model which feds for the UAV as reference(a moveing reference) and it gives the cost value from the initial point of UAV.      
✨ "NMPC_TT" has the code of NMPC.      

## ✈️ UAV follows the target illustration 🔥      

### 📌 Target without angular velocity        
<!--img height="40" width="40" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Tracking1.gif"-->
<img align="middle" height="400" width="500" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Tracking1.gif">           


### 📌 Target with angular velocity                       
<!--img height="40" width="40" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Tracking5.gif"-->
<img align="middle" height="400" width="500" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Tracking5.gif">        


### 🔗 Basic source Code of this paper 📝        
[NMPC-based UAV 3D Target Tracking In The Presence Of Obstacles and Visibility Constraints](https://ieeexplore.ieee.org/document/9476710)
