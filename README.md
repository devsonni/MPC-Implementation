# NMPC Implementation đšī¸     

đ For simple intiution Model predictive control uses the system model to do prediction of future states of system for some predicted optimal inputs in prediction horizon, and applys only one input and does the same process again to compansate the unmeasured noise or disturbance in system.      
đ In this instance our system is unmaned arieal vehicle and it's following the target which is mobile vehicle.       
đ The cost function that sould be minimize for predicted inputs is made by the distance between UAV and target which is moving.   
đ This code uses the casadi framework to code NMPC.

â¨ "State predictive model of target" folder of repository has the targets model which feds for the UAV as reference(a moveing reference) and it gives the cost value from the initial point of UAV.      
â¨ "NMPC_TT" has the code of NMPC.      

## âī¸ UAV tracks the target illustration đĨ      

### đ UAV without Gimbal                
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TrackWTG1.jpg">            
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TrackTWG.jpg">     
   
##          

### đ UAV with 3-DoF Gimbal                
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TargetTrack6.jpg">         
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TargetTrack5.jpg">        
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TargetTrack3.jpg">      

## đĸđī¸ Without Obstacle Avoidance      
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Obstacle2.jpg">

## đĸđī¸ With Obstacle Avoidance       
<img align="right" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/Obstacle3.jpg">         
<img align="left" height="500" width="950" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/obstacle5.jpg">        


### đ This code is source code of this paper đ        
[NMPC-based UAV 3D Target Tracking In The Presence Of Obstacles and Visibility Constraints](https://ieeexplore.ieee.org/document/9476710)       



<!-- 
First           | Second
----------------|----------------------
<img align="left" height="500" width="600" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TargetTrack6.jpg">|<img align="left" height="500" width="600" src="https://github.com/devsonni/MPC-Implementation/blob/main/gif/TargetTrack6.jpg"> -->
