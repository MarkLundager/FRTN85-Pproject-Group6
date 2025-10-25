# Self-balancing Robot:

## Contributors

This project was done by:

Thomás Esteves  
Juan José Quintero Sarmiento  
Mark Peter Lundager  
Gabriel Carreria Figueiredo  
Simon Ziegenbalg  


## Hardware Setup
- A 6 RSS (Revolut-Spherical-Spherical) Stewart Platform.
   - 6 Servo-motors representing the 6 Revolute Joints
   - 2 discs (Base frame and Top frame)
   - 6 horns and 6 arms/links (The spherical-spherical joints)
- IMU-sensor to know the orientation of the platform
- An arduino to control the servo-motors and communicate with the IMU.
- External power source (Arduino cannot power the 6 servo motors by itself)

## Software setup:
- Python simulations to verify inverse kinematics, forward kinematics (currently not used), IMU and outer-loop PI Controller
- C++ version of the code for the arduino.

## Description of our implementation
Please see report in the repository.


