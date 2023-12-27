# Earth-Luna-CR3BP
###### Jack McKinley

## Project Description
Given an initial state vector (multiple are provided), this MATLAB code propagates orbits in the Circular Restricted Earth-Luna system by MATLAB's numeric ordinary differential equation solver, ODE45. The code also numerically solves all its Lagrange points via the Newton-Raphson method. 

## How to Use
All code is written and ran in MATLAB.

Under "Initializations", multiple initial state vectors representing different orbit types have already been provided, as well as an empty vector as a template. Simply comment/uncomment the desired vector and run the code to receive a three dimensional plot of the orbit. The period of time the orbit can be propagated across can be easily adjusted with variable "t_final". Locations of the colinear Lagrange points are also outputted to the user in terms of characteristic length. 

## Credits
Equations of Motion derived from Professor Kathleen Howell's AAE 632 "Advanced Orbital Dynamics" lectures.
