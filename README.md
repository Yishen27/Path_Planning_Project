# Path_Planning_Project
---
This is a Udacity Self-Driving-Car Nano degree program project which aims at safely navigate around a virtual highway with other traffic 
that is driving +-10 MPH of the 50 MPH speed limit. The car's localization and sensor fusion data, and a sparse map list of waypoints 
around the highway are provided. With the path planning code, the car is is capable of utilizing sensor fusion data to calculate 
safety/efficiency/goal cost and generating the path which pass slow traffic by changing lanes while driving as fast as possible 
within the speed limit in complex highway condition. 

## Program performance
The simulated car was able to complete a 6946m highway loop as fast as possible without any accident.

## Main file structure
 
 ### Folder *src*
 This is the folder contains all c++ files constitute the program.
 
 - main.cpp 
 Main code that realizes the path planning. Most of the functions are defined here.
 
 - helpers.h
 All helper functions, include the cost calculations are defined here.
 
 - spline.h
 The spline libiary need for smoother path. 
 
 ### Project_readme.md
 Origianl readme for this project.
