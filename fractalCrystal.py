#!/bin/env python
# @ Guanglin Tang, Texas A&M University, 2016.07.12
# tang2013@tamu.edu
# Shape defined in Yang and Liou 1998. Figure 1(b) left

#            
#         __/ \__
#         \     /
#   __/\__/     \__/\__
#   \                 /
#   /_               _\
#     \             /
#   __/             \__
#   \                 /
#   /_  __       __  _\
#     \/  \     /  \/
#         /_   _\
#           \ /
#            


from numpy import *
import matplotlib.pyplot as plt

beta = 30.*pi/180.
a1 = .5
a2 = .5
L = .2

points1 = array([[a1+a2,0],
          [a1+a2-a1*tan(beta)/sqrt(3.),a1*tan(beta)],
          [a1,a1*tan(beta)],
          [a1+(a2-a1*tan(beta)/sqrt(3.))/2.,a1*tan(beta)+(a2-a1*tan(beta)/sqrt(3.))*sqrt(3.)/2.]])
points = points1.copy()
for th in arange(60.,360.,60.)*pi/180.:
    points = append(points,array([points1[:,0]*cos(th)-points1[:,1]*sin(th),points1[:,0]*sin(th)+points1[:,1]*cos(th)]).T,axis=0)

Np = len(points)
points = points.tolist()
print Np*4
## top face
for np in range(Np):
    print ' '.join(['%15.8e'%tmp for tmp in [0.,0.,L/2.]])
    print ' '.join(['%15.8e'%tmp for tmp in points[mod(np+1,Np)]+[L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[np]+[L/2.,]])
## bottom face
for np in range(Np):
    print ' '.join(['%15.8e'%tmp for tmp in [0.,0.,-L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[np]+[-L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[mod(np+1,Np)]+[-L/2.,]])
## side face
for np in range(Np):
    print ' '.join(['%15.8e'%tmp for tmp in points[np]+[L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[mod(np+1,Np)]+[-L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[np]+[-L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[np]+[L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[mod(np+1,Np)]+[L/2.,]])
    print ' '.join(['%15.8e'%tmp for tmp in points[mod(np+1,Np)]+[-L/2.,]])
