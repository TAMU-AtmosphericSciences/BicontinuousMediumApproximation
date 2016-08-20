#!/bin/env python
# generate graupel shape file
# @ Guanglin Tang, Texas A&M University, 2016.08.05
# tang2013@tamu.edu
# Ref: Wang PaoKuan, 1982
from numpy import *
import sys

def graupel(C,a,lamb,N):
    th = linspace(-1.,1.,N)*pi/2.
    x = zeros(th.shape)
    y = zeros(th.shape)
    print str(N)
    for n in range(len(th)):
        y[n] = C*sin(th[n])
        x[n] = a*(1.0-(y[n]/C)**2)**0.5*arccos(y[n]/lamb/C)
        print ' '.join('%16.8e'%tmp for tmp in [x[n],y[n]])

if __name__=='__main__':
# python graupel.py 1.0 0.5 1.2
    C = 1.0
    a = 0.5
    lamb = 1.2
    N = 100
    if len(sys.argv)>1:
        C = float(sys.argv[1])
    if len(sys.argv)>2:
        a = float(sys.argv[2])
    if len(sys.argv)>3:
        lamb = float(sys.argv[3])
    if len(sys.argv)>4:
        N = int(sys.argv[4])
    graupel(C,a,lamb,N)
