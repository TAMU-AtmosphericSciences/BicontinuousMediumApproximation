#!/bin/env python
# @ Guanglin Tang, Texas A&M University, 2016.07.12
# tang2013@tamu.edu

## convert .dat file to .stl file
from numpy import *
from randMedium import readIgom

def dat2stl(file1):
    tri = readIgom(file1)
    print 'solid Untitled'
    for n in range(len(tri)):
        norm = cross(tri[n,0]-tri[n,1],tri[n,0]-tri[n,2])
        norm = norm/(norm**2).sum()**0.5
        print 'facet normal',' '.join([str(tmp) for tmp in norm])
        print '  outer loop'
        print '    vertex',' '.join([str(tmp) for tmp in tri[n,0]])
        print '    vertex',' '.join([str(tmp) for tmp in tri[n,1]])
        print '    vertex',' '.join([str(tmp) for tmp in tri[n,2]])
        print '  endloop'
        print 'endfacet'
    print 'endsolid Untitled'

if __name__=='__main__':
    if len(sys.argv)>1:
        file1 = sys.argv[1]
    dat2stl(file1)
