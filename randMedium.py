#!/bin/env python
# generate bicontinuous medium from a given IGOM/msh/txt shape file, output msh and IGOM shape files
# @ Guanglin Tang, Texas A&M University, 2016.07.05
# tang2013@tamu.edu
from numpy import *
import sys
import subprocess
import datetime
from scipy.special import erfinv

#--------------------------------------------------------#
def tri2geo(tri): # tri[Ntri,3,3]
    eps = 1.e-5
    Ntri = len(tri)
    points = []
    lines = []
    faces = zeros((Ntri,3),dtype=int)
    # find points, faces
    for nTri in range(Ntri):
        for nPt in range(3):
            found = False
            for nPt1,point in enumerate(points):
                if ((tri[nTri,nPt]-point)**2).sum()<eps:
                    found = True
                    faces[nTri,nPt] = nPt1+1
                    break
            if not found:
                points += [tri[nTri,nPt].copy(),]
                faces[nTri,nPt] = len(points)
    faces0 = faces.copy()
    # find lines, faces
    for nTri in range(Ntri):
        for nPt in range(3):
            found = False
            for nln1,line in enumerate(lines):
                if faces0[nTri,nPt-1]==line[0] and faces0[nTri,nPt]==line[1]:
                    found = True
                    faces[nTri,nPt] = nln1+1
                    break
                if faces0[nTri,nPt-1]==line[1] and faces0[nTri,nPt]==line[0]:
                    found = True
                    faces[nTri,nPt] = -(nln1+1)
                    break
            if not found:
                lines += [[faces0[nTri,nPt-1],faces0[nTri,nPt]],]
                faces[nTri,nPt] = len(lines)
                
    # delete uncessary lines
    faces = faces.tolist()
    nln = 0
    while nln<len(lines):
        indx = array([len(where(abs(array(tmp))==nln+1)[0]) for tmp in faces])
        indx = where(indx==1)[0]
        cross1 = cross(points[lines[abs(faces[indx[0]][0])-1][0]-1]-points[lines[abs(faces[indx[0]][0])-1][1]-1],points[lines[abs(faces[indx[0]][1])-1][0]-1]-points[lines[abs(faces[indx[0]][1])-1][1]-1])
        cross2 = cross(points[lines[abs(faces[indx[1]][0])-1][0]-1]-points[lines[abs(faces[indx[1]][0])-1][1]-1],points[lines[abs(faces[indx[1]][1])-1][0]-1]-points[lines[abs(faces[indx[1]][1])-1][1]-1])
        if (cross(cross1,cross2)**2).sum()/(cross1**2).sum()/(cross2**2).sum()<eps:
            del lines[nln]
            if  nln+1 in faces[indx[0]]:    faces[indx[0]] = faces[indx[0]][faces[indx[0]].index( nln+1)+1:]+faces[indx[0]][:faces[indx[0]].index( nln+1)]
            if -nln-1 in faces[indx[0]]:    faces[indx[0]] = faces[indx[0]][faces[indx[0]].index(-nln-1)+1:]+faces[indx[0]][:faces[indx[0]].index(-nln-1)]
            if  nln+1 in faces[indx[1]]:    faces[indx[1]] = faces[indx[1]][faces[indx[1]].index( nln+1)+1:]+faces[indx[1]][:faces[indx[1]].index( nln+1)]
            if -nln-1 in faces[indx[1]]:    faces[indx[1]] = faces[indx[1]][faces[indx[1]].index(-nln-1)+1:]+faces[indx[1]][:faces[indx[1]].index(-nln-1)]
            faces[indx[0]] += faces[indx[1]]
            del faces[indx[1]]
            for face in faces:
                for nln1 in range(len(face)):
                    if abs(face[nln1]) > nln:
                        face[nln1] = (abs(face[nln1])-1)*sign(face[nln1])
            crt = inf
            if -faces[indx[0]][-1] in faces[indx[0]]:    crt = abs(faces[indx[0]][-1]); faces[indx[0]].remove(crt); faces[indx[0]].remove(-crt); del lines[crt-1]
            if -faces[indx[0]][-2] in faces[indx[0]]:    crt = abs(faces[indx[0]][-2]); faces[indx[0]].remove(crt); faces[indx[0]].remove(-crt); del lines[crt-1]
            for face in faces:
                for nln1 in range(len(face)):
                    if abs(face[nln1]) > crt:
                        face[nln1] = (abs(face[nln1])-1)*sign(face[nln1])
        else:
            nln += 1
    return array(points),array(lines),faces
#--------------------------------------------------------#

def writeGeo(points,lines,faces,fc,Npart,geoFil): # write geo file for gmsh
    f1 = open(geoFil,'w')
    f1.write('lc=' + str(fc) + ';\n' )
    f1.write('\n')
    for n,point in enumerate(points):
        f1.write('Point(' + str(n+1) + ')={' + ','.join(['%15.7e'%tmp for tmp in point]) + ',lc};\n' )
    f1.write('\n')
    for n,line in enumerate(lines):
        f1.write('Line('  + str(n+1) + ')={' + ','.join([    '%d'%tmp for tmp in line ]) + '};\n' )
    f1.write('\n')
    for n,face in enumerate(faces):
        f1.write('Line Loop('+str(n+1)+')={' + ','.join([   '%d'%tmp for tmp in face ]) + '};\n' )
    f1.write('\n')
    for n in range(len(faces)):
        f1.write('Plane Surface('+str(n+1)+')={'+str(n+1)+'};\n')
    f1.write('\n')
    for n in range(Npart):
        f1.write('Surface Loop(' + str(n+1) +')={' + ','.join(['%d'%(tmp+1) for tmp in range(len(faces)/Npart*n,len(faces)/Npart*(n+1))]) + '};\n' )
    f1.write('\n')
    for n in range(Npart):
        f1.write('Volume(' + str(n+1) + ')={' + str(n+1) + '};\n' )

#--------------------------------------------------------#

def writeGeoAxisy(points,fc,geoFil): # write geo file for gmsh for axisymmetric obj.
    N = len(points)
    f1 = open(geoFil,'w')
    f1.write('lc=' + str(fc) + ';\n' )
    f1.write('\n')
    f1.write('Point(' + str(1) + ')={' + ','.join(['%15.7e'%tmp for tmp in 0,0,points[0][1]]) + ',lc};\n' )
    for n,point in enumerate(points[1:-1]):
        f1.write('Point(' + str(1+5*n+1) + ')={' + ','.join(['%15.7e'%tmp for tmp in 0,0,point[1]]) + ',lc};\n' )
        f1.write('Point(' + str(1+5*n+2) + ')={' + ','.join(['%15.7e'%tmp for tmp in point[0],0,point[1]]) + ',lc};\n' )
        f1.write('Point(' + str(1+5*n+3) + ')={' + ','.join(['%15.7e'%tmp for tmp in 0,point[0],point[1]]) + ',lc};\n' )
        f1.write('Point(' + str(1+5*n+4) + ')={' + ','.join(['%15.7e'%tmp for tmp in -point[0],0,point[1]]) + ',lc};\n' )
        f1.write('Point(' + str(1+5*n+5) + ')={' + ','.join(['%15.7e'%tmp for tmp in 0,-point[0],point[1]]) + ',lc};\n' )
    f1.write('Point(' + str(1+5*(n+1)+1) + ')={' + ','.join(['%15.7e'%tmp for tmp in 0,0,points[-1][1]]) + ',lc};\n' )
    f1.write('\n')

    f1.write('Line('  + str(1) + ')={' + ','.join([    '%d'%tmp for tmp in 1,3 ]) + '};\n' )
    f1.write('Line('  + str(2) + ')={' + ','.join([    '%d'%tmp for tmp in 1,4 ]) + '};\n' )
    f1.write('Line('  + str(3) + ')={' + ','.join([    '%d'%tmp for tmp in 1,5 ]) + '};\n' )
    f1.write('Line('  + str(4) + ')={' + ','.join([    '%d'%tmp for tmp in 1,6 ]) + '};\n' )
    for n in range(len(points)-2):
        f1.write('Line('  + str(4+8*n+1) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+2,min(1+5*(n+1)+2,1+5*(N-2)+1) ]) + '};\n' )
        f1.write('Line('  + str(4+8*n+2) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+3,min(1+5*(n+1)+3,1+5*(N-2)+1) ]) + '};\n' )
        f1.write('Line('  + str(4+8*n+3) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+4,min(1+5*(n+1)+4,1+5*(N-2)+1) ]) + '};\n' )
        f1.write('Line('  + str(4+8*n+4) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+5,min(1+5*(n+1)+5,1+5*(N-2)+1) ]) + '};\n' )
        f1.write('Circle('  + str(4+8*n+5) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+2,1+5*n+1,1+5*n+3 ]) + '};\n' )
        f1.write('Circle('  + str(4+8*n+6) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+3,1+5*n+1,1+5*n+4 ]) + '};\n' )
        f1.write('Circle('  + str(4+8*n+7) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+4,1+5*n+1,1+5*n+5 ]) + '};\n' )
        f1.write('Circle('  + str(4+8*n+8) + ')={' + ','.join([    '%d'%tmp for tmp in 1+5*n+5,1+5*n+1,1+5*n+2 ]) + '};\n' )
    f1.write('\n')

    f1.write('Line Loop('+str(1)+')={' + ','.join([   '%d'%tmp for tmp in 1,9,-2 ]) + '};\n' )
    f1.write('Line Loop('+str(2)+')={' + ','.join([   '%d'%tmp for tmp in 2,10,-3 ]) + '};\n' )
    f1.write('Line Loop('+str(3)+')={' + ','.join([   '%d'%tmp for tmp in 3,11,-4 ]) + '};\n' )
    f1.write('Line Loop('+str(4)+')={' + ','.join([   '%d'%tmp for tmp in 4,12,-1 ]) + '};\n' )
    for n in range(len(points)-3):
        f1.write('Line Loop('+str(4+4*n+1)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*n+1,4+8*(n+1)+5,-4-8*n-2,-4-8*n-5 ]) + '};\n' )
        f1.write('Line Loop('+str(4+4*n+2)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*n+2,4+8*(n+1)+6,-4-8*n-3,-4-8*n-6 ]) + '};\n' )
        f1.write('Line Loop('+str(4+4*n+3)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*n+3,4+8*(n+1)+7,-4-8*n-4,-4-8*n-7 ]) + '};\n' )
        f1.write('Line Loop('+str(4+4*n+4)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*n+4,4+8*(n+1)+8,-4-8*n-1,-4-8*n-8 ]) + '};\n' )
    f1.write('Line Loop('+str(4+4*(n+1)+1)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*(n+1)+1,-4-8*(n+1)-2,-4-8*(n+1)-5 ]) + '};\n' )
    f1.write('Line Loop('+str(4+4*(n+1)+2)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*(n+1)+2,-4-8*(n+1)-3,-4-8*(n+1)-6 ]) + '};\n' )
    f1.write('Line Loop('+str(4+4*(n+1)+3)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*(n+1)+3,-4-8*(n+1)-4,-4-8*(n+1)-7 ]) + '};\n' )
    f1.write('Line Loop('+str(4+4*(n+1)+4)+')={' + ','.join([   '%d'%tmp for tmp in 4+8*(n+1)+4,-4-8*(n+1)-1,-4-8*(n+1)-8 ]) + '};\n' )
    f1.write('\n')

    for n in range(4*(N-1)):
        f1.write('Ruled Surface('+str(n+1)+')={'+str(n+1)+'};\n')
    f1.write('\n')

    f1.write('Surface Loop(' + str(1) +')={' + ','.join(['%d'%(tmp+1) for tmp in range(4*(N-1))]) + '};\n' )
    f1.write('\n')

    f1.write('Volume(' + str(1) + ')={' + str(1) + '};\n' )
#--------------------------------------------------------#
def readIgom(fil): # read igom shape file
    Ntri = subprocess.Popen(' '.join(['head -n 1',fil]),stdout=subprocess.PIPE,shell=True).communicate()[0].split()[0]
    Ntri = int(Ntri)
    tri = subprocess.Popen(' '.join(['tail -n',str(3*Ntri),fil]),stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    tri = array(tri,dtype=float).reshape(-1,3,3)
    return tri
#--------------------------------------------------------#

def readMsh(fil): # read msh shape file
    nlines = subprocess.Popen('grep -n Nodes '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = array(nodes,dtype=float).reshape(-1,4)[:,-3:]
    nlines = subprocess.Popen('grep -n Elements '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    faces = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split('\n')[:-1]
    for m1,face in enumerate(faces):
        if face.split()[1]=='2':
            break
    for m2,face in enumerate(faces[::-1]):
        if face.split()[1]=='2':
            break
    faces = array([tmp.split()[-3:] for tmp in faces[m1:len(faces)-m2]],dtype=int).tolist()
    n = 0
    while n<len(faces):
        if faces[n][0]==faces[n][1] or faces[n][0]==faces[n][2] or faces[n][1]==faces[n][2]:
            del faces[n]
        else:
            n += 1
    faces = array(faces)
    tri = nodes[faces-1,:]
    return tri
#--------------------------------------------------------#

def biContinue(fil,fv,kc,zp): # generate bi-continuous medium using msh file
    nlines = subprocess.Popen('grep -n Nodes '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = array(nodes,dtype=float).reshape(-1,4)[:,-3:]
    nlines = subprocess.Popen('grep -n Elements '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    faces = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split('\n')[:-1]
    for m1,face in enumerate(faces):
        if face.split()[1]=='4':
            break
    for m2,face in enumerate(faces[::-1]):
        if face.split()[1]=='4':
            break
    faces = array([tmp.split()[1:] for tmp in faces[m1:len(faces)-m2]],dtype=int)
    x = (nodes[faces[:,-4]-1,0]+nodes[faces[:,-3]-1,0]+nodes[faces[:,-2]-1,0]+nodes[faces[:,-1]-1,0])/4
    y = (nodes[faces[:,-4]-1,1]+nodes[faces[:,-3]-1,1]+nodes[faces[:,-2]-1,1]+nodes[faces[:,-1]-1,1])/4
    z = (nodes[faces[:,-4]-1,2]+nodes[faces[:,-3]-1,2]+nodes[faces[:,-2]-1,2]+nodes[faces[:,-1]-1,2])/4
    N = len(faces)
    Ntri = int(nlines[1])-int(nlines[0])-2-len(faces)
    Nwave = 10000
    sp = zp+1
    lp = sp/kc
    rot=2
    cnorm = sqrt(Nwave)
    level = erfinv(1 - 2*fv)
    twopi = 2*pi
    randmp = random.random((Nwave,3))
    ransk = random.gamma(sp,1/lp,Nwave)
    amp = zeros((N))
    icella = zeros((N))
    for n in range(Nwave):
        sk=ransk[n]
        delta=twopi*randmp[n,0]
        mu=(randmp[n,1]-0.5)*rot
        smu=sqrt(1.0-mu*mu)
        phi=twopi*randmp[n,2]
        knx=sk*smu*cos(phi)
        kny=sk*smu*sin(phi)
        knz=sk*mu
        arg=knx*x+kny*y+knz*z+delta
        amp=amp+cos(arg)
    pa = where(amp/cnorm>level)[0]
    Nice = len(pa)
    icella[pa] = 1
    subprocess.Popen('head -n '+nlines[0]+' '+fil+' > '+fil+'_bic',stdout=subprocess.PIPE,shell=True).communicate()[0]
    f1 = open(fil+'_bic','a')
    f1.write(str(Ntri+Nice)+'\n')
    f1.close()
    subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[0])+1+Ntri)+'p '+fil+' >> '+fil+'_bic',stdout=subprocess.PIPE,shell=True).communicate()[0]
    f1 = open(fil+'_bic','a')
    for n,npa in enumerate(pa):
        f1.write(' '.join([str(tmp) for tmp in [Ntri+1+n,]+list(faces[npa,:])])+'\n')
    f1.write('$EndElements')
    f1.close()
#--------------------------------------------------------#


def msh2igom1(fil): # msh file surface to igom shape file
    nlines = subprocess.Popen('grep -n Nodes '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = array(nodes,dtype=float).reshape(-1,4)[:,-3:]
    nlines = subprocess.Popen('grep -n Elements '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    faces = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split('\n')[:-1]
    for m1,face in enumerate(faces):
        if face.split()[1]=='2':
            break
    for m2,face in enumerate(faces[::-1]):
        if face.split()[1]=='2':
            break
    faces = array([tmp.split()[-3:] for tmp in faces[m1:len(faces)-m2]],dtype=int)
    print str(len(faces))
    for n in faces.reshape(-1):
        print ' '.join(['%15.7e'%tmp for tmp in nodes[n-1]])
#--------------------------------------------------------#

def msh2igom2(fil): # msh file volume to igom shape file
    nlines = subprocess.Popen('grep -n Nodes '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    nodes = array(nodes,dtype=float).reshape(-1,4)[:,-3:]
    nlines = subprocess.Popen('grep -n Elements '+fil+'|cut -d : -f 1',stdout=subprocess.PIPE,shell=True).communicate()[0].split()
    volums = subprocess.Popen('sed -n '+str(int(nlines[0])+2)+','+str(int(nlines[1])-1)+'p '+fil,stdout=subprocess.PIPE,shell=True).communicate()[0].split('\n')[:-1]
    for m1,face in enumerate(volums):
        if face.split()[1]=='4':
            break
    for m2,face in enumerate(volums[::-1]):
        if face.split()[1]=='4':
            break
    volums = array([tmp.split()[-4:] for tmp in volums[m1:len(volums)-m2]],dtype=int)
    faces = []
    faces1 = []
    for volum in volums:
        for newFace in [volum[[0,1,2]],volum[[2,1,3]],volum[[2,3,0]],volum[[3,1,0]]]:
            newFace = roll(newFace,-newFace.argmin()).tolist()
            faces.append(newFace)
            faces1.append(sorted(newFace))
    order = [i[0] for i in sorted(enumerate(faces1), key=lambda x:x[1])]
    faces = array(faces)[order].tolist()
    faces1 = array(faces1)[order].tolist()
    n = 0
    while n<len(faces)-1:
        if faces1[n]==faces1[n+1]:
            del faces[n:n+2]
            del faces1[n:n+2]
        else:
            n += 1
    print str(len(faces))
    for face in faces:
        for n in face:
            print ' '.join(['%15.7e'%tmp for tmp in nodes[n-1]])
#--------------------------------------------------------#
    

if __name__=='__main__':
## python randMedium.py /scratch/user/tang2013/Data/particleShape/hexgon.dat 0.2 1 0.5 100 1 >aaa.dat
    fil = '/scratch/user/tang2013/Data/particleShape/hexgon.dat'
    lc = 0.2
    Npart = 1
    fv = 0.5
    kc = 100.
    zp = 1.
    if len(sys.argv)>1:
        fil = sys.argv[1]
    if len(sys.argv)>2:
        lc = float(sys.argv[2])
    if len(sys.argv)>3:
        Npart = int(sys.argv[3])
    if len(sys.argv)>4:
        fv = float(sys.argv[4])
    if len(sys.argv)>5:
        kc = float(sys.argv[5])
    if len(sys.argv)>6:
        zp = float(sys.argv[6])

    geoFil = datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y%m%d')
    #geoFil = '/scratch/user/tang2013/tmp/'+datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y%m%d%H%M%S%f')
    if fil[-3:]=='dat':
        points,lines,faces = tri2geo(readIgom(fil))
        writeGeo(points,lines,faces,lc,Npart,geoFil+'.geo')
    elif fil[-3:]=='msh':
        points,lines,faces = tri2geo(readMsh(fil))
        writeGeo(points,lines,faces,lc,Npart,geoFil+'.geo')
    elif fil[-3:]=='txt':
        points = loadtxt(fil,skiprows=1)
        writeGeoAxisy(points,lc,geoFil+'.geo')
    subprocess.Popen('gmsh '+geoFil+'.geo -3',stdout=subprocess.PIPE,shell=True).communicate()
    biContinue(geoFil+'.msh',fv,kc,zp)
    #msh2igom1(geoFil+'.msh')
    msh2igom2(geoFil+'.msh_bic')
