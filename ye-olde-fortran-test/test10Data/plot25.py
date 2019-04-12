#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np

pt.rc("font", size=16)
x = np.loadtxt("plot25.dat1")
x = x.reshape(  300,  300, order="F")
y = np.loadtxt("plot25.dat2")
f=open("plot25.dat3","r")
iss = 0
nd =     1
for i in range(nd):
    nn = int(f.readline())
    pt.fill(y[iss:iss+nn,0],y[iss:iss+nn,1],facecolor="w")
    iss=iss+nn

pt.imshow(x, extent=[-1.87, 1.87,-1.87, 1.87])
cb=pt.colorbar(shrink=.9)
pt.savefig("plot25.pdf")
