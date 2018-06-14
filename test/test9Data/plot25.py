#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np

pt.rc("font", size=16)
x = np.loadtxt("plot25.dat1")
x = x.reshape(  200,  200, order="F")
y = np.loadtxt("plot25.dat2")
f=open("plot25.dat3","r")
iss = 0
nd =     4
for i in range(nd):
    nn = int(f.readline())
    pt.fill(y[iss:iss+nn,0],y[iss:iss+nn,1],facecolor="w")
    iss=iss+nn

pt.imshow(x, extent=[-0.20, 1.20,-0.20, 0.70])
cb=pt.colorbar(shrink=.9)
pt.savefig("plot25.pdf")
