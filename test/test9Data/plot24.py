#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np

x = np.loadtxt("plot24.dat1")
x = x.reshape(np.size(x)/2,2)
x1=x[:,0]
x2=x[:,1]

y = np.loadtxt("plot24.dat2")
y = y.reshape(np.size(y)/2,2)
y1=y[:,0]
y2=y[:,1]

pt.plot(x1,x2,"k.")
pt.plot(y1,y2,"gx")
pt.title("a                               ")
pt.axes().set_aspect("equal")
pt.show()
