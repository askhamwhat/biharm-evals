#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np

x = np.loadtxt("plot29.dat")
x = x.reshape(np.size(x)/2,2)
x1=x[:,0]
x2=x[:,1]
pt.plot(x1,x2,"kx")
pt.title("a                               ")
pt.show()
