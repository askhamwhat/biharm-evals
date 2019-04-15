from numpy import *
from pylab import *
import scipy.io as sio

mat_content = sio.loadmat('mat-files/example_annulus_speedtest.mat')
zkvals = mat_content['zkvals']
nchvals = mat_content['nchvals']
r1 = 1
r2 = 1.7
nchovals = ceil(nchvals/r1*r2)+1
n = transpose(nchvals + nchovals)*16
tt = mat_content['tt']

llll = ["{:.2f}".format(zkvals[0,0]),"{:.2f}".format(zkvals[0,4]),
"{:.2f}".format(zkvals[0,9])]

figure(1)
loglog(n,tt[0,:],'k.',markersize=7,label="k="+llll[0])
loglog(n,tt[4,:],'k^',markersize=5,label="k="+llll[1])
loglog(n,tt[9,:],'ks',markersize=5,label="k="+llll[2])
loglog(n,n*log(n)/2500,'k--',label=r"$O(N\log{N})$",linewidth=0.5)
legend(loc='upper left')
xticks([250,500,1000,2000,4000],[r'$2.5 \times 10^{2}$',r'$5.0 \times 10^{2}$',r'$1.0\times10^{3}$',r'$2.0 \times 10^{3}$',r'$4.0\times 10^{3}$'])
savefig('speed_res.pdf',bbox_inches='tight')
