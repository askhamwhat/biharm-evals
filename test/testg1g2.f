
      program testg1g2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        this subroutine tests the greens function
c        in the farkas representation for the helmholtz stokes
c        equation with the formula in the jiang-quaife-kropinski
c        paper
cccccc

      implicit real *8 (a-h,o-z)
      real *8 dtoobig, dtoosmall
      complex *16 eye
      data eye / (0.0d0,1.0d0) /
      complex *16 zk
      real *8 rr,rnsrc(2),rntarg(2),source(2),targ(2)
      complex *16 val,grad(2),hess(4),val2
     
      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n

      done = 1
      pi = atan(done)*4



      zk = 1 + hkrand(0)
      source(1) = hkrand(0)
      source(2) = hkrand(0)

      targ(1) = hkrand(0) + 2.0d0
      targ(2) = hkrand(0)


      thet = 2.0d0*pi*hkrand(0)
      rnsrc(1) = cos(thet) 
      rnsrc(2) = sin(thet) 


      thet = 2.0d0*pi*hkrand(0)
      rntarg(1) = cos(thet) 
      rntarg(2) = sin(thet)

cc      call prin2('rnsrc=*',rnsrc,2)


      call zhfark_ck1(zk,source,targ,rnsrc,pars2,val,
     1   grad,hess)
      
      call helmstokes_fark_g1_jiang(zk,source,targ,rnsrc,val2)
      

cc      call prin2('val=*',val,2)
cc      call prin2('val2=*',val2,2)

      err = abs(val-val2)
      call prin2('error kernel1=*',err,1)


      val = 0
      val2 = 0
      call zhfark_ck2(zk,source,targ,rnsrc,pars2,val,
     1   grad,hess)
      
      call helmstokes_fark_g2_jiang(zk,source,targ,rnsrc,val2)
      

cc      call prin2('val=*',val,2)
cc      call prin2('val2=*',val2,2)

      err = abs(val-val2)
      call prin2('error kernel2=*',err,1)

      

      stop
      end

c------------------------------------------------

      subroutine helmstokes_fark_g2_jiang(zk,src,targ,rnsrc,val)
      implicit real *8 (a-h,o-z)
      complex *16 zk,val,zk2,hanks(0:100),ima,ztmp
      data ima/(0.0d0,1.0d0)/
      real *8 src(2),targ(2),rnsrc(2)


      done = 1
      pi = atan(done)*4

cc      call prin2('src=*',src,2)

      dx = targ(1) - src(1)
      dy = targ(2) - src(2)

      rr = sqrt(dx**2 + dy**2)
      rrn = dx*rnsrc(1) + dy*rnsrc(2)
      
      zk2 = rr*zk

cc      call prin2('zk2=*',zk2,2)
cc      call prin2('zk=*',zk,2)

      ifexpon = 1
      call hanks104(zk2,hanks(0),5,ifexpon)

cc      call prin2('hanks=*',hanks(0),10)

      ztmp = pi*ima/2*(hanks(0)-2/zk2*hanks(1)) + 2/zk2**2
      val = -ztmp*(0.5d0-rrn**2/rr**2)/pi

      return
      end
c---------------------------------      

      subroutine helmstokes_fark_g1_jiang(zk,src,targ,rnsrc,val)
      implicit real *8 (a-h,o-z)
      complex *16 zk,val,zk2,hanks(0:100),ima,ztmp,ztmp2,ztmp3
      data ima/(0.0d0,1.0d0)/
      real *8 src(2),targ(2),rnsrc(2)


      done = 1
      pi = atan(done)*4

cc      call prin2('src=*',src,2)

      dx = targ(1) - src(1)
      dy = targ(2) - src(2)

      rr = sqrt(dx**2 + dy**2)
      rrn = dx*rnsrc(1) + dy*rnsrc(2)
      
      zk2 = rr*zk

cc      call prin2('zk2=*',zk2,2)
cc      call prin2('zk=*',zk,2)

      ifexpon = 1
      call hanks104(zk2,hanks(0),5,ifexpon)

cc      call prin2('hanks=*',hanks(0),10)

      ztmp = pi*ima/2*(hanks(0)-2/zk2*hanks(1)) + 2/zk2**2
      ztmp2 = 3*ztmp + pi*ima/2*zk2*hanks(1) + 0.5d0
      ztmp3 = 4*ztmp + pi*ima/2*zk2*hanks(1)

cc      call prin2('ztmp=*',ztmp,2)
cc      call prin2('ztmp2=*',ztmp2,2)
cc      call prin2('ztmp3=*',ztmp3,2)
      val = -rrn/pi/rr**2*ztmp2 + rrn**3/pi/rr**4*ztmp3 

      return
      end
c---------------------------------      
