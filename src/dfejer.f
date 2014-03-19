      subroutine dfejer(n,x,w)
c
c This is a double-precision version of  fejer.
c
      double precision x,w,dpi,dn,dc1,dc0,dt,dsum,dc2
      dimension x(n),w(n)
cf2py integer intent(in) :: n
cf2py real*8 intent(out),depend(n),dimension(n) :: x
cf2py real*8 intent(out),depend(n),dimension(n) :: w
      dpi=4.d0*datan(1.d0)
      nh=n/2
      np1h=(n+1)/2
      dn=dble(n)
      do 10 k=1,nh
        x(n+1-k)=dcos(.5d0*dble(2*k-1)*dpi/dn)
        x(k)=-x(n+1-k)
   10 continue
      if(2*nh.ne.n) x(np1h)=0.d0
      do 30 k=1,np1h
        dc1=1.d0
        dc0=2.d0*x(k)*x(k)-1.d0
        dt=2.d0*dc0
        dsum=dc0/3.d0
        do 20 m=2,nh
          dc2=dc1
          dc1=dc0
          dc0=dt*dc1-dc2
          dsum=dsum+dc0/dble(4*m*m-1)
   20   continue
        w(k)=2.d0*(1.d0-2.d0*dsum)/dn
        w(n+1-k)=w(k)
   30 continue
      return
      end

