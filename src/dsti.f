c
c
      subroutine dsti(n,ncap,x,w,alpha,beta,ierr,dp0,dp1,dp2)
c
c This is a double-precision version of the routine  sti.
c
      double precision x,w,alpha,beta,dp0,dp1,dp2,dtiny,d1mach,
     *dhuge,dsum0,dsum1,dsum2,dt
      dimension x(ncap),w(ncap),alpha(n),beta(n),dp0(ncap),
     *dp1(ncap),dp2(ncap)
cf2py integer intent(in) :: n
cf2py integer intent(hide),depend(x) :: ncap=len(x)
cf2py real*8 intent(in) :: x
cf2py real*8 intent(in),depend(ncap),check(len(w)>=ncap) :: w
cf2py real*8 intent(out,out=alpha),depend(n),dimension(n) :: alpha
cf2py real*8 intent(out,out=beta),depend(n),dimension(n) :: beta
cf2py real*8 intent(hide),depend(ncap),dimension(ncap) :: dp0
cf2py real*8 intent(hide),depend(ncap),dimension(ncap) :: dp1
cf2py real*8 intent(hide),depend(ncap),dimension(ncap) :: dp2
cf2py integer intent(out) :: ierr
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
      dsum0=0.d0
      dsum1=0.d0
      do 10 m=1,ncap
        dsum0=dsum0+w(m)
        dsum1=dsum1+w(m)*x(m)
   10 continue
      alpha(1)=dsum1/dsum0
      beta(1)=dsum0
      if(n.eq.1) return
      do 20 m=1,ncap
        dp1(m)=0.d0
        dp2(m)=1.d0
   20 continue
      do 40 k=1,nm1
        dsum1=0.d0
        dsum2=0.d0
        do 30 m=1,ncap
          if(w(m).eq.0.d0) goto 30
          dp0(m)=dp1(m)
          dp1(m)=dp2(m)
          dp2(m)=(x(m)-alpha(k))*dp1(m)-beta(k)*dp0(m)
          if(dabs(dp2(m)).gt.dhuge .or. dabs(dsum2).gt.dhuge) then
            ierr=k
            return
          end if
          dt=w(m)*dp2(m)*dp2(m)
          dsum1=dsum1+dt
          dsum2=dsum2+dt*x(m)
   30   continue
        if(dabs(dsum1).lt.dtiny) then
          ierr=-k
          return
        end if
        alpha(k+1)=dsum2/dsum1
        beta(k+1)=dsum1/dsum0
        dsum0=dsum1
   40 continue
      return
      end

