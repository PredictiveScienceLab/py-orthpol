c
c
      subroutine dlancz(n,ncap,x,w,alpha,beta,ierr,dp0,dp1)
c
c This is a double-precision version of the routine  lancz.
c
      double precision x(ncap),w(ncap),alpha(n),beta(n),
     *dp0(ncap),dp1(ncap),dpi,dgam,dsig,dt,xlam,drho,dtmp,
     *dtsig,dtk
cf2py integer intent(in) :: n
cf2py integer intent(hide),depend(x) :: ncap=len(x)
cf2py real*8 intent(in) :: x
cf2py real*8 intent(in),depend(ncap),check(len(w)>=ncap) :: w
cf2py real*8 intent(out,out=alpha),depend(n),dimension(n) :: alpha
cf2py real*8 intent(out,out=beta),depend(n),dimension(n) :: beta
cf2py real*8 intent(hide),depend(ncap),dimension(ncap) :: dp0
cf2py real*8 intent(hide),depend(ncap),dimension(ncap) :: dp1
cf2py integer intent(out) :: ierr
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        dp0(i)=x(i)
        dp1(i)=0.d0
   10 continue
      dp1(1)=w(1)
      do 30 i=1,ncap-1
        dpi=w(i+1)
        dgam=1.d0
        dsig=0.d0
        dt=0.d0
        xlam=x(i+1)
        do 20 k=1,i+1
          drho=dp1(k)+dpi
          dtmp=dgam*drho
          dtsig=dsig
          if(drho.le.0.d0) then
            dgam=1.d0
            dsig=0.d0
          else
            dgam=dp1(k)/drho
            dsig=dpi/drho
          end if
          dtk=dsig*(dp0(k)-xlam)-dgam*dt
          dp0(k)=dp0(k)-(dtk-dt)
          dt=dtk
          if(dsig.le.0.d0) then
            dpi=dtsig*dp1(k)
          else
            dpi=(dt**2)/dsig
          end if
          dtsig=dsig
          dp1(k)=dtmp
   20   continue
   30 continue
      do 40 k=1,n
        alpha(k)=dp0(k)
        beta(k)=dp1(k)
   40 continue
      return
      end

