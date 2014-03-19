c
c
      subroutine slancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)
c
c This routine carries out the same task as the routine  sti, but
c uses the more stable Lanczos method. The meaning of the input
c and output parameters is the same as in the routine  sti. (This
c routine is adapted from the routine RKPW in W.B. Gragg and
c W.J. Harrod,The numerically stable reconstruction of Jacobi
c matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)
c
      dimension x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap)
cf2py integer intent(in) :: n
cf2py integer intent(hide),depend(x) :: ncap=len(x)
cf2py real intent(in) :: x
cf2py real intent(in),depend(ncap),check(len(w)>=ncap) :: w
cf2py real intent(out),depend(n),dimension(n) :: alpha
cf2py real intent(out),depend(n),dimension(n) :: beta
cf2py real intent(hide),depend(ncap),dimension(ncap) :: p0
cf2py real intent(hide),depend(ncap),dimension(ncap) :: p1
cf2py integer intent(out) :: ierr
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        p0(i)=x(i)
        p1(i)=0.
   10 continue
      p1(1)=w(1)
      do 30 i=1,ncap-1
        pi=w(i+1)
        gam=1.
        sig=0.
        t=0.
        xlam=x(i+1)
        do 20 k=1,i+1
          rho=p1(k)+pi
          tmp=gam*rho
          tsig=sig
          if(rho.le.0.) then
            gam=1.
            sig=0.
          else
            gam=p1(k)/rho
            sig=pi/rho
          end if
          tk=sig*(p0(k)-xlam)-gam*t
          p0(k)=p0(k)-(tk-t)
          t=tk
          if(sig.le.0.) then
            pi=tsig*p1(k)
          else
            pi=(t**2)/sig
          end if
          tsig=sig
          p1(k)=tmp
   20   continue
   30 continue
      do 40 k=1,n
        alpha(k)=p0(k)
        beta(k)=p1(k)
   40 continue
      return
      end

