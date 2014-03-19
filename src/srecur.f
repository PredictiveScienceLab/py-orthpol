c
c
      subroutine srecur(n,ipoly,al,be,a,b,ierr)
c
c This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
c in the recurrence relation
c
c       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                            k=0,1,...,n-1,
c
c       p(-1)(x)=0,  p(0)(x)=1,
c
c for some classical (monic) orthogonal polynomials, and sets  b(0)
c equal to the total mass of the weight distribution. The results are
c stored in the arrays  a,b,  which hold, respectively, the coefficients
c a(k-1),b(k-1), k=1,2,...,n.
c
c       Input:  n - - the number of recursion coefficients desired
c               ipoly-integer identifying the polynomial as follows:
c                     1=Legendre polynomial on (-1,1)
c                     2=Legendre polynomial on (0,1)
c                     3=Chebyshev polynomial of the first kind
c                     4=Chebyshev polynomial of the second kind
c                     5=Jacobi polynomial with parameters  al=-.5,be=.5
c                     6=Jacobi polynomial with parameters  al,be
c                     7=generalized Laguerre polynomial with
c                       parameter  al
c                     8=Hermite polynomial
c               al,be-input parameters for Jacobi and generalized
c                     Laguerre polynomials
c
c       Output: a,b - arrays containing, respectively, the recursion
c                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
c               ierr -an error flag, equal to  0  on normal return,
c                     equal to  1  if  al  or  be  are out of range
c                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0)
c                     overflows when  ipoly=6  or  ipoly=7, equal to  3
c                     if  n  is out of range, and equal to  4  if  ipoly
c                     is not an admissible integer. In the case  ierr=2,
c                     the coefficient  b(0)  is set equal to the largest
c                     machine-representable number.
c
c The subroutine calls for the function subroutines  r1mach,gamma  and
c alga. The routines  gamma  and  alga , which are included in this
c file, evaluate respectively the gamma function and its logarithm for
c positive arguments. They are used only in the cases  ipoly=6  and
c ipoly=7.
c
      external gamma
      dimension a(n),b(n)
cf2py integer intent(in) :: n
cf2py integer intent(in) :: ipoly
cf2py real optional,intent(in) :: al=-.5
cf2py real optional,intent(in) :: be=.5
cf2py real intent(out),depend(n),dimension(n) :: a
cf2py real intent(out),depend(n),dimension(n) :: b
cf2py integer intent(out) :: ierr
      if(n.lt.1) then
        ierr=3
        return
      end if
      almach=alog(r1mach(2))
      ierr=0
      do 10 k=1,n
        a(k)=0.
   10 continue
      if(ipoly.eq.1) then
        b(1)=2.
        if (n.eq.1) return
        do 20 k=2,n
          fkm1=real(k-1)
          b(k)=1./(4.-1./(fkm1*fkm1))
   20   continue
        return
      else if (ipoly.eq.2) then
        a(1)=.5
        b(1)=1.
        if(n.eq.1) return
        do 30 k=2,n
          a(k)=.5
          fkm1=real(k-1)
          b(k)=.25/(4.-1./(fkm1*fkm1))
   30   continue
        return
      else if(ipoly.eq.3) then
        b(1)=4.*atan(1.)
        if(n.eq.1) return
        b(2)=.5
        if(n.eq.2) return
        do 40 k=3,n
          b(k)=.25
   40   continue
        return
      else if(ipoly.eq.4) then
        b(1)=2.*atan(1.)
        if(n.eq.1) return
        do 50 k=2,n
          b(k)=.25
   50   continue
        return
      else if(ipoly.eq.5) then
        b(1)=4.*atan(1.)
        a(1)=.5
        if(n.eq.1) return
        do 60 k=2,n
          b(k)=.25
   60   continue
        return
      else if(ipoly.eq.6) then
        if(al.le.-1. .or. be.le.-1.) then
          ierr=1
          return
        else
          alpbe=al+be
          a(1)=(be-al)/(alpbe+2.)
          t=(alpbe+1.)*alog(2.)+alga(al+1.)+alga(be+1.)-
     *      alga(alpbe+2.)
          if(t.gt.almach) then
            ierr=2
            b(1)=r1mach(2)
          else
            b(1)=exp(t)
          end if
          if(n.eq.1) return
          al2=al*al
          be2=be*be
          a(2)=(be2-al2)/((alpbe+2.)*(alpbe+4.))
          b(2)=4.*(al+1.)*(be+1.)/((alpbe+3.)*(alpbe+2.)**2)
          if(n.eq.2) return
          do 70 k=3,n
            fkm1=real(k-1)
            a(k)=.25*(be2-al2)/(fkm1*fkm1*(1.+.5*alpbe/fkm1)*
     *        (1.+.5*(alpbe+2.)/fkm1))
            b(k)=.25*(1.+al/fkm1)*(1.+be/fkm1)*(1.+alpbe/fkm1)/
     *      ((1.+.5*(alpbe+1.)/fkm1)*(1.+.5*(alpbe-1.)/fkm1)
     *      *(1.+.5*alpbe/fkm1)**2)
   70     continue
          return
        end if
      else if(ipoly.eq.7) then
        if(al.le.-1.) then
          ierr=1
          return
        else
          a(1)=al+1.
          b(1)=gamma(al+1.,ierr)
          if(ierr.eq.2) b(1)=r1mach(2)
          if(n.eq.1) return
          do 80 k=2,n
            fkm1=real(k-1)
            a(k)=2.*fkm1+al+1.
            b(k)=fkm1*(fkm1+al)
   80     continue
          return
        end if
      else if(ipoly.eq.8) then
        b(1)=sqrt(4.*atan(1.))
        if(n.eq.1) return
        do 90 k=2,n
          b(k)=.5*real(k-1)
   90   continue
        return
      else
        ierr=4
      end if
      end

      function alga(x)
c
c This is an auxiliary function subroutine (not optimized in any
c sense) evaluating the logarithm of the gamma function for positive
c arguments  x. It is called by the subroutine  gamma. The integer  m0
c in the first executable statement is the smallest integer  m  such
c that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
c largest machine-representable number. The routine is based on a
c rational approximation valid on [.5,1.5] due to W.J. Cody and
c K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the
c case  n=7  in Table II. For the computation of  m0  it calls upon the
c function subroutines  t  and  r1mach. The former, appended below,
c evaluates the inverse function  t = t(y)  of  y = t ln t.
c
      dimension cnum(8),cden(8)
      data cnum/4.120843185847770,85.68982062831317,243.175243524421,
     *-261.7218583856145,-922.2613728801522,-517.6383498023218,
     *-77.41064071332953,-2.208843997216182/,
     *cden/1.,45.64677187585908,377.8372484823942,951.323597679706,
     *846.0755362020782,262.3083470269460,24.43519662506312,
     *.4097792921092615/
c
c The constants in the statement below are  exp(1.)  and  .5*alog(8.).
c
      m0=2.71828*t((alog(r1mach(2))-1.03972)/2.71828)
      xi=aint(x)
      if(x-xi.gt..5) xi=xi+1.
      m=ifix(xi)-1
c
c Computation of log gamma on the standard interval (1/2,3/2]
c
      xe=x-real(m)
      snum=cnum(1)
      sden=cden(1)
      do 10 k=2,8
        snum=xe*snum+cnum(k)
        sden=xe*sden+cden(k)
   10 continue
      alga=(xe-1.)*snum/sden
c
c Computation of log gamma on (0,1/2]
c
      if(m.eq.-1) then
        alga=alga-alog(x)
        return
      else if(m.eq.0) then
        return
      else
c
c Computation of log gamma on (3/2,5/2]
c
        p=xe
        if(m.eq.1) then
          alga=alga+alog(p)
          return
        else
c
c Computation of log gamma for arguments larger than 5/2
c
          mm1=m-1
c
c The else-clause in the next statement is designed to avoid possible
c overflow in the computation of  p  in the if-clause, at the expense
c of computing many logarithms.
c
          if(m.lt.m0) then
            do 20 k=1,mm1
              p=(xe+real(k))*p
   20       continue
            alga=alga+alog(p)
            return
          else
            alga=alga+alog(xe)
            do 30 k=1,mm1
              alga=alga+alog(xe+real(k))
   30       continue
            return
          end if
        end if
      end if
      end

      function gamma(x,ierr)
c
c This evaluates the gamma function for real positive  x, using the
c function subroutines  alga  and  r1mach. In case of overflow, the
c routine returns the largest machine-representable number and the
c error flag  ierr=2.
c
      almach=alog(r1mach(2))
      ierr=0
      t=alga(x)
      if(t.ge.almach) then
        ierr=2
        gamma=r1mach(2)
        return
      else
        gamma=exp(t)
        return
      end if
      end

      function t(y)
c
c This evaluates the inverse function  t = t(y)  of y = t ln t  for
c nonnegative  y  to an accuracy of about one percent. For the
c approximation used, see pp. 51-52 in W. Gautschi,Computational
c aspects of three-term recurrence relations'', SIAM Rev. 9, 1967,
c 24-82.
c
      if(y.le.10.) then
        p=.000057941*y-.00176148
        p=y*p+.0208645
        p=y*p-.129013
        p=y*p+.85777
        t=y*p+1.0125
      else
        z=alog(y)-.775
        p=(.775-alog(z))/(1.+z)
        p=1./(1.+p)
        t=y*p/z
      end if
      return
      end

