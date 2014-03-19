c
c
      subroutine ssti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)
c
c This routine applies Stieltjes's procedure'' (cf. Section 2.1 of
c W. Gautschi,On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow)
c or  +k  (in the case of overflow), where  k  is the recursion index
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close
c to  ncap, since there is a distinct possibility of numerical
c instability developing. (See W. Gautschi,Is the recurrence relation
c for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
c In that case, the routine  lancz  should be used.
c
c The routine uses the function subroutine  r1mach.
c
      dimension x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap),
     *p2(ncap)
cf2py integer intent(in) :: n
cf2py integer intent(hide),depend(x) :: ncap=len(x)
cf2py real intent(in) :: x
cf2py real intent(in),depend(ncap),check(len(w)>=ncap) :: w
cf2py real intent(out),depend(n),dimension(n) :: alpha
cf2py real intent(out),depend(n),dimension(n) :: beta
cf2py real intent(hide),depend(ncap),dimension(ncap) :: p0
cf2py real intent(hide),depend(ncap),dimension(ncap) :: p1
cf2py real intent(hide),depend(ncap),dimension(ncap) :: p2
cf2py integer intent(out) :: ierr
      tiny=10.*r1mach(1)
      huge=.1*r1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
c
c Compute the first alpha- and beta-coefficient.
c
      sum0=0.
      sum1=0.
      do 10 m=1,ncap
        sum0=sum0+w(m)
        sum1=sum1+w(m)*x(m)
   10 continue
      alpha(1)=sum1/sum0
      beta(1)=sum0
      if(n.eq.1) return
c
c Compute the remaining alpha- and beta-coefficients.
c
      do 20 m=1,ncap
        p1(m)=0.
        p2(m)=1.
   20 continue
      do 40 k=1,nm1
        sum1=0.
        sum2=0.
        do 30 m=1,ncap
c
c The following statement is designed to avoid an overflow condition
c in the computation of  p2(m)  when the weights  w(m)  go to zero
c faster (and underflow) than the  p2(m)  grow.
c
          if(w(m).eq.0.) goto 30
          p0(m)=p1(m)
          p1(m)=p2(m)
          p2(m)=(x(m)-alpha(k))*p1(m)-beta(k)*p0(m)
c
c Check for impending overflow.
c
          if(abs(p2(m)).gt.huge .or. abs(sum2).gt.huge) then
            ierr=k
            return
          end if
          t=w(m)*p2(m)*p2(m)
          sum1=sum1+t
          sum2=sum2+t*x(m)
   30   continue
c
c Check for impending underflow.
c
        if(abs(sum1).lt.tiny) then
          ierr=-k
          return
        end if
        alpha(k+1)=sum2/sum1
        beta(k+1)=sum1/sum0
        sum0=sum1
   40 continue
      return
      end

