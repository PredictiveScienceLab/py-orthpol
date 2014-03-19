c
c
      subroutine slob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)
c
c Given  n  and a measure  dlambda, this routine generates the
c (n+2)-point Gauss-Lobatto quadrature formula
c
c   integral over supp(dlambda) of f(x)dlambda(x)
c
c      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))
c
c              + w(n+1)f(x(n+1)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
c =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
c dlambda. The nodes and weights are computed in terms of the
c eigenvalues and first component of the normalized eigenvectors of
c a slightly modified Jacobi matrix of order  n+2. The routine calls
c upon the subroutine  gauss  and the function subroutine  r1mach.
c
c   Input:  n - -  the number of interior points in the Gauss-Lobatto
c                  formula; type integer
c           alpha,beta - arrays of dimension  n+2  to be supplied with
c                  the recursion coefficients  alpha(k-1), beta(k-1),
c                  k=1,2,...,n+2, of the underlying measure; the
c                  routine does not use  alpha(n+2), beta(n+2)
c           aleft,right - the prescribed left and right endpoints
c                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula
c
c   Output: zero - an array of dimension  n+2  containing the nodes (in
c                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
c           weight-an array of dimension  n+2  containing the weights
c                  weight(k)=w(k), k=0,1,...,n,n+1
c           ierr - an error flag inherited from the routine  gauss
c
c The arrays  e,a,b  are needed for working space.
c
      dimension alpha(*),beta(*),zero(*),weight(*),e(*),a(*),b(*)
c
c The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
c dimension  n+2.
c
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real intent(in) :: alpha
cf2py real intent(in),depend(n),dimension(n) :: beta
cf2py real intent(in) :: aleft
cf2py real intent(in) :: right
cf2py real intent(out),depend(n),dimension(n+2) :: zero
cf2py real intent(out),depend(n),dimension(n+2) :: weight
cf2py real intent(hide),depend(n),dimension(n+2) :: e
cf2py real intent(hide),depend(n),dimension(n+2) :: a
cf2py real intent(hide),depend(n),dimension(n+2) :: b
cf2py integer intent(out) :: ierr
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        a(k)=alpha(k)
        b(k)=beta(k)
   10 continue
      p0l=0.
      p0r=0.
      p1l=1.
      p1r=1.
      do 20 k=1,np1
        pm1l=p0l
        p0l=p1l
        pm1r=p0r
        p0r=p1r
        p1l=(aleft-a(k))*p0l-b(k)*pm1l
        p1r=(right-a(k))*p0r-b(k)*pm1r
   20 continue
      det=p1l*p0r-p1r*p0l
      a(np2)=(aleft*p1l*p0r-right*p1r*p0l)/det
      b(np2)=(right-aleft)*p1l*p1r/det
      call sgauss(np2,a,b,epsma,zero,weight,ierr,e)
      return
      end

