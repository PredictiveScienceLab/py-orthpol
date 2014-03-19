c
c
      subroutine sradau(n,alpha,beta,end,zero,weight,ierr,e,a,b)
c
c Given  n  and a measure  dlambda, this routine generates the
c (n+1)-point Gauss-Radau quadrature formula
c
c   integral over supp(dlambda) of f(t)dlambda(t)
c
c     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
c =w(k), k=0,1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure
c dlambda. The nodes and weights are computed as eigenvalues and
c in terms of the first component of the respective normalized
c eigenvectors of a slightly modified Jacobi matrix of order  n+1.
c To do this, the routine calls upon the subroutine  gauss. It also
c uses the function subroutine  r1mach.
c
c    Input:  n - -  the number of interior points in the Gauss-Radau
c                   formula; type integer
c            alpha,beta - arrays of dimension  n+1  to be supplied with
c                   the recursion coefficients  alpha(k-1), beta(k-1),
c                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not
c                   used by the routine
c            end -  the prescribed endpoint  x(0)  of the Gauss-Radau
c                   formula; type real
c
c    Output: zero - array of dimension  n+1  containing the nodes (in
c                   increasing order)  zero(k)=x(k), k=0,1,2,...,n
c            weight-array of dimension  n+1  containing the weights
c                   weight(k)=w(k), k=0,1,2,...,n
c            ierr - an error flag inherited from the routine  gauss
c
c The arrays  e,a,b  are needed for working space.
c
      dimension alpha(*),beta(*),zero(*),weight(*),e(*),a(*),b(*)
c
c The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
c dimension  n+1.
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real intent(in) :: alpha
cf2py real intent(in),depend(n),dimension(n) :: beta
cf2py real intent(in) :: end
cf2py real intent(out),depend(n),dimension(n+1) :: zero
cf2py real intent(out),depend(n),dimension(n+1) :: weight
cf2py real intent(hide),depend(n),dimension(n+1) :: e
cf2py real intent(hide),depend(n),dimension(n+1) :: a
cf2py real intent(hide),depend(n),dimension(n+1) :: b
c
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      np1=n+1
      do 10 k=1,np1
        a(k)=alpha(k)
        b(k)=beta(k)
   10 continue
      p0=0.
      p1=1.
      do 20 k=1,n
        pm1=p0
        p0=p1
        p1=(end-a(k))*p0-b(k)*pm1
   20 continue
      a(np1)=end-b(np1)*p0/p1
      call sgauss(np1,a,b,epsma,zero,weight,ierr,e)
      return
      end

