c
c
      subroutine knum(n,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
c
c This routine generates
c
c   rho(k)(z)=integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,
c
c where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
c with respect to the measure  dlambda(t), and the integral is extended
c over the support of  dlambda. It is assumed that  z  is a complex
c number outside the smallest interval containing the support of
c dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
c members of the minimal solution of the basic three-term recurrence
c relation
c
c      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,
c
c satisfied by the orthogonal polynomials  pi(k)(z).
c
c   Input:  n  - -  the largest integer  k  for which  rho(k)  is 
c                   desired
c           nu0  -  an estimate of the starting backward recurrence 
c                   index; if no better estimate is known, set 
c                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
c                   weight functions, estimates of  nu0  are generated
c                   respectively by the routines  nu0jac,nu0lag  and
c                   nu0her
c           numax - an integer larger than  n  cutting off backward 
c                   recursion in case of nonconvergence; if  nu0  
c                   exceeds  numax, then the routine aborts with the 
c                   error flag  ierr  set equal to  nu0
c           z - - - the variable in  rho(k)(z); type complex
c           eps - - the relative accuracy to which the  rho(k)  are
c                   desired
c           a,b - - arrays of dimension  numax  to be supplied with the
c                   recurrence coefficients  a(k-1), b(k-1), k=1,2,...,
c                   numax.
c
c   Output: rho - - an array of dimension  n+1  containing the results
c                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
c           nu  - - the starting backward recurrence index that yields
c                   convergence
c           ierr  - an error flag equal to zero on normal return, equal 
c                   to  nu0  if  nu0 > numax, and equal to  numax in 
c                   case of nonconvergence.
c
c The complex array  rold  of dimension  n+1  is used for working space.
c            
      complex z,rho,rold,r
      dimension a(numax),b(numax),rho(*),rold(*)
c
c The arrays  rho,rold  are assumed to have dimension  n+1.
c
      ierr=0
      np1=n+1
      if(nu0.gt.numax) then
        ierr=nu0
        return
      end if
      if(nu0.lt.np1) nu0=np1
      nu=nu0-5
      do 10 k=1,np1
        rho(k)=(0.,0.)
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
        rold(k)=rho(k)
   30 continue
      r=(0.,0.)
      do 40 j=1,nu
        j1=nu-j+1
        r=cmplx(b(j1),0.)/(z-cmplx(a(j1),0.)-r)
        if(j1.le.np1) rho(j1)=r
   40 continue
      do 50 k=1,np1
        if(cabs(rho(k)-rold(k)).gt.eps*cabs(rho(k))) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
        rho(k)=rho(k)*rho(k-1)
   70 continue
      return
      end 

