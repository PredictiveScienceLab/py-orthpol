c
c
      subroutine dradau(n,alpha,beta,end,zero,weight,ierr,de,
     *da,db)
c
c This is a double-precision version of the routine  radau.
c
      double precision end,depsma,dp0,dp1,dpm1,alpha(*),beta(*),
     *zero(*),weight(*),de(*),da(*),db(*),d1mach
c
c The arrays  alpha,beta,zero,weight,de,da,db  are assumed to have
c dimension  n+1.
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real intent(in) :: alpha
cf2py real intent(in),depend(n),dimension(n) :: beta
cf2py real intent(in) :: end
cf2py real intent(out),depend(n),dimension(n+1) :: zero
cf2py real intent(out),depend(n),dimension(n+1) :: weight
cf2py real intent(hide),depend(n),dimension(n+1) :: de
cf2py real intent(hide),depend(n),dimension(n+1) :: da
cf2py real intent(hide),depend(n),dimension(n+1) :: db
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      do 10 k=1,np1
        da(k)=alpha(k)
        db(k)=beta(k)
   10 continue
      dp0=0.d0
      dp1=1.d0
      do 20 k=1,n
        dpm1=dp0
        dp0=dp1
        dp1=(end-da(k))*dp0-db(k)*dpm1
   20 continue
      da(np1)=end-db(np1)*dp0/dp1
      call dgauss(np1,da,db,depsma,zero,weight,ierr,de)
      return
      end

