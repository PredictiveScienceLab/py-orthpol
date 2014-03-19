c
c

      subroutine dlob(n,alpha,beta,aleft,right,zero,weight,
     *ierr,de,da,db)
c
c This is a double-precision version of the routine  lob.
c
      double precision aleft,right,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l,
     *dpm1r,ddet,alpha(*),beta(*),zero(*),weight(*),de(*),da(*),
     *db(*),d1mach
c
c The arrays  alpha,beta,zero,weight,de,da,db  are assumed to have
c dimension  n+2.
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real intent(in) :: alpha
cf2py real intent(in),depend(n),dimension(n) :: beta
cf2py real intent(in) :: aleft
cf2py real intent(in) :: right
cf2py real intent(out),depend(n),dimension(n+2) :: zero
cf2py real intent(out),depend(n),dimension(n+2) :: weight
cf2py real intent(hide),depend(n),dimension(n+2) :: de
cf2py real intent(hide),depend(n),dimension(n+2) :: da
cf2py real intent(hide),depend(n),dimension(n+2) :: db
cf2py integer intent(out) :: ierr
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        da(k)=alpha(k)
        db(k)=beta(k)
   10 continue
      dp0l=0.d0
      dp0r=0.d0
      dp1l=1.d0
      dp1r=1.d0
      do 20 k=1,np1
        dpm1l=dp0l
        dp0l=dp1l
        dpm1r=dp0r
        dp0r=dp1r
        dp1l=(aleft-da(k))*dp0l-db(k)*dpm1l
        dp1r=(right-da(k))*dp0r-db(k)*dpm1r
   20 continue
      ddet=dp1l*dp0r-dp1r*dp0l
      da(np2)=(aleft*dp1l*dp0r-right*dp1r*dp0l)/ddet
      db(np2)=(right-aleft)*dp1l*dp1r/ddet
      call dgauss(np2,da,db,depsma,zero,weight,ierr,de)
      return
      end

