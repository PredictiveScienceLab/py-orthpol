c
c
      subroutine dknum(n,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *ierr,droldr,droldi)
c
c This is a double-precision version of the routine  knum.
c
      double precision dx,dy,deps,da(numax),db(numax),drhor(*),
     *drhoi(*),droldr(*),droldi(*),drr,dri,dden,dt
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
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
        drhor(k)=0.d0
        drhoi(k)=0.d0
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
        droldr(k)=drhor(k)
        droldi(k)=drhoi(k)
   30 continue
      drr=0.d0
      dri=0.d0
      do 40 j=1,nu
        j1=nu-j+1
        dden=(dx-da(j1)-drr)**2+(dy-dri)**2
        drr=db(j1)*(dx-da(j1)-drr)/dden
        dri=-db(j1)*(dy-dri)/dden
        if(j1.le.np1) then
          drhor(j1)=drr
          drhoi(j1)=dri
        end if
   40 continue
      do 50 k=1,np1
        if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
     *    deps*(drhor(k)**2+drhoi(k)**2)) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
        dt=drhor(k)*drhor(k-1)-drhoi(k)*drhoi(k-1)
        drhoi(k)=drhor(k)*drhoi(k-1)+drhoi(k)*drhor(k-1)
        drhor(k)=dt
   70 continue
      return
      end 

