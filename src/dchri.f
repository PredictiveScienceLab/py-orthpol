c
c
      subroutine dchri(n,iopt,da,db,dx,dy,dhr,dhi,dalpha,dbeta,ierr)
c
c This is a double-precision version of the routine  chri.
c
      double precision da,db,dx,dy,dhr,dhi,dalpha,dbeta,deps,d1mach,
     *de,dq,ds,dt,deio,dd,der,dei,deroo,deioo,dso,dero,deoo,deo,du,dc,
     *dc0,dgam,dcm1,dp2
      dimension da(*),db(*),dalpha(n),dbeta(n)
c
c The arrays  da,db  are assumed to have dimension  n+1.
c
      deps=5.d0*d1mach(3)
      ierr=0
      if(n.lt.2) then
        ierr=1
        return
      end if
      if(iopt.eq.1) then
        de=0.d0
        do 10 k=1,n
          dq=da(k)-de-dx
          dbeta(k)=dq*de
          de=db(k+1)/dq
          dalpha(k)=dx+dq+de
   10   continue
        dbeta(1)=db(1)*(da(1)-dx)
        return
      else if(iopt.eq.2) then
        ds=dx-da(1)
        dt=dy
        deio=0.d0
        do 20 k=1,n
          dd=ds*ds+dt*dt
          der=-db(k+1)*ds/dd
          dei=db(k+1)*dt/dd
          ds=dx+der-da(k+1)
          dt=dy+dei
          dalpha(k)=dx+dt*der/dei-ds*dei/dt
          dbeta(k)=dt*deio*(1.d0+(der/dei)**2)
          deio=dei
   20   continue
        dbeta(1)=db(1)*(db(2)+(da(1)-dx)**2+dy*dy)
        return
      else if(iopt.eq.3) then
        dt=dy
        deio=0.d0
        do 30 k=1,n
          dei=db(k+1)/dt
          dt=dy+dei
          dalpha(k)=0.d0
          dbeta(k)=dt*deio
          deio=dei
   30   continue
        dbeta(1)=db(1)*(db(2)+dy*dy)
        return
      else if(iopt.eq.4) then
        dalpha(1)=dx-db(1)/dhr
        dbeta(1)=-dhr
        dq=-db(1)/dhr
        do 40 k=2,n
          de=da(k-1)-dx-dq
          dbeta(k)=dq*de
          dq=db(k)/de
          dalpha(k)=dq+de+dx
   40   continue
        return
      else if(iopt.eq.5) then
        nm1=n-1
        dd=dhr*dhr+dhi*dhi
        deroo=da(1)-dx+db(1)*dhr/dd
        deioo=-db(1)*dhi/dd-dy
        dalpha(1)=dx+dhr*dy/dhi
        dbeta(1)=-dhi/dy
        dalpha(2)=dx-db(1)*dhi*deroo/(dd*deioo)+dhr*deioo/dhi
        dbeta(2)=dy*deioo*(1.d0+(dhr/dhi)**2)
        if(n.eq.2) return
        dso=db(2)/(deroo**2+deioo**2)
        dero=da(2)-dx-dso*deroo
        deio=dso*deioo-dy
        dalpha(3)=dx+deroo*deio/deioo+dso*deioo*dero/deio
        dbeta(3)=-db(1)*dhi*deio*(1.d0+(deroo/deioo)**2)/dd
        if(n.eq.3) return
        do 50 k=3,nm1
          ds=db(k)/(dero**2+deio**2)
          der=da(k)-dx-ds*dero
          dei=ds*deio-dy
          dalpha(k+1)=dx+dero*dei/deio+ds*deio*der/dei
          dbeta(k+1)=dso*deioo*dei*(1.d0+(dero/deio)**2)
          deroo=dero
          deioo=deio
          dero=der
          deio=dei
          dso=ds
   50   continue
        return
      else if(iopt.eq.6) then
        nm1=n-1
        deoo=-db(1)/dhi-dy
        deo=db(2)/deoo-dy
        dalpha(1)=0.d0
        dbeta(1)=-dhi/dy
        dalpha(2)=0.d0
        dbeta(2)=dy*deoo
        if(n.eq.2) return
        dalpha(3)=0.d0
        dbeta(3)=-db(1)*deo/dhi
        if(n.eq.3) return
        do 60 k=3,nm1
          de=db(k)/deo-dy
          dbeta(k+1)=db(k-1)*de/deoo
          dalpha(k+1)=0.d0
          deoo=deo
          deo=de
   60   continue
        return
      else if(iopt.eq.7) then
        du=0.d0
        dc=1.d0
        dc0=0.d0
        do 70 k=1,n
          dgam=da(k)-dx-du
          dcm1=dc0
          dc0=dc
          if(dabs(dc0).gt.deps) then
            dp2=(dgam**2)/dc0
          else
            dp2=dcm1*db(k)
          end if
          if(k.gt.1) dbeta(k)=ds*(dp2+db(k+1))
          ds=db(k+1)/(dp2+db(k+1))
          dc=dp2/(dp2+db(k+1))
          du=ds*(dgam+da(k+1)-dx)
          dalpha(k)=dgam+du+dx
   70   continue
        dbeta(1)=db(1)*(db(2)+(dx-da(1))**2)
        return
      else
        ierr=2
        return
      end if
      end

