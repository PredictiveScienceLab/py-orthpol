c
c
      subroutine drecur(n,ipoly,dal,dbe,da,db,iderr)
c
c This is a double-precision version of the routine  recur.
c
      external dgamma
      double precision dal,dbe,da,db,dlmach,d1mach,dkm1,dalpbe,dt,
     *dlga,dal2,dbe2,dgamma
      dimension da(n),db(n)
cf2py integer intent(in) :: n
cf2py integer intent(in) :: ipoly
cf2py real optional,intent(in) :: dal=-.5
cf2py real optional,intent(in) :: dbe=.5
cf2py real intent(out,out=a),depend(n),dimension(n) :: da
cf2py real intent(out,out=b),depend(n),dimension(n) :: db
cf2py integer intent(out,out=ierr) :: iderr
      if(n.lt.1) then
        iderr=3
        return
      end if
      dlmach=dlog(d1mach(2))
      iderr=0
      do 10 k=1,n
        da(k)=0.d0
   10 continue
      if(ipoly.eq.1) then
        db(1)=2.d0
        if (n.eq.1) return
        do 20 k=2,n
          dkm1=dble(k-1)
          db(k)=1.d0/(4.d0-1.d0/(dkm1*dkm1))
   20   continue
        return
      else if(ipoly.eq.2) then
        da(1)=.5d0
        db(1)=1.d0
        if(n.eq.1) return
        do 30 k=2,n
          da(k)=.5d0
          dkm1=dble(k-1)
          db(k)=.25d0/(4.d0-1.d0/(dkm1*dkm1))
   30   continue
        return
      else if(ipoly.eq.3) then
        db(1)=4.d0*datan(1.d0)
        if(n.eq.1) return
        db(2)=.5d0
        if(n.eq.2) return
        do 40 k=3,n
          db(k)=.25d0
   40   continue
        return
      else if(ipoly.eq.4) then
        db(1)=2.d0*datan(1.d0)
        if(n.eq.1) return
        do 50 k=2,n
          db(k)=.25d0
   50   continue
        return
      else if(ipoly.eq.5) then
        db(1)=4.d0*datan(1.d0)
        da(1)=.5d0
        if(n.eq.1) return
        do 60 k=2,n
          db(k)=.25d0
   60   continue
        return
      else if(ipoly.eq.6) then
        if(dal.le.-1.d0 .or. dbe.le.-1.d0) then
          iderr=1
          return
        else
          dalpbe=dal+dbe
          da(1)=(dbe-dal)/(dalpbe+2.d0)
          dt=(dalpbe+1.d0)*dlog(2.d0)+dlga(dal+1.d0)+dlga(dbe+1.d0)-
     *      dlga(dalpbe+2.d0)
          if(dt.gt.dlmach) then
            iderr=2
            db(1)=d1mach(2)
          else
            db(1)=dexp(dt)
          end if
          if(n.eq.1) return
          dal2=dal*dal
          dbe2=dbe*dbe
          da(2)=(dbe2-dal2)/((dalpbe+2.d0)*(dalpbe+4.d0))
          db(2)=4.d0*(dal+1.d0)*(dbe+1.d0)/((dalpbe+3.d0)*(dalpbe+
     *      2.d0)**2)
          if(n.eq.2) return
          do 70 k=3,n
            dkm1=dble(k-1)
            da(k)=.25d0*(dbe2-dal2)/(dkm1*dkm1*(1.d0+.5d0*dalpbe/dkm1)
     *        *(1.d0+.5d0*(dalpbe+2.d0)/dkm1))
            db(k)=.25d0*(1.d0+dal/dkm1)*(1.d0+dbe/dkm1)*(1.d0+dalpbe/
     *        dkm1)/((1.d0+.5d0*(dalpbe+1.d0)/dkm1)*(1.d0+.5d0*(dalpbe
     *      -1.d0)/dkm1)*(1.d0+.5d0*dalpbe/dkm1)**2)
   70     continue
          return
        end if
      else if(ipoly.eq.7) then
        if(dal.le.-1.d0) then
          iderr=1
          return
        else
          da(1)=dal+1.d0
          db(1)=dgamma(dal+1.d0,iderr)
          if(iderr.eq.2) db(1)=d1mach(2)
          if(n.eq.1) return
          do 80 k=2,n
            dkm1=dble(k-1)
            da(k)=2.d0*dkm1+dal+1.d0
            db(k)=dkm1*(dkm1+dal)
   80     continue
          return
        end if
      else if(ipoly.eq.8) then
        db(1)=dsqrt(4.d0*datan(1.d0))
        if(n.eq.1) return
        do 90 k=2,n
          db(k)=.5d0*dble(k-1)
   90   continue
        return
      else
        iderr=4
      end if
      end

      double precision function dlga(dx)
      double precision dbnum,dbden,dx,d1mach,dc,dp,dy,dt,ds
      dimension dbnum(8),dbden(8)
c
c This routine evaluates the logarithm of the gamma function by a
c combination of recurrence and asymptotic approximation.
c
c The entries in the next data statement are the numerators and
c denominators, respectively, of the quantities B[16]/(16*15),
c B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
c numbers.
c
      data dbnum/-3.617d3,1.d0,-6.91d2,1.d0,-1.d0,1.d0,-1.d0,1.d0/,
     *     dbden/1.224d5,1.56d2,3.6036d5,1.188d3,1.68d3,1.26d3,3.6d2,
     *1.2d1/
c
c The quantity  dprec  in the next statement is the number of decimal
c digits carried in double-precision floating-point arithmetic.
c
      dprec=-alog10(sngl(d1mach(3)))
      dc=.5d0*dlog(8.d0*datan(1.d0))
      dp=1.d0
      dy=dx
      y=sngl(dy)
c
c The quantity  y0  below is the threshold value beyond which asymptotic
c evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
c and I.A. Stegun,Handbook of Mathematical Functions''. The constants
c are .12118868... = ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
c
      y0=exp(.121189*dprec+.053905)
   10 if(y.gt.y0) goto 20
      dp=dy*dp
      dy=dy+1.d0
      y=sngl(dy)
      goto 10
   20 dt=1.d0/(dy*dy)
c
c The right-hand side of the next assignment statement is B[18]/(18*17).
c
      ds=4.3867d4/2.44188d5
      do 30 i=1,8
        ds=dt*ds+dbnum(i)/dbden(i)
   30 continue
      dlga=(dy-.5d0)*dlog(dy)-dy+dc+ds/dy-dlog(dp)
      return
      end

      double precision function dgamma(dx,iderr)
c
c This evaluates the gamma function for real positive  dx, using the
c function subroutine  dlga.
c
      double precision dx,dlmach,d1mach,dt,dlga
      dlmach=dlog(d1mach(2))
      iderr=0
      dt=dlga(dx)
      if(dt.ge.dlmach) then
        iderr=2
        dgamma=d1mach(2)
        return
      else
        dgamma=dexp(dt)
        return
      end if
      end

