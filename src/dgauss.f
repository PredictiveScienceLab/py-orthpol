c
c
      subroutine dgauss(n,alpha,beta,deps,zero,weigh,ierr,de)
c
c This is a double-precision version of the routine  gauss.
c
      double precision alpha,beta,deps,zero,weigh,de,dp,dg,dr,
     *ds,dc,df,db
      dimension alpha(n),beta(n),zero(n),weigh(n),de(n)
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real*8 intent(in) :: alpha
cf2py real*8 intent(in),depend(n),dimension(n) :: beta
cf2py real*8 optional, intent(in) :: eps=1e-6
cf2py real*8 intent(out),depend(n),dimension(n) :: zero
cf2py real*8 intent(out),depend(n),dimension(n) :: weight
cf2py integer intent(out) :: ierr
cf2py real*8 intent(hide),depend(n),dimension(n) :: de
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.d0) then
        ierr=-2
        return
      end if
      weigh(1)=beta(1)
      if (n.eq.1) return
      weigh(1)=1.d0
      de(n)=0.d0
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.d0) then
          ierr=-2
          return
        end if
        de(k-1)=dsqrt(beta(k))
        weigh(k)=0.d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(dabs(de(m)).le.deps*(dabs(zero(m))+dabs(zero(m+1))))
     *      goto 120
  110   continue
  120   dp=zero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(zero(l+1)-dp)/(2.d0*de(l))
        dr=dsqrt(dg*dg+1.d0)
        dg=zero(m)-dp+de(l)/(dg+dsign(dr,dg))
        ds=1.d0
        dc=1.d0
        dp=0.d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*de(i)
          db=dc*de(i)
          if(dabs(df).lt.dabs(dg)) goto 150
          dc=dg/df
          dr=dsqrt(dc*dc+1.d0)
          de(i+1)=df*dr
          ds=1.d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=dsqrt(ds*ds+1.d0)
          de(i+1)=dg*dr
          dc=1.d0/dr
          ds=ds*dc
  160     dg=zero(i+1)-dp
          dr=(zero(i)-dg)*ds+2.d0*dc*db
          dp=ds*dr
          zero(i+1)=dg+dp
          dg=dc*dr-db
          df=weigh(i+1)
          weigh(i+1)=ds*weigh(i)+dc*df
          weigh(i)=dc*weigh(i)-ds*df
  200   continue
        zero(l)=zero(l)-dp
        de(l)=dg
        de(m)=0.d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.dp) goto 260
          k=j
          dp=zero(j)
  260   continue
        if(k.eq.i) goto 300
        zero(k)=zero(i)
        zero(i)=dp
        dp=weigh(i)
        weigh(i)=weigh(k)
        weigh(k)=dp
  300 continue
      do 310 k=1,n
        weigh(k)=beta(1)*weigh(k)*weigh(k)
  310 continue
      return
  400 ierr=l
      return
      end

