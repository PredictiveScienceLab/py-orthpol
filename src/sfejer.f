c
c
      subroutine sfejer(n,x,w)
c
c This routine generates the n-point Fejer quadrature rule.
c
c         input:   n   - the number of quadrature nodes
c         output:  x,w - arrays of dimension  n  holding the quadrature
c                        nodes and weights, respectively; the nodes
c                        are ordered increasingly
c
      dimension x(n),w(n)
cf2py integer intent(in) :: n
cf2py real intent(out),depend(n),dimension(n) :: x
cf2py real intent(out),depend(n),dimension(n) :: w
      pi=4.*atan(1.)
      nh=n/2
      np1h=(n+1)/2
      fn=real(n)
      do 10 k=1,nh
        x(n+1-k)=cos(.5*real(2*k-1)*pi/fn)
        x(k)=-x(n+1-k)
   10 continue
      if(2*nh.ne.n) x(np1h)=0.
      do 30 k=1,np1h
        c1=1.
        c0=2.*x(k)*x(k)-1.
        t=2.*c0
        sum=c0/3.
        do 20 m=2,nh
          c2=c1
          c1=c0
          c0=t*c1-c2
          sum=sum+c0/real(4*m*m-1)
   20   continue
        w(k)=2.*(1.-2.*sum)/fn
        w(n+1-k)=w(k)
   30 continue
      return
      end

