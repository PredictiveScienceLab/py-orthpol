c
c
      subroutine chri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)
c
c This subroutine implements the Christoffel or generalized Christoffel
c theorem. In all cases except  iopt=7, it uses nonlinear recurrence
c algorithms described in W. Gautschi,An algorithmic implementation
c of the generalized Christoffel theorem'', Numerical Integration
c (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
c iopt=7  incorporates a QR step with shift  x  in the manner of
c J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
c Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
c Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
c Problem'', Clarendon Press, Oxford, 1965. Given the recursion
c coefficients  a(k),b(k), k=0,1,...,n, for the (monic) orthogonal
c polynomials with respect to some measure  dlambda(t), it generates
c the recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
c measure
c
c              (t-x)dlambda(t)               if  iopt=1
c              [(t-x)**2+y**2]dlambda(t)     if  iopt=2
c              (t**2+y**2)dlambda(t) with    if  iopt=3
c                dlambda(t) and supp(dlambda)
c                symmetric  with respect to
c                the origin
c              dlambda(t)/(t-x)              if  iopt=4
c              dlambda(t)/[(t-x)**2+y**2]    if  iopt=5
c              dlambda(t)/(t**2+y**2) with   if  iopt=6
c                dlambda(t) and supp(dlambda)
c                symmetric with respect to
c                the origin
c              [(t-x)**2]dlambda(t)          if  iopt=7
c
c
c      Input:   n  - - - the number of recurrence coefficients
c                        desired; type integer
c               iopt - - an integer selecting the desired weight
c                        distribution
c               a,b  - - arrays of dimension  n+1  containing the
c                        recursion coefficients a(k-1),b(k-1),k=1,2,
c                        ...,n+1, of the polynomials orthogonal with
c                        respect to the given measure  dlambda(t)
c               x,y  - - real parameters defining the linear and 
c                        quadratic factors, or divisors, of  dlambda(t)
c               hr,hi  - the real and imaginary part, respectively, of
c                        the integral of dlambda(t)/(z-t), where z=x+iy;
c                        the parameter  hr  is used only if  iopt=4 or 
c                        5, the parameter  hi  only if  iopt=5 or 6
c
c      Output:  alpha,beta - - arrays of dimension  n  containing the
c                         desired recursion coefficients  alpha(k-1),
c                         beta(k-1), k=1,2,...,n
c
c It is assumed that  n  is larger than or equal to 2. Otherwise, the
c routine exits immediately with the error flag  ierr  set equal to 1.
c If  iopt  is not between 1 and 7, the routine exits with  ierr=2. 
c 
c The routine uses the function subroutine  r1mach  to evaluate the
c constant  eps , which is used only if  iopt=7.
c
      dimension a(*),b(*),alpha(n),beta(n)
c
c The arrays  a,b  are assumed to have dimension  n+1.
c
      eps=5.*r1mach(3)
c
c The quantity  eps  is a constant slightly larger than the machine
c precision.
c
      ierr=0
      if(n.lt.2) then
        ierr=1
        return
      end if
c
c What follows implements Eq. (3.7) of W. Gautschi, op. cit.
c
      if (iopt.eq.1) then
        e=0.
        do 10 k=1,n
          q=a(k)-e-x
          beta(k)=q*e
          e=b(k+1)/q
          alpha(k)=x+q+e
   10   continue
c
c Set the first beta-coefficient as discussed in Section 5.1 of the
c companion paper.
c
        beta(1)=b(1)*(a(1)-x)
        return
c
c What follows implements Eq. (4.7) of W. Gautschi, op. cit.
c
      else if(iopt.eq.2) then
        s=x-a(1)
        t=y
        eio=0.
        do 20 k=1,n
          d=s*s+t*t
          er=-b(k+1)*s/d
          ei=b(k+1)*t/d
          s=x+er-a(k+1)
          t=y+ei
          alpha(k)=x+t*er/ei-s*ei/t
          beta(k)=t*eio*(1.+(er/ei)**2)
          eio=ei
   20   continue
c
c Set the first beta-coefficient.
c
        beta(1)=b(1)*(b(2)+(a(1)-x)**2+y*y)
        return
c
c What follows implements Eq. (4.8) of W. Gautschi, op. cit.
c
      else if(iopt.eq.3) then
        t=y
        eio=0.
        do 30 k=1,n
          ei=b(k+1)/t
          t=y+ei
          alpha(k)=0.
          beta(k)=t*eio
          eio=ei
   30   continue
c
c Set the first beta-coefficient.
c
        beta(1)=b(1)*(b(2)+y*y)
        return
c
c What follows implements Eqs. (5.1),(5.2) of W. Gautschi, op. cit.
c
      else if(iopt.eq.4) then
        alpha(1)=x-b(1)/hr
        beta(1)=-hr
        q=-b(1)/hr
        do 40 k=2,n
          e=a(k-1)-x-q
          beta(k)=q*e
          q=b(k)/e
          alpha(k)=q+e+x
   40   continue
        return
c
c What follows implements Eq. (5.8) of W. Gautschi, op. cit.
c
      else if(iopt.eq.5) then
        nm1=n-1
        d=hr*hr+hi*hi
        eroo=a(1)-x+b(1)*hr/d
        eioo=-b(1)*hi/d-y
        alpha(1)=x+hr*y/hi
        beta(1)=-hi/y
        alpha(2)=x-b(1)*hi*eroo/(d*eioo)+hr*eioo/hi
        beta(2)=y*eioo*(1.+(hr/hi)**2)
        if(n.eq.2) return
        so=b(2)/(eroo**2+eioo**2)
        ero=a(2)-x-so*eroo
        eio=so*eioo-y
        alpha(3)=x+eroo*eio/eioo+so*eioo*ero/eio
        beta(3)=-b(1)*hi*eio*(1.+(eroo/eioo)**2)/d
        if(n.eq.3) return
        do 50 k=3,nm1
          s=b(k)/(ero**2+eio**2)
          er=a(k)-x-s*ero
          ei=s*eio-y
          alpha(k+1)=x+ero*ei/eio+s*eio*er/ei
          beta(k+1)=so*eioo*ei*(1.+(ero/eio)**2)
          eroo=ero
          eioo=eio
          ero=er
          eio=ei
          so=s
   50   continue
        return
c
c What follows implements Eq. (5.9) of W. Gautschi, op. cit.
c
      else if(iopt.eq.6) then
        nm1=n-1
        eoo=-b(1)/hi-y
        eo=b(2)/eoo-y
        alpha(1)=0.
        beta(1)=-hi/y
        alpha(2)=0.
        beta(2)=y*eoo
        if(n.eq.2) return
        alpha(3)=0.
        beta(3)=-b(1)*eo/hi
        if(n.eq.3) return
        do 60 k=3,nm1
          e=b(k)/eo-y
          beta(k+1)=b(k-1)*e/eoo
          alpha(k+1)=0.
          eoo=eo
          eo=e
   60   continue
        return
c
c What follows implements a QR step with shift  x.
c
      else if(iopt.eq.7) then
        u=0.
        c=1.
        c0=0.
        do 70 k=1,n
          gamma=a(k)-x-u
          cm1=c0
          c0=c
          if(abs(c0).gt.eps) then
            p2=(gamma**2)/c0
          else
            p2=cm1*b(k)
          end if
          if(k.gt.1) beta(k)=s*(p2+b(k+1))
          s=b(k+1)/(p2+b(k+1))
          c=p2/(p2+b(k+1))
          u=s*(gamma+a(k+1)-x)
          alpha(k)=gamma+u+x
   70   continue
        beta(1)=b(1)*(b(2)+(x-a(1))**2)
        return
      else
        ierr=2
        return
      end if
      end

