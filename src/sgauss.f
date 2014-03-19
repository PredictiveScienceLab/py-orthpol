c
c
      subroutine sgauss(n,alpha,beta,eps,zero,weight,ierr,e)
c
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k) and the weights as
c weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
c dlambda. The routine computes the nodes as eigenvalues, and the
c weights in term of the first component of the respective normalized
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature
c                      formula; type integer
c                alpha,beta - - arrays of dimension  n  to be filled
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: zero- array of dimension  n  containing the Gaussian
c                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
c                      ...,n
c                weight - array of dimension  n  containing the
c                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
c                ierr- an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is
c                      negative.
c
c The array  e  is needed for working space.
c
      dimension alpha(n),beta(n),zero(n),weight(n),e(n)
cf2py integer intent(hide),depend(alpha) :: n=len(alpha)
cf2py real intent(in) :: alpha
cf2py real intent(in),depend(n),dimension(n) :: beta
cf2py real optional, intent(in) :: eps=1e-6
cf2py real intent(out),depend(n),dimension(n) :: zero
cf2py real intent(out),depend(n),dimension(n) :: weight
cf2py integer intent(out) :: ierr
cf2py real intent(hide),depend(n),dimension(n) :: e
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.) then
        ierr=-2
        return
      end if
      weight(1)=beta(1)
      if (n.eq.1) return
      weight(1)=1.
      e(n)=0.
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.) then
          ierr=-2
          return
        end if
        e(k-1)=sqrt(beta(k))
        weight(k)=0.
  100 continue
      do 240 l=1,n
        j=0
c
c Look for a small subdiagonal element.
c
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(abs(e(m)).le.eps*(abs(zero(m))+abs(zero(m+1)))) goto 120
  110   continue
  120   p=zero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
c
c Form shift.
c
        g=(zero(l+1)-p)/(2.*e(l))
        r=sqrt(g*g+1.)
        g=zero(m)-p+e(l)/(g+sign(r,g))
        s=1.
        c=1.
        p=0.
        mml=m-l
c
c For i=m-1 step -1 until l do ...
c
        do 200 ii=1,mml
          i=m-ii
          f=s*e(i)
          b=c*e(i)
          if(abs(f).lt.abs(g)) goto 150
          c=g/f
          r=sqrt(c*c+1.)
          e(i+1)=f*r
          s=1./r
          c=c*s
          goto 160
  150     s=f/g
          r=sqrt(s*s+1.)
          e(i+1)=g*r
          c=1./r
          s=s*c
  160     g=zero(i+1)-p
          r=(zero(i)-g)*s +2.*c*b
          p=s*r
          zero(i+1)=g+p
          g=c*r-b
c
c Form first component of vector.
c
          f=weight(i+1)
          weight(i+1)=s*weight(i)+c*f
          weight(i)=c*weight(i)-s*f
  200   continue
        zero(l)=zero(l)-p
        e(l)=g
        e(m)=0.
        goto 105
  240 continue
c
c Order eigenvalues and eigenvectors.
c
      do 300 ii=2,n
        i=ii-1
        k=i
        p=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.p) goto 260
          k=j
          p=zero(j)
  260   continue
        if(k.eq.i) goto 300
        zero(k)=zero(i)
        zero(i)=p
        p=weight(i)
        weight(i)=weight(k)
        weight(k)=p
  300 continue
      do 310 k=1,n
        weight(k)=beta(1)*weight(k)*weight(k)
  310 continue
      return
c
c Set error - no convergence to an eigenvalue after 30 iterations.
c
  400 ierr=l
      return
      end

