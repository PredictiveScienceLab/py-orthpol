c
c
      subroutine kern(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
c
c This routine generates the kernels in the Gauss quadrature remainder
c term, namely
c
c           K(k)(z)=rho(k)(z)/pi(k)(z), k=0,1,2,...,n,
c
c where  rho(k)  are the output quantities of the routine  knum, and 
c pi(k)  the (monic) orthogonal polynomials. The results are returned 
c in the array  ker  as ker(k)=K(k-1)(z), k=1,2,...,n+1. All the other 
c input and output parameters have the same meaning as in the routine 
c knum.
c
      complex z,ker,rold,p0,p,pm1
      dimension a(numax),b(numax),ker(*),rold(*)
c
c The arrays  ker,rold  are assumed to have dimension  n+1.
c
      call knum(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
      p0=(0.,0.)
      p=(1.,0.)
      do 10 k=1,n
        pm1=p0
        p0=p
        p=(z-a(k))*p0-b(k)*pm1
        ker(k+1)=ker(k+1)/p
   10 continue
      return
      end

