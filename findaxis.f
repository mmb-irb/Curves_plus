      subroutine findaxis(A, eps, dta, m, gamma, x1, it)
c***************************************************************
c*  Find the axis of a rotation through an angle pi by finding *
c*  the eigenvector associated with the greatest eigenvalue of *
c*  the real symmetric matrix Q + 2I by the power method       *
c* ----------------------------------------------------------- *
c* INPUTS:                                                     *
c*         eps   : smallest number in double precision         *
c*         dta   : required precision                          *
c*         m     : maximum number of iterations                *
c*         A     : real square matrix A(3,3)                   *
c* OUTPUTS:                                                    *
c*         it    : error indicator: -1=no convergence,         *
c*                 0=method cannot be applied,                 *
c*                 1=convergence ok.                           *
c*         gamma : greatest eigenvalue (in absolute value)     *
c*                 of input matrix A(n,n)                      *
c*         X1    : associated eigenvector                      *
c***************************************************************
      integer i, it, m, n
      real*8 gamma
      real*8 A(3,3), X(3)
      real*8 eps, dta, X1(3)
      real*8 phi, s, X0(3)
 
      n = 3
      do i=1,n
      X0(i)=1.d0/DSQRT(DFLOAT(i))
      enddo
      it=-1
      l=1
      do while (it==-1.and.l<=m)
      gamma=0.d0
         do i=1,n
         X1(i)=0.d0
         do j=1,n
      X1(i)=X1(i)+A(i,j)*X0(j)
         enddo
         if(dabs(X1(i)).gt.dabs(gamma)) gamma=X1(i)
         enddo
      if(dabs(gamma).lt. eps) then
      it=0
      else
      do i=1,n
      X1(i)=X1(i)/gamma
      enddo
      phi=0.d0
         do i=1,n
	 s=dabs(X1(i)-X0(i))
	 if(s.gt.phi) phi=s
        enddo
      if(phi<dta) then
      it=1
      else
         do i=1,n
         X0(i)=X1(i)
         enddo
         l=l+1
         endif
      endif
      enddo
      return
      end
