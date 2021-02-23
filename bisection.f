      subroutine bisection(a,b,r)
      include 'curves_data.inc'
      integer*4 a,b,dtol,m,r,iter,iterMax
      real*8 dotdelta,dx,fa,fm
      parameter(iterMax=1000,dtol=1)
      fa=dotdelta(a)
      fm=dotdelta(b)
c begin bisection
      if(fa*fm.gt.0) return
      if(fa.lt.0) then
         r=a
         dx=b-a
      else
         r=b
         dx=a-b
      endif
c bisections loop
      do iter=1,iterMax
         dx=0.5d0*dx
         m =nint(dx+float(r))
         fm=dotdelta(m)
         if(fm.le.0) r=m
c this check seldom saves 1 iteration
         if(fm.eq.0) return
c decide which edge of step with root
         if(abs(dx).lt.dtol) then
            m =r+sign(1.d0,dx)*dtol
            fm=dotdelta(m)
            fa=dotdelta(r)
            if(abs(fm).lt.abs(fa)) r=m
            return
         endif
      enddo
      write(6,*) 'Too many iterations in routine "bisection"'
      end
