      subroutine lslig(lb,dcor,rms,nind)
      include 'curves_data.inc'
      real*8 h(3,3),k(3,3)
      character*4 mnam,munit,bln,liga
      character*20 ilnam
      logical*2 nind
      dimension cg1(3),cg2(3),u(3,3),w(21),v(6,6),dcor(40,3),nind(40)
      common/lgd/ilnam(20),paxe(20),blef(20,40,3),parl(20,9),
     1 ilref(20),iluni(20,2),ilfrm(20,-3:20),bln(20,40),nlig,liga(20)
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      rms=0.d0
      sq2=sqrt(2.d0)
      do i=1,3
      cg1(i)=0.
      cg2(i)=0.
      do j=1,3
      u(i,j)=0.
      enddo
      enddo
c------------------------------------select atoms of real base & find cg
      kb=0
      do i=1,ilref(lb)
      if(nind(i)) then
      kb=kb+1
      do j=1,3
      cg1(j)=cg1(j)+dcor(i,j)
      cg2(j)=cg2(j)+blef(lb,i,j)
      enddo
      endif
      enddo
      do j=1,3
      cg1(j)=cg1(j)/kb
      cg2(j)=cg2(j)/kb
      enddo
c-------------------------------------------------------find ls rotation
      do l=1,ilref(lb)
      if(nind(l)) then
      do i=1,3
      do j=1,3
      u(i,j)=u(i,j)+(blef(lb,l,i)-cg2(j))*(dcor(l,j)-cg1(j))/kb
      enddo
      enddo
      endif
      enddo
      det=u(1,1)*(u(2,2)*u(3,3)-u(2,3)*u(3,2))
     1   -u(1,2)*(u(2,1)*u(3,3)-u(2,3)*u(3,1))
     1   +u(1,3)*(u(2,1)*u(3,2)-u(2,2)*u(3,1))
      if(abs(det).lt.1.d-9) return
      m=0
      do i=1,6
      do j=1,i
      m=m+1
      if(i.gt.3.and.j.le.3) then
      w(m)=u(j,i-3)
      else
      w(m)=0.
      endif
      enddo
      enddo
      call eigen(w,v,6)
      if(det.lt.0.0.and.abs(w(3)-w(6)).lt.1.e-6) return
      do i=1,3
      do j=1,3
      h(i,j)=sq2*v(i,j)
      k(i,j)=sq2*v(i+3,j)
      enddo
      enddo
      sn=h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))
     1  +h(2,3)*(h(3,1)*h(1,2)-h(1,1)*h(3,2))
     1  +h(3,3)*(h(1,1)*h(2,2)-h(2,1)*h(1,2))
      if(sn.lt.0.) then
      do i=1,3
      h(i,3)=-h(i,3)
      k(i,3)=-k(i,3)
      enddo
      endif
      do i=1,3
      do j=1,3
      u(i,j)=k(i,1)*h(j,1)+k(i,2)*h(j,2)+sign(1.d0,det)*k(i,3)*h(j,3)
      enddo
      enddo
c-------------------------------------------------------make ls rotation
      rms=0.
      do l=1,ilref(lb)
      if(nind(l)) then
      x0=blef(lb,l,1)-cg2(1)
      y0=blef(lb,l,2)-cg2(2)
      z0=blef(lb,l,3)-cg2(3)
      xnew=u(1,1)*x0+u(1,2)*y0+u(1,3)*z0+cg1(1)
      ynew=u(2,1)*x0+u(2,2)*y0+u(2,3)*z0+cg1(2)
      znew=u(3,1)*x0+u(3,2)*y0+u(3,3)*z0+cg1(3)
      dx=xnew-dcor(l,1)
      dy=ynew-dcor(l,2)
      dz=znew-dcor(l,3)
      rms=rms+(dx*dx+dy*dy+dz*dz)
         dcor(l,1)=xnew
         dcor(l,2)=ynew
         dcor(l,3)=znew
      endif
      enddo
      rms=sqrt(rms/kb)
      return
      end
