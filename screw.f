      subroutine screw(r1,r2,key)
      include 'curves_data.inc'
      character*4 name
      character*1 na,nt
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      dimension rsc(3),r1(4,3),r2(4,3),t(3,3),q(3,3),w(3),v(3),
     & scw(3),a(3,3)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
cjm---------------------------------------------------------------------
cjm  this file commented version of file original_screw.f that was
cjm  written by Maher Moakher, with comments starting with cjm added by
cjm  John Maddocks
c-----------------------------------------------------------------------
cjm  r1 and r2 4x3 arrays first three lines orthonormal vectors of frame
cjm  and fourth line coordinates of frame origin
cjm---------------------------------------------------------------------
cjm  first compute absolute coordinates of translation vector in v(j)
 
      do j=1,3
      v(j)=r2(4,j)-r1(4,j)
      enddo
      dist=sqrt(v(1)**2+v(2)**2+v(3)**2)
c-------------------------------------------------------------screw axis
cjm 1) compute matrix Q = R_1^T R_2 with R_2 =  R_1 Q, absolute coordinates
cjm here R_1 and R_2 contain the first 3 lines of r1 and r2
 
cjm---------------------------------------------------------------------
cjm   note rows of r1(3,3) block are columns d-i of usual rotation matrix
      do i=1,3
      do j=1,3
      q(i,j)=0.d0
        do k=1,3
        q(i,j)=q(i,j)+r1(k,i)*r2(k,j)
        enddo
      enddo
      enddo
cjm---------------------------------------------------------------------
cjm 2) use TrQ = 1 + 2 cos theta  to get theta in [0,pi]
cjm---------------------------------------------------------------------
      ct = (q(1,1)+q(2,2)+q(3,3)-1.d0)/2.d0
      if(abs(ct).gt.1.)ct=sign(1.d0,ct)
      theta = dacos(ct)
cjm---------------------------------------------------------------------
cjm 3) deal with special case theta close to zero and Q close to identity
cjm---------------------------------------------------------------------
      if (theta .lt. 1.d-6) then
           if(dist.lt.1.d-3) then
           do j=1,6
           var(j)=0.
           enddo
           return
           endif
 
        usc(1)= v(1)
        usc(2)= v(2)
        usc(3)= v(3)
 
cjm   special case angle close to pi
      else if (pi-theta.lt.1.d-6) then
        do i=1,3
           do j=1,3
              a(i,j) = q(i,j)
           enddo
           a(i,i) = a(i,i) + 1.d0
        enddo
        call findaxis(a, 1.d-10, 1.d-6, 40, gamma, usc, it)
cjm   this is the main case, usc is screw vector
      else
        usc(1)=-q(3,2)+q(2,3)
        usc(2)=-q(1,3)+q(3,1)
        usc(3)=-q(2,1)+q(1,2)
      endif
      unorm=dsqrt(usc(1)**2+usc(2)**2+usc(3)**2)
      do j=1,3
 
      usc(j)=usc(j)/unorm
      scw(j)=usc(j)*theta*crd
      enddo
 
      hdot=usc(1)*v(1)+usc(2)*v(2)+usc(3)*v(3)
cjm   after here usc(1,2,3) holds a unit vector with usc.v ge 0,
cjm   and theta in -pi,pi and axial vector is 2 sin theta usc
 
cjm   compute three junction angle variables wrt R1 frame
        do k=1,3
        dot=0.d0
        do j=1,3
        dot=dot+scw(j)*r1(k,j)
        enddo
        var(k+3)=dot
        enddo
 
c------------------------------------------------------------------point
cjm   w is half of projection of vector v perpendicular to axis vector
 
      w(1)=(v(1)-hdot*usc(1))/2.d0
      w(2)=(v(2)-hdot*usc(2))/2.d0
      w(3)=(v(3)-hdot*usc(3))/2.d0
 
cjm   now compute centre of circle using (e1+e2)/2 plus (usc x w)/tan(theta/2)
      tn=1/tan(theta/2.d0)
      usc(4)=(r2(4,1)+r1(4,1))/2.d0+tn*(-w(2)*usc(3)+w(3)*usc(2))
      usc(5)=(r2(4,2)+r1(4,2))/2.d0+tn*(-w(3)*usc(1)+w(1)*usc(3))
      usc(6)=(r2(4,3)+r1(4,3))/2.d0+tn*(-w(1)*usc(2)+w(2)*usc(1))
 
c--------------------------------------------------------------transform
      if(key.lt.0) return
      st=sin(theta/2)
      ct=cos(theta/2)
      cm=ct-1
      urot(4,1)=(r1(4,1)+r2(4,1))/2
      urot(4,2)=(r1(4,2)+r2(4,2))/2
      urot(4,3)=(r1(4,3)+r2(4,3))/2
 
cjm     compute halfway frame using rodrigues formula and half angle
        t(1,1)= ct       -cm*usc(1)*usc(1)
        t(1,2)=-st*usc(3)-cm*usc(1)*usc(2)
        t(1,3)= st*usc(2)-cm*usc(1)*usc(3)
        t(2,1)= st*usc(3)-cm*usc(2)*usc(1)
        t(2,2)= ct       -cm*usc(2)*usc(2)
        t(2,3)=-st*usc(1)-cm*usc(2)*usc(3)
        t(3,1)=-st*usc(2)-cm*usc(3)*usc(1)
        t(3,2)= st*usc(1)-cm*usc(3)*usc(2)
        t(3,3)= ct       -cm*usc(3)*usc(3)
      do i=1,3
      dot=0.d0
      do j=1,3
        urot(i,j)=0.d0
        do k=1,3
        urot(i,j)=urot(i,j)+r1(i,k)*t(j,k)
        enddo
      dot=dot+v(j)*urot(i,j)
      enddo
      var(i)=dot
      enddo
      return
      end
