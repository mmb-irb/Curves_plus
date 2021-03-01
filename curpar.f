      subroutine curpar
      include 'curves_data.inc'
      character*1 na,nt
      character*4 snam,sunit
      character*128 file,ftop,lis,lib,ibld,sol,back
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 lpa,traj
      parameter (cfac=40.0d0)
      dimension r1(4,3),r2(4,3)
      common/cha/file,ftop,lis,lib,ibld,sol,back
      common/dat/wback,wbase,rvfac,isym,itst,itnd,itdel,itbkt,
     1 naxlim,circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
      common/sto/bpdat(n3,n5),traj
c-------------------------------------------------------single point O/P
      if(.not.traj) then
      kas=0
      khs=0
      write(6,10)
10    format(/2x,'(F) Curvature analysis',
     1      //5x,'N',4x,'BP step',8x,'Cur',6x,'Rad',6x,'Reg'
     1 /2x,'---------------------------------------------')
      endif
c----------------------------------------------------------for each step
      risum=0.d0
      cusum=0.d0
      cusum2=0.d0
      if(circ) then
      ido=1
      iup=nlev
         else
         ido=2
         iup=nlev-2
         endif
      do i=ido,iup
      i1=i-1
      i2=i
      i3=i+1
      i4=i+2
      if(i1.eq.0) i1=nlev
      if(i3.gt.nlev) i3=i3-nlev
      if(i4.gt.nlev) i4=i4-nlev
c-------------------------------fit circles to 3 pts and calculate radii
      ax=uvw(i2,4,1)-uvw(i1,4,1)
      ay=uvw(i2,4,2)-uvw(i1,4,2)
      az=uvw(i2,4,3)-uvw(i1,4,3)
      a2=ax*ax+ay*ay+az*az
      a=sqrt(a2)
      bx=uvw(i3,4,1)-uvw(i2,4,1)
      by=uvw(i3,4,2)-uvw(i2,4,2)
      bz=uvw(i3,4,3)-uvw(i2,4,3)
      b2=bx*bx+by*by+bz*bz
      b=sqrt(b2)
      cx=uvw(i4,4,1)-uvw(i3,4,1)
      cy=uvw(i4,4,2)-uvw(i3,4,2)
      cz=uvw(i4,4,3)-uvw(i3,4,3)
      c2=cx*cx+cy*cy+cz*cz
      c=sqrt(c2)
         dx=uvw(i3,4,1)-uvw(i1,4,1)
         dy=uvw(i3,4,2)-uvw(i1,4,2)
         dz=uvw(i3,4,3)-uvw(i1,4,3)
         d2=dx*dx+dy*dy+dz*dz
         d=sqrt(d2)
         ex=uvw(i4,4,1)-uvw(i2,4,1)
         ey=uvw(i4,4,2)-uvw(i2,4,2)
         ez=uvw(i4,4,3)-uvw(i2,4,3)
         e2=ex*ex+ey*ey+ez*ez
         e=sqrt(e2)
c!mp  ang1=acos((a2+b2-d2)/(2*a*b))
c!mp  ang2=acos((b2+c2-e2)/(2*b*c))
c!mp  rc1=d/(2*sin(ang1))
c!mp  rc2=e/(2*sin(ang2))
c!mp  mp Rewrite formula to replace 4 trigonometric function calls
c!mp  mp with 2 calls to sqrt.
      ca1=(a2+b2-d2)/(2*a*b)
      ca2=(b2+c2-e2)/(2*b*c)
      rc1=0.5*d/sqrt(1-ca1*ca1)
      rc2=0.5*e/sqrt(1-ca2*ca2)
      rho=sqrt(rc1*rc2)
      curv=cfac/rho
      bpdat(i2,20)=curv
      cusum=cusum+curv
      cusum2=cusum2+curv*curv
      risum=risum+bpdat(i2,18)
c-------------------------------------find vectors to centers of circles
      px=a2*bx+b2*ax
      py=a2*by+b2*ay
      pz=a2*bz+b2*az
      qx=ay*bz-az*by
      qy=az*bx-ax*bz
      qz=ax*by-ay*bx
      q2=qx*qx+qy*qy+qz*qz
      q=sqrt(q2)
         v1x=(qy*pz-qz*py)/(2*q2)
         v1y=(qz*px-qx*pz)/(2*q2)
         v1z=(qx*py-qy*px)/(2*q2)

      px=b2*cx+c2*bx
      py=b2*cy+c2*by
      pz=b2*cz+c2*bz
      tx=by*cz-bz*cy
      ty=bz*cx-bx*cz
      tz=bx*cy-by*cx
      t2=tx*tx+ty*ty+tz*tz
      t=sqrt(t2)
         v2x=(ty*pz-tz*py)/(2*t2)
         v2y=(tz*px-tx*pz)/(2*t2)
         v2z=(tx*py-ty*px)/(2*t2)
c---------------------------------------------------------------register
c!mp moved v calculation after screw to ensure use of mid-step
c!mp frame urot(3,:).
            do m=1,4
            do j=1,3
            r1(m,j)=uvw(i2,m,j)
            r2(m,j)=uvw(i3,m,j)
            enddo
            enddo
            call screw(r1,r2,1)
        vx=v1x/rc1+v2x/rc2
        vy=v1y/rc1+v2y/rc2
        vz=v1z/rc1+v2z/rc2
        dot=vx*urot(3,1)+vy*urot(3,2)+vz*urot(3,3)
        vxp=vx-dot*urot(3,1)
        vyp=vy-dot*urot(3,2)
        vzp=vz-dot*urot(3,3)
        rp=sqrt(vxp*vxp+vyp*vyp+vzp*vzp)
        vxp=vxp/rp
        vyp=vyp/rp
        vzp=vzp/rp
      dot=vxp*urot(2,1)+vyp*urot(2,2)+vzp*urot(2,3)
      reg=acos(dot)
      ax=urot(2,2)*vzp-urot(2,3)*vyp
      ay=urot(2,3)*vxp-urot(2,1)*vzp
      az=urot(2,1)*vyp-urot(2,2)*vxp
      dot=ax*urot(3,1)+ay*urot(3,2)+az*urot(3,3)
      if(dot.lt.0.) reg=-reg
      bpdat(i2,21)=reg*crd 
c-------------------------------------------------------------vector O/P
      if(.not.traj) then
      write(6,15) i2,na(i2,1),nu(i2,1),na(i3,1),nu(i3,1),curv,
     1 rho,reg*crd
15    format(2x,i3,') ',a1,i4,'/',a1,i4,2x,3f9.3)
      x0=urot(4,1)
      y0=urot(4,2)
      z0=urot(4,3)
      pfac=rvfac*curv
         kas=kas+1
         cors(kas,1)=x0
         cors(kas,2)=y0
         cors(kas,3)=z0
         snam(kas)='C'
         sunit(kas)='ROC'
         nunis(kas)=i
         kas=kas+1
         cors(kas,1)=x0+pfac*vx
         cors(kas,2)=y0+pfac*vy
         cors(kas,3)=z0+pfac*vz
         snam(kas)='C'
         sunit(kas)='ROC'
         nunis(kas)=i
         mats(khs+1)=80000
         mats(khs+2)=kas-1
         mats(khs+3)=kas
         khs=khs+3
         endif
      enddo
c-------------------------------------------------------single point O/P
      if(.not.traj) then
      if(circ) then
      cuave=cusum/nlev
      cudev=sqrt(cusum2/nlev-cuave*cuave)
      cucir=cfac*2*pi/risum
      write(6,30) cuave,cudev,cucir
30    format(/2x,'<curvature> = ',f7.3,' std. dev. = ',f7.3,
     1 ' curvature eqv. circle = ',f7.3)
         else
         cuave=cusum/(nlev-3)
         cudev=sqrt(cusum2/(nlev-3)-cuave*cuave)
         write(6,32) cuave,cudev
32       format(/2x,'<curvature> = ',f7.3,' std. dev. = ',f7.3)
         endif
      call pdbout('_R',1)
      endif
      return
      end
