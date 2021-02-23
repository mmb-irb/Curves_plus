      subroutine intaxe
      include 'curves_data.inc'
      character*1 na,nt
      character*4 snam,sunit,inam
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,lpa
      dimension r1(4,3),r2(4,3),t(3,3),dr(3),v(3)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/ion/pari(n1,3),ilis(n1,2),klis(40),ilib(40),kion(5),
     1 kisa,nion,nspl,inam(40)
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
      n=0
      iup=nlev
      if(circ) iup=nlev+1
      do i=2,iup
      iu=i
      il=iu-1
      if(i.gt.nlev) iu=1
            do m=1,4
            do j=1,3
            r1(m,j)=uvw(il,m,j)
            r2(m,j)=uvw(iu,m,j)
            enddo
            enddo
            do j=1,3
            v(j)=r2(4,j)-r1(4,j)
            enddo
            call screw(r1,r2,0)
            proj=v(1)*usc(1)+v(2)*usc(2)+v(3)*usc(3)
 
      mup=nspl-1
      if(i.eq.iup) mup=mup+1
      do k=0,mup
      fth=theta*k/nspl
      fle=proj*k/nspl
      st=sin(fth)
      ct=cos(fth)
      cm=ct-1
      t(1,1)= ct       -cm*usc(1)*usc(1)
      t(1,2)=-st*usc(3)-cm*usc(1)*usc(2)
      t(1,3)= st*usc(2)-cm*usc(1)*usc(3)
      t(2,1)= st*usc(3)-cm*usc(2)*usc(1)
      t(2,2)= ct       -cm*usc(2)*usc(2)
      t(2,3)=-st*usc(1)-cm*usc(2)*usc(3)
      t(3,1)=-st*usc(2)-cm*usc(3)*usc(1)
      t(3,2)= st*usc(1)-cm*usc(3)*usc(2)
      t(3,3)= ct       -cm*usc(3)*usc(3)
         do j=1,3
         dr(j)=r1(4,j)-usc(j+3)
         enddo
      n=n+1
      do j=1,3
      uint(n,1,j)=r1(1,1)*t(j,1)+r1(1,2)*t(j,2)+r1(1,3)*t(j,3)
      uint(n,2,j)=r1(2,1)*t(j,1)+r1(2,2)*t(j,2)+r1(2,3)*t(j,3)
      uint(n,3,j)=r1(3,1)*t(j,1)+r1(3,2)*t(j,2)+r1(3,3)*t(j,3)
      uint(n,4,j)=  dr(1)*t(j,1)  +dr(2)*t(j,2)  +dr(3)*t(j,3)
     1            +usc(j+3)+fle*usc(j)
      enddo
         snam(n)='H'
         sunit(n)='AXIS'
         nunis(n)=il
         mats(n)=n
         do j=1,3
         cors(n,j)=uint(n,4,j)
         enddo
      enddo ! k points
      enddo ! i levels
         kas=n
         khs=n
      return
      end
