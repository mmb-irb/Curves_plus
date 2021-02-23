      subroutine params(totbend)
      include 'curves_data.inc'
      character*1 na,nt
      character*4 mnam,munit,snam,sunit,nback,inam
      character*8 sugt(10),puck,axang
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      character clin*56,strand(4)*3
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,lpa,
     1 traj,invert(n3)
      dimension r1(4,3),r2(4,3),sums(8),sumc(8),jtran(7),
     1 vaf(n3,2),mad(18)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/ion/pari(n1,3),ilis(n1,2),klis(40),ilib(40),kion(5),
     1 kisa,nion,nspl,inam(40)
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
      common/sto/bpdat(n3,n5),traj
      data strand/'1st','2nd','3rd','4th'/
      data sugt/'C3''-endo','C4''-exo ','O1''-endo','C1''-exo ',
     1'C2''-endo','C3''-exo ','C4''-endo','O1''-exo ',
     1'C1''-endo','C2''-exo '/
      data jtran/13,9,10,7,8,11,12/
      data mad/1,2,80000,1,3,80000,1,4,80000,1,5,6,80000,5,7,80000,5,8/
      nspl=200
c---------------------------------------------------------dyad inversion
      do k=2,nst
      do i=1,nlev
      if(ni(i,k).ne.0) then
      if(idr(k).lt.0) then
      do m=2,3
      do j=1,3
      ref(i,k,m,j)=-ref(i,k,m,j)
      enddo
      enddo
         else
         do m=1,2
         do j=1,3
         ref(i,k,m,j)=-ref(i,k,m,j)
         enddo
         enddo
         endif
      endif
      enddo
      enddo
c--------------------------------------------------------------setup upm
      do i=1,nlev
      if(ni(i,1).eq.0.or.ni(i,2).eq.0) then
      mj=0
      do j=1,nst
      if(ni(i,j).gt.0.and.mj.eq.0) mj=j !... strand with base present
      enddo
      do k=1,4
      do j=1,3
      upm(i,k,j)=ref(i,mj,k,j)
      enddo
      enddo
         else !... bases exist in strands 1-2
         do m=1,4
         do j=1,3
         r1(m,j)=ref(i,1,m,j)
         r2(m,j)=ref(i,2,m,j)
         enddo
         enddo
         call screw(r2,r1,0)
         do k=1,4
         do j=1,3
         upm(i,k,j)=urot(k,j)
         enddo
         enddo
         endif
      enddo
c-------------------------------------calculate helical axis (needs upm)
      call axis
      call smooth
c---------------------------------------------helical axis inversion (Z)
      iup=nlev
      if(circ) iup=nlev+1
      do i=2,iup
      il=i-1
      iu=i
      if(iu.gt.nlev) iu=1
      dx=upm(iu,4,1)-upm(il,4,1)
      dy=upm(iu,4,2)-upm(il,4,2)
      dz=upm(iu,4,3)-upm(il,4,3)
      dot=dx*upm(iu,3,1)+dy*upm(iu,3,2)+dz*upm(iu,3,3)
      invert(il)=.false.
      if(dot.lt.0) invert(il)=.true.
      enddo
      invert(iu)=invert(il)
c---------------------------------------------------------bp-axis params
      if(.not.traj) write(6,3)
3     format(/2x,'(A) BP-Axis',8x,
     1 'Xdisp',3x,'Ydisp',3x,'Inclin',4x,'Tip',2x,'Ax-bend'/)
         do j=1,6
         sums(j)=0.d0
         sumc(j)=0.d0
         enddo
      ic=0
      do i=1,nlev
      ic=ic+1
c-------------------------------simplified screw axis for bp-axis params
         dx=upm(i,4,1)-uvw(i,4,1)
         dy=upm(i,4,2)-uvw(i,4,2)
         dz=upm(i,4,3)-uvw(i,4,3)
         ct=upm(i,3,1)*uvw(i,3,1)+upm(i,3,2)*uvw(i,3,2)
     1     +upm(i,3,3)*uvw(i,3,3)
         theta=acos(ct)
         if(abs(theta).gt.1.d-4) then
         st=sin(theta)
         ux=upm(i,3,2)*uvw(i,3,3)-upm(i,3,3)*uvw(i,3,2)
         uy=upm(i,3,3)*uvw(i,3,1)-upm(i,3,1)*uvw(i,3,3)
         uz=upm(i,3,1)*uvw(i,3,2)-upm(i,3,2)*uvw(i,3,1)
         r=sqrt(ux*ux+uy*uy+uz*uz)
         ux=ux/r
         uy=uy/r
         uz=uz/r
         xx=upm(i,1,1)
         yy=upm(i,1,2)
         zz=upm(i,1,3)
         uvw(i,1,1)=(ux*ux+(1-ux*ux)*ct)*xx+(ux*uy*(1-ct)-uz*st)*yy+
     1              (ux*uz*(1-ct)+uy*st)*zz
         uvw(i,1,2)=(ux*uy*(1-ct)+uz*st)*xx+(uy*uy+(1-uy*uy)*ct)*yy+
     1              (uy*uz*(1-ct)-ux*st)*zz
         uvw(i,1,3)=(ux*uz*(1-ct)-uy*st)*xx+(uy*uz*(1-ct)+ux*st)*yy+
     1              (uz*uz+(1-uz*uz)*ct)*zz
         uvw(i,2,1)=uvw(i,3,2)*uvw(i,1,3)-uvw(i,3,3)*uvw(i,1,2)
         uvw(i,2,2)=uvw(i,3,3)*uvw(i,1,1)-uvw(i,3,1)*uvw(i,1,3)
         uvw(i,2,3)=uvw(i,3,1)*uvw(i,1,2)-uvw(i,3,2)*uvw(i,1,1)
            else
            do m=1,2
            do j=1,3
            uvw(i,m,j)=upm(i,m,j)
            enddo
            enddo
            endif
         theta=-theta*crd
         ux=theta*ux
         uy=theta*uy
         uz=theta*uz
         do j=1,3
         var(j)=  dx*uvw(i,j,1)+dy*uvw(i,j,2)+dz*uvw(i,j,3)
         var(j+3)=ux*uvw(i,j,1)+uy*uvw(i,j,2)+uz*uvw(i,j,3)
         enddo
c-----------------------------------------------------------------------
         if(invert(i)) then
         var(1)=-var(1)
         var(4)=-var(4)
         var(5)=var(5)-180.
         endif
         if(abs(var(4)).gt.180.0) var(4)=var(4)-sign(360.d0,var(4))
         if(abs(var(5)).gt.180.0) var(5)=var(5)-sign(360.d0,var(5))
         do j=1,3
         sums(j)=sums(j)+var(j)
         enddo
         do j=4,6
         sums(j)=sums(j)+sin(var(j)*cdr)
         sumc(j)=sumc(j)+cos(var(j)*cdr)
         enddo
      il=i-1
      if(i.eq.1.and.circ) il=nlev
      if(i.gt.1.or.circ) then
      dot=0.
      do j=1,3
      dot=dot+uvw(i,3,j)*uvw(il,3,j)
      enddo
      if(abs(dot).gt.1) dot=sign(1.d0,dot)
      stpbend=acos(dot)*crd
      bpdat(i,11)=stpbend
      write(axang,'(f8.1)') stpbend
         else
         stpbend=0.
         axang='     ---'
         endif
      bpdat(i,7)= var(1)
      bpdat(i,8)= var(2)
      bpdat(i,9)= var(4)
      bpdat(i,10)=var(5)
      if(.not.traj) then
      if(nst.eq.1) then
      write(6,29) i,na(i,1),nu(i,1),var(1),var(2),
     1 var(4),var(5),axang
29    format(2x,i3,') ',a1,i4,6x,2f8.2,2f8.1,a8)
          else
          write(6,30) i,na(i,1),nu(i,1),na(i,2),nu(i,2),var(1),var(2),
     1    var(4),var(5),axang
30        format(2x,i3,') ',a1,i4,'-',a1,i4,2f8.2,2f8.1,a8)
          endif
      endif
      enddo ! loop i
 
      do j=1,3
      sums(j)=sums(j)/ic
      enddo
      do j=4,6
      r=sqrt(sums(j)**2+sumc(j)**2)
      angle=acos(sumc(j)/r)*crd
      if(asin(sums(j)/r).lt.0.) angle=-angle
      sums(j)=angle
      enddo
      dot=0.
      do j=1,3
      dot=dot+uvw(naxlim+1,3,j)*uvw(nlev-naxlim,3,j)
      enddo
      if(abs(dot).gt.1) dot=sign(1.d0,dot)
      totbend=acos(dot)*crd
      if(.not.traj) write(6,14) sums(1),sums(2),sums(4),sums(5),
     1 totbend,naxlim+1,nlev-naxlim
14    format(/7x,'Average:',3x,2f8.2,2f8.1,'  Total bend =',f8.1,
     1 ' (',i3,' to ',i3,')')
 
c------------------------------------------------inter axis frame params
            do j=7,8
            sums(j)=0.d0
            sumc(j)=0.d0
            enddo
         ic=0
         iup=nlev
         if(circ) iup=nlev+1
         do i=2,iup
         iu=i
         if(i.gt.nlev) iu=1
         il=i-1
         ic=ic+1
            do m=1,4
            do j=1,3
            r1(m,j)=uvw(il,m,j)
            r2(m,j)=uvw(iu,m,j)
            enddo
            enddo
            call screw(r1,r2,0)
            if(invert(il)) then
            var(3)=-var(3)
            var(6)=-var(6)
            endif
         if(abs(var(4)).gt.180.0) var(4)=var(4)-sign(360.d0,var(4))
         if(abs(var(5)).gt.180.0) var(5)=var(5)-sign(360.d0,var(5))
         if(abs(var(6)).gt.180.0) var(6)=var(6)-sign(360.d0,var(6))
         vaf(il,1)=var(3)
         vaf(il,2)=var(6)
         sums(7)=sums(7)+var(3)
         sums(8)=sums(8)+sin(var(6)*cdr)
         sumc(8)=sumc(8)+cos(var(6)*cdr)
         enddo
 
         sums(7)=sums(7)/ic
         r=sqrt(sums(8)**2+sumc(8)**2)
         angle=acos(sumc(8)/r)*crd
         if(asin(sums(8)/r).lt.0.) angle=-angle
         sums(8)=angle
c--------------------------------------------------------intra bp params
      if(nst.gt.1) then
      if(.not.traj) write(6,11)
11    format(/2x,'(B) Intra-BP parameters')
      do ks=2,nst
      if(.not.traj) write(6,1) ks
1     format(/2x,'Strands 1-',i1,7x,'Shear',2x,'Stretch',1x,
     1 'Stagger',2x,'Buckle',2x,'Propel',1x,'Opening'/)
         do j=1,6
         sums(j)=0.d0
         sumc(j)=0.d0
         enddo
      ic=0
      do i=1,nlev
      if(ni(i,1).ne.0.and.ni(i,ks).ne.0) then
      ic=ic+1
         do m=1,4
         do j=1,3
         r1(m,j)=ref(i,1,m,j)
         r2(m,j)=ref(i,ks,m,j)
         enddo
         enddo
         call screw(r2,r1,0)
         if(invert(i)) then
         var(1)=-var(1)
         var(4)=-var(4)
         endif
         if(abs(var(4)).gt.180.0) var(4)=var(4)-sign(360.d0,var(4))
         if(abs(var(5)).gt.180.0) var(5)=var(5)-sign(360.d0,var(5))
         if(abs(var(6)).gt.180.0) var(6)=var(6)-sign(360.d0,var(6))
         do j=1,3
         sums(j)=sums(j)+var(j)
         enddo
         do j=4,6
         sums(j)=sums(j)+sin(var(j)*cdr)
         sumc(j)=sumc(j)+cos(var(j)*cdr)
         enddo
      do j=1,6
      bpdat(i,j)=var(j)
      enddo
      if(.not.traj) then
      write(6,10) i,na(i,1),nu(i,1),na(i,2),nu(i,2),(var(j),j=1,6)
10    format(2x,i3,') ',a1,i4,'-',a1,i4,3f8.2,3f8.1,2x,f8.2,f8.1)
      endif
      endif
      enddo
 
      do j=1,3
      sums(j)=sums(j)/ic
      enddo
      do j=4,6
      r=sqrt(sums(j)**2+sumc(j)**2)
      angle=acos(sumc(j)/r)*crd
      if(asin(sums(j)/r).lt.0.) angle=-angle
      sums(j)=angle
      enddo
      if(ic.ge.1.and..not.traj) write(6,12)(sums(j),j=1,6)
12    format(/7x,'Average:',3x,3f8.2,3f8.1,f8.2,f8.1)
      enddo
      endif
c--------------------------------------------------------inter bp params
      if(.not.traj) write(6,2)
2     format(/2x,'(C) Inter-BP',7x,'Shift',3x,'Slide',4x,'Rise',4x,
     1 'Tilt',4x,'Roll',3x,'Twist',3x,'H-Ris',3x,'H-Twi'/)
         do j=1,6
         sums(j)=0.d0
         sumc(j)=0.d0
         enddo
      ic=0
      iup=nlev
      if(circ) iup=nlev+1
      do i=2,iup
      iu=i
      il=i-1
      if(i.gt.nlev) iu=1
      ic=ic+1
         do m=1,4
         do j=1,3
         r1(m,j)=upm(il,m,j)
         r2(m,j)=upm(iu,m,j)
         enddo
         enddo
         call screw(r1,r2,0)
         if(invert(il)) then
         var(3)=-var(3)
         var(6)=-var(6)
         endif
         if(abs(var(4)).gt.180.0) var(4)=var(4)-sign(360.d0,var(4))
         if(abs(var(5)).gt.180.0) var(5)=var(5)-sign(360.d0,var(5))
         if(abs(var(6)).gt.180.0) var(6)=var(6)-sign(360.d0,var(6))
         do j=1,3
         sums(j)=sums(j)+var(j)
         enddo
         do j=4,6
         sums(j)=sums(j)+sin(var(j)*cdr)
         sumc(j)=sumc(j)+cos(var(j)*cdr)
         enddo
      do j=1,6
      bpdat(i-1,11+j)=var(j)
      enddo
      bpdat(i-1,18)=vaf(il,1)
      bpdat(i-1,19)=vaf(il,2)
      if(.not.traj) then
      write(6,20) i-1,na(il,1),nu(il,1),na(iu,1),nu(iu,1),
     1 (var(j),j=1,6),vaf(il,1),vaf(il,2)
20    format(2x,i3,') ',a1,i4,'/',a1,i4,3f8.2,3f8.1,f8.2,f8.1)
      endif
      enddo
      do j=1,3
      sums(j)=sums(j)/ic
      enddo
      do j=4,6
      r=sqrt(sums(j)**2+sumc(j)**2)
      angle=acos(sumc(j)/r)*crd
      if(asin(sums(j)/r).lt.0.) angle=-angle
      sums(j)=angle
      enddo
      if(.not.traj) write(6,12)(sums(j),j=1,8)
c-----------------------------------------------------------other params
      call axref
      if(.not.traj.and.klig.eq.0.and.kion(1).eq.0) call pdbout('_X',1)
      call backbo
      call manta(xdsp)
      if(klig.gt.0.or.kion(1).gt.0.or.ibld.ne.' ') call intaxe
         if(kion(1).gt.0) then
         call ionpar
         if(.not.traj) call pdbout('_X',1)
         endif
            if(klig.gt.0) then
            call ligpar
            if(.not.traj) call pdbout('_X',1)
            endif
      return
      end
