      subroutine ligloc(itj)
      include 'curves_data.inc'
      character*1 na,nt,base,type,first
      character*4 mnam,munit,name,ban,bln,liga,nback
      character*8 ibnam
      character*20 ilnam
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,nind,
     1 traj,start,need(20)
      dimension dcor(40,3),ds(9),ss(3),a(3,3),nind(40)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/lgd/ilnam(20),paxe(20),blef(20,40,3),parl(20,9),
     1 ilref(20),iluni(20,2),ilfrm(20,-3:20),bln(20,40),nlig,liga(20)
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/sto/bpdat(n3,n5),traj
      start=.true.
      kl=0
c----------------------------------------------find ligand reference pts
      do k=1,kcen
      is=ncen(k-1)+1
      ie=ncen(k)
      do n=1,nlig
      if(munit(is)(:3).eq.liga(n).and.
     1 ((nunit(is).ge.iluni(n,1).and.nunit(is).le.iluni(n,2)).or.
     1 (iluni(n,1).eq.0.and.iluni(n,2).eq.0))) then
         km=0
         do mm=1,ilref(n)
         nind(mm)=.false.
         need(mm)=.false.
         enddo
         do mm=1,ilfrm(n,0)
         need(ilfrm(n,mm))=.true.
         enddo
      do i=is,ie
      name=mnam(i)
         do mm=1,ilref(n)
         if(name.eq.bln(n,mm).and..not.nind(mm)) then
         nind(mm)=.true.
         do j=1,3
         dcor(mm,j)=corm(i,j)
         enddo
         km=km+1
         endif
         enddo
      enddo
      do mm=1,ilref(n)
      if(.not.nind(mm)) then
      if(need(mm)) then
      write(6,12) bln(n,mm),'(Necessary for analysis)'
12    format(/2x,' Ligand atom ',a4,' missing ',a)
      stop
         else
         if(itj.eq.itst) write(6,12) bln(n,mm),
     1   '(not essential for analysis)'
         endif
      endif
      enddo
c------------------------------------------------ls fit of standard base
      if(fit.and.km.ge.3) then
      call lslig(n,dcor,rms,nind)
         if(test) then
         if(start) write(6,50)
50       format(/2x,'LS fitting of ligands ...',
     1        //2x,'Ligand',8x,'Rms (ang)'/)
         write(6,60) munit(is)(:3),nunit(is),rms,ilnam(n)
60       format(2x,a4,1x,i3,4x,f7.3,3x,a)
         start=.false.
         endif
      endif
c-------------------------------------------setup ligand reference frame
      kl=kl+1
      ilig(kl)=is
      if(abs(paxe(n)).gt.nlev)
     1 stop '---- Ligand paxe cannot be no. of levels ----'
      plig(kl)=paxe(n)
      m1=ilfrm(n,-3)
      m2=ilfrm(n,-2)
      m3=ilfrm(n,-1)
      m4=ilfrm(n,0)
      x0=0.
      y0=0.
      z0=0.
      do mm=1,m4
      j=ilfrm(n,mm)
      x0=x0+dcor(j,1)
      y0=y0+dcor(j,2)
      z0=z0+dcor(j,3)
      enddo
      rel(kl,4,1)=x0/m4
      rel(kl,4,2)=y0/m4
      rel(kl,4,3)=z0/m4
      ax=dcor(m2,1)-dcor(m1,1)
      ay=dcor(m2,2)-dcor(m1,2)
      az=dcor(m2,3)-dcor(m1,3)
      ra=sqrt(ax*ax+ay*ay+az*az)
      ax=ax/ra
      ay=ay/ra
      az=az/ra
      rel(kl,2,1)=ax
      rel(kl,2,2)=ay
      rel(kl,2,3)=az
      bx=dcor(m3,1)-dcor(m1,1)
      by=dcor(m3,2)-dcor(m1,2)
      bz=dcor(m3,3)-dcor(m1,3)
      dot=ax*bx+ay*by+az*bz
      cx=bx-ax*dot
      cy=by-ay*dot
      cz=bz-az*dot
      rc=sqrt(cx*cx+cy*cy+cz*cz)
      cx=cx/rc
      cy=cy/rc
      cz=cz/rc
      rel(kl,1,1)=cx
      rel(kl,1,2)=cy
      rel(kl,1,3)=cz
      rel(kl,3,1)=cy*az-cz*ay
      rel(kl,3,2)=cz*ax-cx*az
      rel(kl,3,3)=cx*ay-cy*ax
c-----------------------------------------------------------------------
      endif
      enddo
      enddo
      klig=kl
      return
      end
