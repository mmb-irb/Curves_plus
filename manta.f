      subroutine manta
      include 'curves_data.inc'
      parameter (nbac=10) !... note: nbac must be even
      character*1 na,nt,nke,chain(8)
      character*4 mnam,munit,snam,sunit,nat,nback,ext(3)
      character*6 nbo,stg(n3,8)
      character*7 string
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,locp,
     1 rhelx,lpa,traj
      dimension bkc(n3*nbac,3,4),bpc(n3,3,4),bpt(n3,3,4),
     1 pl(3),pu(3),ul(3),uu(3),f(4,3),wid(n2*5,n2*5),
     1 wex(2),dex(2),iex(2),jex(2)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/sto/bpdat(n3,n5),traj
      common/str/nbo(40),nat(40,8),ish(40,8),ntor,nke(40)
      data wmax/18.0/,ext/'.min','.dia','.maj'/,
     1 chain/'A','B','C','D','E','F','G','H'/
      sq2=sqrt(2.d0)
      kwid=100
      rhelx=.true.
      joff=19+nst*(ntor-1)
c-----------------------------------------------------------------------
      if(.not.traj) write(6,1)
1     format(/2x,'(E) Groove parameters'/)
      if(.not.circ) then
      ntot=nbac*(nlev-2)+1
      delr=float(nlev-2)/(ntot-1)
         else
         ntot=nbac*nlev
         delr=float(nlev)/ntot
         endif
c--------------------------------------------------find phosphorus atoms
      do ks=1,nst
      l=1
      ms=1
      me=nlev
      if(.not.circ) then
      if(idr(ks).lt.0) then
      me=nlev-1
      else
      ms=2
      endif
         else
         if(idr(ks).lt.0) l=nlev
      endif
      do m=ms,me
      locp=.false.
      mu=abs(ni(m,ks))
      if(mu.ne.0) then
      do i=ncen(mu-1)+1,ncen(mu)
         if(mnam(i).eq.nback(ks)) then
         do j=1,3
         bpc(l,j,ks)=corm(i,j)
         enddo
         locp=.true.
         endif
      enddo
      l=l+1
      if(l.gt.nlev) l=1
      endif
         if(.not.locp.and..not.traj) then
         write(6,5) nback(ks),ks,m
5        format(4x,'No groove analysis: ',a,' absent in strand = ',i1,
     $   ' at level = ',i3)
         return
         endif
      enddo
c------------------------------------------------build backbone tangents
      lstart=2
      lend=nlev-2
        if(circ) then
        lstart=1
        lend=nlev
        endif
      do l=lstart,lend
      ldo=l-1
      if(ldo.lt.1) ldo=ldo+(nlev-1)
      lup=l+1
      if(lup.gt.nlev-1) lup=lup-(nlev-1)
      ax=bpc(l,1,ks)-bpc(ldo,1,ks) !... y-x
      ay=bpc(l,2,ks)-bpc(ldo,2,ks)
      az=bpc(l,3,ks)-bpc(ldo,3,ks)
      ra2=ax*ax+ay*ay+az*az
      ra=sqrt(ra2)
      bx=bpc(lup,1,ks)-bpc(l,1,ks) !... z-y
      by=bpc(lup,2,ks)-bpc(l,2,ks)
      bz=bpc(lup,3,ks)-bpc(l,3,ks)
      rb2=bx*bx+by*by+bz*bz
      rb=sqrt(rb2)
      cx=bpc(lup,1,ks)-bpc(ldo,1,ks) !... z-x
      cy=bpc(lup,2,ks)-bpc(ldo,2,ks)
      cz=bpc(lup,3,ks)-bpc(ldo,3,ks)
      rc2=cx*cx+cy*cy+cz*cz
      rc=sqrt(rc2)
      th=acos((bx*cx+by*cy+bz*cz)/(rb*rc))
      r=ra/(2*sin(th))
         ex=ax/ra
         ey=ay/ra
         ez=az/ra
         gam=ex*bx+ey*by+ez*bz
         epx=bx-gam*ex
         epy=by-gam*ey
         epz=bz-gam*ez
         ep=sqrt(epx*epx+epy*epy+epz*epz)
         epx=epx/ep
         epy=epy/ep
         epz=epz/ep
         bpt(l,1,ks)=cos(th)*ex+sin(th)*epx
         bpt(l,2,ks)=cos(th)*ey+sin(th)*epy
         bpt(l,3,ks)=cos(th)*ez+sin(th)*epz
            if(.not.circ) then
            if(l.eq.2) then
            bpt(1,1,ks)=cos(th)*ex-sin(th)*epx
            bpt(1,2,ks)=cos(th)*ey-sin(th)*epy
            bpt(1,3,ks)=cos(th)*ez-sin(th)*epz
               else if(l.eq.nlev-2) then
               x=bpt(l,1,ks)
               y=bpt(l,2,ks)
               z=bpt(l,3,ks)
               bx=sq2*bx/rb
               by=sq2*by/rb
               bz=sq2*bz/rb
               bpt(nlev-1,1,ks)=x*(bx*bx-1)+y* bx*by   +z* bx*bz
               bpt(nlev-1,2,ks)=x* by*bx   +y*(by*by-1)+z* by*bz
               bpt(nlev-1,3,ks)=x* bz*bx   +y* bz*by   +z*(bz*bz-1)
               endif
            endif
      enddo
c---------------------------------------------build backbone spline axe
      do i=1,ntot
      r=1+(i-1)*delr
      l=int(r)
      if(.not.circ.and.l.eq.nlev-1) l=nlev-2
      if(circ.and.l.eq.nlev+1) l=nlev
         lup=l+1
         if(lup.gt.nlev) lup=lup-nlev
         do j=1,3
         pl(j)=bpc(l,j,ks)
         pu(j)=bpc(lup,j,ks)
         ul(j)=bpt(l,j,ks)
         uu(j)=bpt(lup,j,ks)
         enddo
         gl=sqrt((pu(1)-pl(1))**2+(pu(2)-pl(2))**2+(pu(3)-pl(3))**2)
         do j=1,3
         f(1,j)=pl(j)
         f(2,j)=ul(j)
         f(3,j)=( 3/gl**2)*(pu(j)-pl(j))-(1/gl   )*(uu(j)+2*ul(j))
         f(4,j)=(-2/gl**3)*(pu(j)-pl(j))+(1/gl**2)*(uu(j)+  ul(j))
         enddo
      t=gl*(r-l)
      do j=1,3
      bkc(i,j,ks)=f(1,j)+f(2,j)*t+f(3,j)*t**2+f(4,j)*t**3
      enddo
      enddo ! ntot pts. along backbone
      enddo ! ks strands
c-------------------------------------------------------------pdb output
         ns=0
         lm=0
         do ks=1,nst
         do n=1,ntot
         ns=ns+1
         lm=lm+1
         mats(lm)=ns
         snam(ns)='O'
         sunit(ns)='BAC'//chain(ks)
         nunis(ns)=ks
         do j=1,3
         cors(ns,j)=bkc(n,j,ks)
         enddo
         enddo
            if(ks.lt.nst) then
            lm=lm+1
            mats(lm)=80000
            endif
         enddo
      if(nst.eq.1) then
      kas=ns
      khs=lm
      if(.not.traj) call pdbout('_B',1)
      return
      endif
c-----------------------------------------------------------------planes
      kup=nst
      if(nst.eq.2) kup=1
      do kst=1,kup !...loop over strand pairs, unless duplex
      ksl=kst
      ksu=kst+1
      if(ksu.gt.nst) ksu=1
 
      do i=1,ntot
      x1=bkc(i,1,ksl)
      y1=bkc(i,2,ksl)
      z1=bkc(i,3,ksl)
         do j=1,ntot
         x2=bkc(j,1,ksu)
         y2=bkc(j,2,ksu)
         z2=bkc(j,3,ksu)
         cx=x1-x2
         cy=y1-y2
         cz=z1-z2
         rc2=cx*cx+cy*cy+cz*cz
         wid(i,j)=sqrt(rc2)
         enddo
      enddo
 
c-------------------------------------------------test O/P of wid matrix
      if(test.and..not.traj) then
      write(string,'(''_'',i1,'':'',i1)') ksl,ksu
      open(4,file='wid'//string(:4)//'.txt',status='unknown')
      do i=1,ntot
      write(4,'(1000f7.2)') (wid(i,j),j=1,ntot)
      enddo
      close(4)
      endif
c---------------------------------------------------------detect grooves
      kop=0
      do k=2,2*ntot,nbac
      kop=kop+1
      lev=1.5+(kop-1)*0.5
      if(lev.gt.nlev) lev=lev-nlev
c----------------------------------------------setup axis point and dyad
         if(mod(kop,2).eq.0) then ! we are at a base level
         dax=uvw(lev,4,1)
         day=uvw(lev,4,2)
         daz=uvw(lev,4,3)
         uax=uvw(lev,1,1)
         uay=uvw(lev,1,2)
         uaz=uvw(lev,1,3)
         dbx=(ref(lev,ksl,5,1)+ref(lev,ksu,5,1))/2
         dby=(ref(lev,ksl,5,2)+ref(lev,ksu,5,2))/2
         dbz=(ref(lev,ksl,5,3)+ref(lev,ksu,5,3))/2
            else ! we are halfway between levels
            levu=lev+1
            if(circ.and.levu.gt.nlev) levu=1
            dax=(uvw(lev,4,1)+uvw(levu,4,1))/2
            day=(uvw(lev,4,2)+uvw(levu,4,2))/2
            daz=(uvw(lev,4,3)+uvw(levu,4,3))/2
            uax=(uvw(lev,1,1)+uvw(levu,1,1))
            uay=(uvw(lev,1,2)+uvw(levu,1,2))
            uaz=(uvw(lev,1,3)+uvw(levu,1,3))
            ru=sqrt(uax*uax+uay*uay+uaz*uaz)
            uax=uax/ru
            uay=uay/ru
            uaz=uaz/ru
            dbx=(ref(lev,ksl,5,1) +ref(lev,ksu,5,1)
     &         +ref(levu,ksl,5,1)+ref(levu,ksu,5,1))/4
            dby=(ref(lev,ksl,5,2) +ref(lev,ksu,5,2)
     &         +ref(levu,ksl,5,2)+ref(levu,ksu,5,2))/4
            dbz=(ref(lev,ksl,5,3) +ref(lev,ksu,5,3)
     &         +ref(levu,ksl,5,3)+ref(levu,ksu,5,3))/4
            endif
c-------------------------------------------------------find best minima
      px=(bkc(k/2,1,ksl)+bkc(k/2,1,ksu))/2
      py=(bkc(k/2,2,ksl)+bkc(k/2,2,ksu))/2
      pz=(bkc(k/2,3,ksl)+bkc(k/2,3,ksu))/2
      refx=px-dax
      refy=py-day
      refz=pz-daz
      rd=sqrt(refx*refx+refy*refy+refz*refz)
      if(nst.eq.2.or.rd.lt.0.5) then
      refx=-uax
      refy=-uay
      refz=-uaz
      endif
            do i=1,2
            iex(i)=0
            jex(i)=0
            wex(i)=0.
            dex(i)=0.
            enddo
      wmip=100.
      wmin=100.
      do il=(k-kwid)/2,(k+kwid)/2
      i=il
      j=k-i
      if(circ) then
      if(i.lt.1) i=i+ntot
      if(j.lt.1) j=j+ntot
      if(i.gt.ntot) i=i-ntot
      if(j.gt.ntot) j=j-ntot
      endif
      if(circ.or.
     1   (i.gt.2.and.i.lt.ntot-1.and.j.gt.2.and.j.lt.ntot-1)) then
      w=wid(i,j)
      ido=i-1
      iup=i+1
      jdo=j-1
      jup=j+1
         if(circ) then
         if(ido.lt.1) ido=ido+ntot
         if(jdo.lt.1) jdo=jdo+ntot
         if(iup.gt.ntot) iup=iup-ntot
         if(jup.gt.ntot) jup=jup-ntot
         endif
      if(w.lt.wid(ido,jup).and.w.lt.wid(iup,jdo)) then
      ind=ind+1
      px=(bkc(i,1,ksl)+bkc(j,1,ksu))/2
      py=(bkc(i,2,ksl)+bkc(j,2,ksu))/2
      pz=(bkc(i,3,ksl)+bkc(j,3,ksu))/2
      dx=px-dbx
      dy=py-dby
      dz=pz-dbz
      rd=sqrt(dx*dx+dy*dy+dz*dz)
      dot=(px-dax)*refx+(py-day)*refy+(pz-daz)*refz
c        if(dot.gt.0.and.w.lt.wmip) then
         if(dot.gt.0.and.w.lt.wmip.and.w-2*wback.lt.wmax) then
         wmip=w
         iex(1)=i
         jex(1)=j
         wex(1)=w-2*wback
         dex(1)=rd-wbase
         endif
c        if(dot.lt.0.and.w.lt.wmin) then
         if(dot.lt.0.and.w.lt.wmin.and.w-2*wback.lt.wmax) then
         wmin=w
         iex(2)=i
         jex(2)=j
         wex(2)=w-2*wback
         dex(2)=rd-wbase
         endif
      endif
      endif
      enddo
c------------------------------------------------------extra PDB output
         mup=2
         if(nst.gt.2) mup=1
         do m=1,mup
         mch=kst+nst
         if(mup.eq.2) mch=m+nst
         if(iex(m).ne.0) then
         i=iex(m)
         j=jex(m)
         ns=ns+1
         snam(ns)='B'
         sunit(ns)='GRV'//chain(mch)
         nunis(ns)=mch
         do l=1,3
         cors(ns,l)=bkc(i,l,ksl)
         enddo
         ns=ns+1
         snam(ns)='B'
         sunit(ns)='GRV'//chain(mch)
         nunis(ns)=mch
         do l=1,3
         cors(ns,l)=bkc(j,l,ksu)
         enddo
         lm=lm+1
         mats(lm)=80000
         lm=lm+1
         mats(lm)=ns-1
         lm=lm+1
         mats(lm)=ns
         endif
         enddo
c----------------------------------------------------------create output
      if(nst.eq.2) then ! need two grooves
      if(iex(1).gt.0) then
      write(stg(kop,1),'(f6.1)')  wex(1)
      write(stg(kop,2),'(f6.1)') dex(1)
         if(mod(kop,2).eq.0) then
         bpdat(lev,joff+1)=wex(1)
         bpdat(lev,joff+2)=dex(1)
         endif
      else
      stg(kop,1)='      '
      stg(kop,2)='      '
      endif
 
      if(iex(2).gt.0) then
      write(stg(kop,3),'(f6.1)')  wex(2)
      write(stg(kop,4),'(f6.1)') dex(2)
         if(mod(kop,2).eq.0) then
         bpdat(lev,joff+3)=wex(2)
         bpdat(lev,joff+4)=dex(2)
         endif
      else
      stg(kop,3)='      '
      stg(kop,4)='      '
      endif
 
      else ! nst > 2 only need one groove
         if(iex(1).gt.0) then
         write(stg(kop,2*kst-1),'(f6.1)') wex(1)
         write(stg(kop,2*kst),'(f6.1)') dex(1)
            if(mod(kop,2).eq.0) then
            bpdat(lev,joff+2*kst-1)=wex(1)
            bpdat(lev,joff+2*kst)=dex(1)
            endif
         else
         stg(kop,2*kst-1)='      '
         stg(kop,2*kst)  ='      '
         endif
      endif
      enddo !... K points along diagonal of wid matrix
      enddo !... pairs of backbones
c---------------------------------------------------------output grooves
         kas=ns
         khs=lm
         if(.not.traj) call pdbout('_B',1)
 
      if(.not.traj) then
      write(6,35) ('W',ks,ks+1,'D',ks,ks+1,ks=1,nst-1),
     1 'W',nst,1,'D',nst,1
35    format(2x,' Level',8x,8(:,3x,a1,2i1,2x))
      write(6,*)
      endif
      kop=0
      do k=2,2*ntot,nbac
      kop=kop+1
      lev=1+kop/2
      if(lev.gt.nlev) lev=lev-nlev
      r=1.5+float(k-2)/(2*nbac)
      if(mod(kop,2).eq.0) then
      write(string,'(2x,a1,i4)') na(lev,1),nu(lev,1)
         else
         string='       '
         endif
      if(.not.traj) write(6,40) r,string,(stg(kop,m),m=1,2*nst)
40    format(2x,f5.1,a7,8(2x,a6))
      enddo
 
      return
      end
