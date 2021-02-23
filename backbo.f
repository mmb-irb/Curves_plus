      subroutine backbo
      include 'curves_data.inc'
      parameter (mwid=10) !... number of torsion per line
      parameter (bondmax=2.5) !... maximum allowed bond length
      character*1 nke,na,nt,base,type
      character*4 nat,ban,mnam,munit,name
      character*6 nbo
      character*7 tnam,sugt(10),tor(40)
      character*8 ibnam
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 start,traj,nind(5)
      dimension is(3),ie(3),ita(5),store(5),tnam(40)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/sto/bpdat(n3,n5),traj
      common/str/nbo(40),nat(40,8),ish(40,8),ntor,nke(40)
      data sugt/'  C3''en','  C4''ex','  O1''en','  C1''ex','  C2''en',
     1          '  C3''ex','  C4''en','  O1''ex','  C1''en','  C2''ex'/
      if(.not.traj) write(6,41)
41    format(/2x,'(D) Backbone Parameters')
      do k=1,nst
      joff=(k-1)*(ntor-1)+19
      ms=1
      me=min(mwid,ntor)
5     start=.true.
      do i=1,nlev
      if(ni(i,k).ne.0) then
c-----------------------------------------------------setup unit indices
      iu=abs(ni(i,k))
      is(2)=ncen(iu-1)+1
      ie(2)=ncen(iu)
      iup=0
      ido=0
      if(idr(k).gt.0) then
      if(i.lt.nlev) iup=abs(ni(i+1,k))
      if(i.eq.nlev.and.circ) iup=abs(ni(1,k))
      if(i.gt.1) ido=abs(ni(i-1,k))
      if(i.eq.1.and.circ) ido=abs(ni(nlev,k))
         else
         if(i.gt.1) iup=abs(ni(i-1,k))
         if(i.eq.1.and.circ) iup=abs(ni(nlev,k))
         if(i.lt.nlev) ido=abs(ni(i+1,k))
         if(i.eq.nlev.and.circ) ido=abs(ni(1,k))
         endif
      if(iup.gt.0) then
      is(3)=ncen(iup-1)+1
      ie(3)=ncen(iup)
         else
         is(3)=0
         ie(3)=0
         endif
      if(ido.gt.0) then
      is(1)=ncen(ido-1)+1
      ie(1)=ncen(ido)
         else
         is(1)=0
         ie(1)=0
         endif
c---------------------------------------------------find necessary atoms
      do m=ms,me
      tnam(m)=nbo(m)
      write(tor(m),'(a7)') '   ----'
c---------------------------------------------------sugar pseudorotation
      if(nke(m).eq.'S') then
      jc=0
         do j=1,5
         nind(j)=.false.
         enddo
      do j=1,5
      name=nat(m,j)
         do l=is(2),ie(2)
         if(name.eq.mnam(l).and..not.nind(j)) then
         nind(j)=.true.
         ita(j)=l
         jc=jc+1
         goto 15
         endif
         enddo
15    continue
      enddo
      if(jc.lt.5) goto 50
      store(1)=torp(ita(1),ita(2),ita(3),ita(4))
      store(2)=torp(ita(2),ita(3),ita(4),ita(5))
      store(3)=torp(ita(3),ita(4),ita(5),ita(1))
      store(4)=torp(ita(4),ita(5),ita(1),ita(2))
      store(5)=torp(ita(5),ita(1),ita(2),ita(3))
      a=0.
      b=0.
      do j=1,5
      a=a+store(j)*cos(cdr*(144.*(j-1)))
      b=b-store(j)*sin(cdr*(144.*(j-1)))
      enddo
      a=a*0.4d0
      b=b*0.4d0
      amp=sqrt(a*a+b*b)
         if(amp.gt.0.) then
         cp=a/amp
         sp=b/amp
         if(abs(cp).gt.1.) cp=sign(1.d0,cp)
         pha=acos(cp)*crd
         if(sp.lt.0.) pha=-pha
            else
            pha=0.
            endif
      bpdat(i,joff+m-2)=pha
      bpdat(i,joff+m-1)=amp
      write(tor(m-2),'(f7.1)') pha
      write(tor(m-1),'(f7.1)') amp
      if(pha.lt.0) pha=360+pha
      tor(m)=sugt(1+int(pha/36.))
c---------------------------------------------------------other torsions
         else if(nke(m).ne.'-') then
         jc=0
         jd=0
         if(nke(m).eq.'B'.and.nt(i,k).eq.'Y') jd=4
            do j=1,4
            nind(j)=.false.
            enddo
         do j=1,4
         id=ish(m,j+jd)
         name=nat(m,j+jd)
            if(is(2+id).eq.0) goto 20
            do l=is(2+id),ie(2+id)
            if(name.eq.mnam(l).and..not.nind(j)) then
            nind(j)=.true.
            ita(j)=l
            jc=jc+1
            goto 20
            endif
            enddo
20       continue
         enddo
         if(jc.lt.4) goto 50
c---------------------------------------------------------distance check
         do j=2,4
         jl=ita(j-1)
         ju=ita(j)
         dx=corm(ju,1)-corm(jl,1)
         dy=corm(ju,2)-corm(jl,2)
         dz=corm(ju,3)-corm(jl,3)
         r=sqrt(dx*dx+dy*dy+dz*dz)
         if(r.gt.bondmax) goto 50
         enddo
c-----------------------------------------------------------------------
         tors=torp(ita(1),ita(2),ita(3),ita(4))
         bpdat(i,joff+m)=tors
         write(tor(m),'(f7.1)') tors
         endif
50    continue
      enddo !...m torsions for nucleotide
c-----------------------------------------------------------------output
      if(start.and..not.traj) then
      start=.false.
      write(6,24) k,(tnam(m),m=ms,me)
24    format(/2x,' Strand ',i1,5x,10a7/)
      endif
      if(.not.traj) write(6,26) i,na(i,k),nu(i,k),(tor(m),m=ms,me)
26    format(2x,i3,') ',a1,i4,2x,10a7)
c-----------------------------------------------------------------------
      endif
      enddo !... i levels
         if(me.lt.ntor) then
         ms=ms+mwid
         me=min(ntor,me+mwid)
         goto 5
         endif
      enddo !... k strands
      return
      end
