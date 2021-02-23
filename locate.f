      subroutine locate(itj,skip)
      include 'curves_data.inc'
      character*1 na,nt,base,type,first
      character*4 mnam,munit,name,ban,nback
      character*8 ibnam
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 nind,traj,skip
      dimension dcor(40,3),ds(9),ss(3),a(3,3),nind(40)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/sto/bpdat(n3,n5),traj
 
crl   Find necessary atoms (base and backbone), optiionally LS fit standard base
crl   Generate base ref. system ref(level,strand,vec,xyz)
crl   vec1 - dyad (into major), vec2 - perp. dyad to 5'-3', vec3 - normal, vec4
      skip=.false. 
      refmax=0.
      refmin=0.
      rmsx=0.d0
c------------------------------------------------find base reference pts
      do n=1,nst
      do m=1,nlev
      iu=abs(ni(m,n))
         if(iu.eq.0) then !...no nucleotide present
         na(m,n)='-'
         nu(m,n)=0
         goto 20
         endif
 
crl   ncen defines limits of subunits iu (nucleotides)
 
      is=ncen(iu-1)+1
      ie=ncen(iu)
      numb=nunit(is)
         ifir=ichar(munit(is)(:1))
         if(ifir.ge.97.and.ifir.le.122) munit(is)(:1)=char(ifir-32)
         first=munit(is)
      na(m,n)=first
      nu(m,n)=numb
      do lb=1,nbas
      if(first.eq.base(lb)) goto 10
      enddo
      if(itj.eq.itst) write(6,8) m,n,munit(is)
8     format(/2x,'Base not recognized at level ',i3,' strand ',i1,
     1 ': ',a4)
      na(m,n)='-'
      nu(m,n)=0
      goto 20
10    kb=0
      nt(m,n)=type(lb)
         do mm=1,ibref(lb)
         nind(mm)=.false.
         enddo
      do i=is,ie
      name=mnam(i)
         do mm=1,ibref(lb)
         if(name.eq.ban(lb,mm).and..not.nind(mm)) then
         nind(mm)=.true.
         do j=1,3
         dcor(mm,j)=corm(i,j)
         enddo
         kb=kb+1
         endif
         enddo
      enddo
      do mm=1,ibref(lb)
      if(.not.nind(mm)) then
      if(mm.le.4) then
      write(6,12) ban(lb,m),n,m,'(Necessary for analysis)'
12    format(/2x,' Base atom ',a4,' missing in strand ',
     1 i3,' level ',i3,a)
      stop
         else
         if(itj.eq.itst) write(6,12) ban(lb,m),n,m,
     1   '(not essential for analysis)'
         endif
      endif
      enddo
c------------------------------------------------ls fit of standard base
      if(fit.and.kb.ge.3) then
      call lsfit(lb,dcor,rms,nind)
      if(rms.gt.rmsx) rmsx=rms
         if(test) then
         if(n.eq.1.and.m.eq.1) write(6,50)
50       format(/2x,'LS fitting of standard bases ...',
     1        //3x,'Str',3x,'Pos',2x,'Base',10x,'Rms (ang)'/)
         write(6,60) n,m,munit(is),nunit(is),rms
60       format(2x,i3,' : ',i3,')  ',a4,1x,i3,4x,f7.3)
         endif
      endif
c------------------------------------------B-S bond and normal direction
      x0=dcor(2,1)
      y0=dcor(2,2)
      z0=dcor(2,3)
      ax=dcor(1,1)-x0
      ay=dcor(1,2)-y0
      az=dcor(1,3)-z0
      ra=sqrt(ax*ax+ay*ay+az*az)
      ax=ax/ra
      ay=ay/ra
      az=az/ra
         bx=dcor(3,1)-x0
         by=dcor(3,2)-y0
         bz=dcor(3,3)-z0
         cx=ay*bz-az*by
         cy=az*bx-ax*bz
         cz=ax*by-ay*bx
c-----------------------------------------------------------------normal
         do j=1,9
         if(j.le.3) ss(j)=0.d0
         ds(j)=0.d0
         enddo
      do i=1,ibref(lb)
      if(nind(i)) then
      do j=1,3
      ss(j)=ss(j)+dcor(i,j)
      enddo
      endif
      enddo
         do j=1,3
         ss(j)=ss(j)/kb
         enddo
      j=0
      do mm=1,3
      do l=1,mm
      j=j+1
         do i=1,ibref(lb)
         if(nind(i)) ds(j)=ds(j)+(dcor(i,mm)-ss(mm))*(dcor(i,l)-ss(l))
         enddo
      enddo
      enddo
      call eigen(ds,a,3)
      ra=sqrt(a(1,3)**2+a(2,3)**2+a(3,3)**2)
      rx=a(1,3)/ra
      ry=a(2,3)/ra
      rz=a(3,3)/ra
         dot=rx*cx+ry*cy+rz*cz
         if(dot.lt.0.) then
         rx=-rx
         ry=-ry
         rz=-rz
         endif
      ref(m,n,3,1)=rx
      ref(m,n,3,2)=ry
      ref(m,n,3,3)=rz
c----------------------------------------------construct reference point
 
crl  the distance dis and the angles th1 and th2 place the base reference point
crl  using the Tsukuba convention, or using Curves the convention if refo=.t.
 
      ca=cos(cdr*(th1))
      sa=sin(cdr*(th1))
      xx=ax*dis
      yy=ay*dis
      zz=az*dis
      ref(m,n,4,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1             (rx*rz*(1-ca)+ry*sa)*zz+x0
      ref(m,n,4,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1             (ry*rz*(1-ca)-rx*sa)*zz+y0
      ref(m,n,4,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1             (rz*rz+(1-rz*rz)*ca)*zz+z0
c----------------------------------------------------construct perp dyad
      cb=cos(cdr*(th2))
      sb=sin(cdr*(th2))
      xx=ax
      yy=ay
      zz=az
      dx=(rx*rx+(1-rx*rx)*cb)*xx+(rx*ry*(1-cb)-rz*sb)*yy+
     1   (rx*rz*(1-cb)+ry*sb)*zz
      dy=(rx*ry*(1-cb)+rz*sb)*xx+(ry*ry+(1-ry*ry)*cb)*yy+
     1   (ry*rz*(1-cb)-rx*sb)*zz
      dz=(rx*rz*(1-cb)-ry*sb)*xx+(ry*rz*(1-cb)+rx*sb)*yy+
     1   (rz*rz+(1-rz*rz)*cb)*zz
      ref(m,n,2,1)=dx
      ref(m,n,2,2)=dy
      ref(m,n,2,3)=dz
c---------------------------------------------------------construct dyad
      ref(m,n,1,1)=dy*rz-dz*ry
      ref(m,n,1,2)=dz*rx-dx*rz
      ref(m,n,1,3)=dx*ry-dy*rx
c--------------------------------------------------atom for groove depth
      ref(m,n,5,1)=dcor(4,1)
      ref(m,n,5,2)=dcor(4,2)
      ref(m,n,5,3)=dcor(4,3)
 
      do i=1,4
      do j=1,3
      if(ref(m,n,i,j).gt.refmax) refmax=ref(m,n,i,j)
      if(ref(m,n,i,j).lt.refmin) refmin=ref(m,n,i,j)
      enddo
      enddo
 
      if(frames) write(8,22) n,m,
     1 ((nint(ref(m,n,i,j)*1000),j=1,3),i=1,4)
22    format(14i7)
20    continue
      enddo
      enddo
c---------------------------------------------------------------------LS
      if(fit.and..not.test.and..not.traj) write(6,70) rmsx
70    format(/2x,'LS fitting of standard bases ...RMS max = ',f6.3)
 
      if(refmax.gt.9999.or.refmin.lt.-999) then
      if(traj) then
      write(6,23) itj
23    format(/2x,'---- I/P coord. error, snapshot ',i7,' skipped')
      skip=.true.
         else
         write(6,24)
24       format(/2x,'---- I/P coord. error ----')
         stop
         endif
      endif
      return
      end
