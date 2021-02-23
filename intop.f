      subroutine intop(eto,kto,inx)
      include 'curves_data.inc'
      character*3 eto
      character*4 name,mnam,munit,mnx,namer(n2)
      character*80 line
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      dimension inx(n1p),ncx(n1p),mnx(n1p),nd(10)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
c----------------------------------------------------------read pdb file
      if(eto.eq.'pdb') then
      i=0
      k=0
12    read(2,5,end=50) line
      if(line(:4).eq.'ATOM') then
      k=k+1
         if(index(line(13:16),'H').ne.0) then
         inx(k)=0
         goto 12
         endif
      i=i+1
      inx(k)=i
      read(line,11) mnam(i),munit(i),nunit(i)
11    format(12x,a4,1x,a3,2x,i4)
      endif
      goto 12
50    kam=i
      kto=k
         k=1
         ncen(0)=0
         name=munit(1)
         numb=nunit(1)
         do i=2,kam
         if(munit(i).ne.name.or.nunit(i).ne.numb) then
         ncen(k)=i-1
         k=k+1
         name=munit(i)
         numb=nunit(i)
         endif
         enddo
         ncen(k)=kam
         kcen=k
         goto 150
      endif
c-----------------------------------------------------------read new top
      read(2,5) line
      if(line(:1).eq."%") then
10    read(2,5,end=100) line
5     format(a)
      if(index(line,' POINTERS').ne.0) then
      read(2,5) line
      i1=index(line,'(')
      i2=index(line,')')
      read(2,line(i1:i2)) kam,(nd(i),i=1,10),kcen
      endif
      if(index(line,'ATOM_NAME').ne.0) then
      read(2,5) line
      i1=index(line,'(')
      i2=index(line,')')
      read(2,line(i1:i2)) (mnx(i),i=1,kam)
      endif
      if(index(line,'RESIDUE_LABEL').ne.0) then
      read(2,5) line
      i1=index(line,'(')
      i2=index(line,')')
      read(2,line(i1:i2)) (namer(i),i=1,kcen)
      endif
      if(index(line,'RESIDUE_POINTER').ne.0) then
      read(2,5) line
      i1=index(line,'(')
      i2=index(line,')')
      read(2,line(i1:i2)) (ncx(i),i=1,kcen)
      endif
      goto 10
c-----------------------------------------------------------read old top
      else
      read(2,*) kam,(nd(i),i=1,10),kcen
      read(2,*)
      read(2,*)
      read(2,'(20a4)') (mnx(i),i=1,kam)
      nlin1=kam/5
      if(nlin1*5.ne.kam) nlin1=nlin1+1
      nlin2=kam/12
      if(nlin2*12.ne.kam) nlin2=nlin2+1
      nd2=nd(1)*nd(1)
      nlin3=nd2/12
      if(nlin3*12.ne.nd2) nlin3=nlin3+1
      do i=1,2*(nlin1+nlin2)+nlin3
      read(2,*)
      enddo
      read(2,'(20a4)') (namer(i),i=1,kcen)
      read(2,'(12i6)') (ncx(i),i=1,kcen)
      endif
c---------------------------------------------remove H and compress mnam
100   i=0
      ncen(0)=0
      ncx(kcen+1)=kam+1
      do n=1,kcen
      do k=ncx(n),ncx(n+1)-1
      if(index(mnx(k),'H').eq.0) then
      i=i+1
      inx(k)=i
      mnam(i)=mnx(k)
         else
         inx(k)=0
         endif
      enddo
      ncen(n)=i
      enddo
      kto=kam
      kam=i
c------------------------------------------------------identify subunits
      do i=1,kcen
      do j=ncen(i-1)+1,ncen(i)
      nunit(j)=i
      munit(j)=namer(i)
      enddo
      enddo
c-----------------------------------------------------------cleanup data
150   do i=1,kam
20    if(munit(i)(:1).eq.' '
     1 .or.(munit(i)(:1).eq.'D'.and.munit(i)(2:2).ne.'M')
     1 .or.munit(i)(:1).eq.'R') then
      munit(i)=munit(i)(2:)
      goto 20
      endif
30    if(mnam(i)(:1).eq.' '.or.(ichar(mnam(i)(:1)).ge.48
     1 .and.ichar(mnam(i)(:1)).le.57)) then
      mnam(i)=mnam(i)(2:)
      goto 30
      endif
40    iq=index(mnam(i),'''')
      if(iq.ne.0) then
      mnam(i)(iq:iq)='*'
      goto 40
      endif
      enddo
      return
      end
