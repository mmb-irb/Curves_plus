      subroutine setup(ions)
      include 'curves_data.inc'
      character*1 base,type,nke
      character*4 name,ban,nat,liga,bln,inam
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      character line*80,nbo*6,ibnam*8,ilnam*20
      logical*2 ions,there,start
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/ion/pari(n1,3),ilis(n1,2),klis(40),ilib(40),kion(5),
     1 kisa,nion,nspl,inam(40)
      common/lgd/ilnam(20),paxe(20),blef(20,40,3),parl(20,9),
     1 ilref(20),iluni(20,2),ilfrm(20,-3:20),bln(20,40),nlig,liga(20)
      common/str/nbo(40),nat(40,8),ish(40,8),ntor,nke(40)
      kfi=index(lib,' ')-1
c------------------------------------------------------------input bases
      inquire(file=lib(:kfi)//'_b.lib',exist=there)
      if(.not.there) then
      write(6,600) '  ---- Base lib file missing ----'
600   format(/2x,a,/)
      stop
      endif
      open(unit=4,file=lib(:kfi)//'_b.lib',status='old')
      k=0
10    k=k+1
      read(4,'(a)',end=100) line
         if(line(:1).eq.'#') then
         k=k-1
         goto 10
         endif
      base(k)=line(1:1)
      type(k)=line(3:3)
      read(line(4:),*) ibref(k),ibnam(k)
      do j=1,ibref(k)
      read(4,*) (bref(k,j,m),m=1,3),ban(k,j)
      enddo
      goto 10
100   nbas=k-1
      close(4)
c-----------------------------------------------------------input strand
      inquire(file=lib(:kfi)//'_s.lib',exist=there)
      if(.not.there) then
      write(6,600) '  ---- Strand lib file missing ----'
      stop
      endif
      open(unit=4,file=lib(:kfi)//'_s.lib',status='old')
      k=0
20    k=k+1
      read(4,'(a)',end=200) line
         if(line(:1).eq.'#') then
         k=k-1
         goto 20
         endif
      nke(k)=line(:1)
      if(line(:1).eq.'S') then
      nbo(k)='Phase'
      nke(k)='-'
      k=k+1
      nbo(k)='Ampli'
      nke(k)='-'
      k=k+1
      nbo(k)='Puckr'
      nke(k)='S'
      read(line(2:),*) (nat(k,j),j=1,5)
         else
         jup=4
         if(line(:1).eq.'B') jup=8
         read(line(2:),*) (nat(k,j),j=1,jup),nbo(k)
         do j=1,jup
         if(nat(k,j)(:1).eq.'+') then
         ish(k,j)=1
         nat(k,j)=nat(k,j)(2:)
            else if(nat(k,j)(:1).eq.'-') then
            ish(k,j)=-1
            nat(k,j)=nat(k,j)(2:)
               else
               ish(k,j)=0
               endif
         enddo
         endif
      goto 20
200   ntor=k-1
      close(4)
c-------------------------------------------------------------input ions
      if(ions) then
      inquire(file=lib(:kfi)//'_i.lib',exist=there)
      if(.not.there) then
      write(6,600) '  ---- Ion lib file missing ----'
      stop
      endif
      open(unit=4,file=lib(:kfi)//'_i.lib',status='old')
      k=0
40    k=k+1
      read(4,'(a)',end=45) line
         if(line(:1).eq.'#') then
         k=k-1
         goto 40
         endif
      read(line,*) inam(k),ilib(k)
      goto 40
45    nion=k-1
      close(4)
      endif
c-----------------------------------------------------------input ligand
      nlig=0
      if(lig.ne.' ') then
      kfi=index(lig,' ')-1
      inquire(file=lig(:kfi)//'.lig',exist=there)
      if(.not.there) then
      write(6,600) '  ---- Ligand lib file missing ----'
      stop
      endif
      open(unit=4,file=lig(:kfi)//'.lig',status='old')
      k=0
30    k=k+1
      read(4,'(a)',end=300) line
         if(line(:1).eq.'#') then
         k=k-1
         goto 30
         endif
      knam=index(line,' ')-1
      liga(k)=line(1:knam)
      read(line(knam+1:),*) ilref(k),iluni(k,1),iluni(k,2),
     1 paxe(k),ilnam(k)
      read(4,'(a)') line
         m=0
         n=0
         start=.false.
6        m=m+1
         if(.not.start.and.line(m:m).ne.' ') then
         start=.true.
         ms=m
         endif
         if(start.and.line(m:m).eq.' ') then
         start=.false.
         n=n+1
         read(line(ms:m-1),*) ilfrm(k,n)
         endif
         if(m.lt.80) goto 6
         ilfrm(k,0)=n
      read(4,*) (ilfrm(k,j),j=-3,-1)
      do j=1,ilref(k)
      read(4,*) (blef(k,j,m),m=1,3),bln(k,j)
      enddo
      goto 30
300   nlig=k-1
      close(4)
      endif
      return
      end
