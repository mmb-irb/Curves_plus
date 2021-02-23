      subroutine input(ext,kto)
      include 'curves_data.inc'
      character*1 base,type,aloc,acha
      character*3 ext,auni,p3(20),m3(20),atran(20,2)
      character*4 mnam,munit,name,anam,ban
      character*8 ibnam
      character*80 lini,mcode
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      data p3/
     1 'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
     1 'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'/
      data m3/
     1 'ala','cys','asp','glu','phe','gly','his','ile','lys','leu',
     1 'met','asn','pro','gln','arg','ser','thr','val','trp','tyr'/
c--------------------------------------------------------read mac format
      if(ext.eq.'mac'.or.ext.eq.'MAC') then
      read(1,5) mcode
9     read(1,5) lini
      if(lini(:1).eq.'#') goto 9
      read(lini,6) kam,kcen
      i=0
      do k=1,kam
      read(1,5) lini
      if(index(lini(1:4),'H').eq.0) then
      i=i+1
      read(lini,7) mnam(i),corm(i,1),corm(i,2),corm(i,3),
     1 munit(i),nunit(i)
      endif
      enddo
      kam=i
5     format(a)
6     format(2i5)
7     format(a4,3f10.5,15x,a4,i4)
8     format(7x,6i5,8x,6i5)
c--------------------------------------------------------read pdb format
      else if(ext.eq.'pdb'.or.ext.eq.'PDB'.or.ext.eq.'pgp') then
      i=0
      kto=0
      ktr=0
100   read(1,5,end=110) lini
      if(index(lini,'ENDMDL').ne.0) goto 110
      if(lini(:6).eq.'TRANSL') then
      ktr=ktr+1
      inx=index(lini(8:),' ')+7
      iny=index(lini(inx+1:),' ')+inx
      atran(ktr,1)=lini(8:inx-1)
      atran(ktr,2)=lini(inx+1:iny-1)
      goto 100
      endif
      if(lini(:4).ne.'ATOM'.and.lini(:6).ne.'HETATM') goto 100
      read(lini,106) anam,aloc,auni,acha,nuni,x,y,z
106   format(12x,a4,a1,a3,1x,a1,i4,4x,3f8.3)
         if(index(anam,'H').ne.0) then
         kto=kto+1
         goto 100
         endif
      if(lini(:6).eq.'HETATM'.and.auni.ne.sol(:3)) goto 100
      i=i+1
      mnam(i)=  anam
107      if(auni(:1).eq.' ') then
         auni=auni(2:)
         goto 107
         endif
         if(auni(:1).eq.'D') then
         do m=1,nbas
         if(auni(2:2).eq.base(m)) then
         auni=auni(2:)
         goto 108
         endif
         enddo
         endif
108      do k=1,ktr
         if(auni.eq.atran(k,1)) auni=atran(k,2)
         enddo
         if(auni.ne.sol(:3).and.
     1   (auni.eq.'WAT'.or.auni.eq.'H2O'.or.auni.eq.'HOH')) then
         i=i-1
         goto 100
         endif
         do k=1,20
         if(auni.eq.p3(k).or.auni.eq.m3(k)) then
         i=i-1
         goto 100
         endif
         enddo
      munit(i)= auni//acha
      nunit(i)= nuni
      corm(i,1)=x
      corm(i,2)=y
      corm(i,3)=z
      goto 100
110   kam=i
      kto=kto+kam
      if(kam.eq.0) return
c-------------------------------------------------------------------trap
      else
      write(6,*) '  ---- Unknown geometry file type ----'
      stop
      endif
      if(ext.ne.'pgp') close(1)
c---------------left justify atom and subunit names and replace ! with *
      do i=1,kam
10    if(mnam(i)(:1).eq.' '.or.(ichar(mnam(i)(:1)).ge.48
     1 .and.ichar(mnam(i)(:1)).le.57)) then
      mnam(i)=mnam(i)(2:)
      goto 10
      endif
11    iq=index(mnam(i),'''')
      if(iq.ne.0) then
      mnam(i)(iq:iq)='*'
      goto 11
      endif
12    if(munit(i)(:1).eq.' ') then
      munit(i)=munit(i)(2:)
      goto 12
      endif
      enddo
c----------------------------------------------------------find subunits
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
      return
      end
