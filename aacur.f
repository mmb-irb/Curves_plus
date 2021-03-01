      program aacur
c***********************************************************************
c                                                                      *
c     -------------------- Curves+ (3.0nc)  --------------------       *
c                                                                      *
c     Helical analysis of irregular nucleic acids                      *
c                                                                      *
c     Cambridge helical parameter definitions                          *
c     and Tsukuba base reference conventions                           *
c                                                                      *
c     Simplified output - no conflict between local and global params. *
c                                                                      *
c     Up to 4 strands of unequal lengths with interruptions            *
c                                                                      *
c     New helical axis definition without minimization                 *
c                                                                      *
c     Improved groove geometry measurements                            *
c                                                                      *
c     Authors:  R.Lavery, K. Zakrzewska                                *
c     (Screw axes: J. Maddocks, M. Moakher, D. Petkeviciute)           *
c     (Helicoidal volumes for ion densities: J.H. Maddocks)            *
c     (Bisection algorithm for ion posn: M. Pasi, J.H. Maddocks)       *
c     (NetCDF, XTC input: M. Pasi)                                     *
c                                                                      *
c                                                                      *
c     Ver 1.1 Change allow for unit numbers >999                 10/09 *
c     Ver 1.2 Tidy code, clip paths in O/P, auto remove water    01/11 *
c             Change _b and _x O/P files to suit Chimiera              *
c             Correct graphic groove O/P                               *
c     Ver 1.3 Automatically eliminate protein, water and hetero  04/11 *
c             at imput. Automatically choose first copy of any         *
c             atoms which occur as multiple copies within a given      *
c             residue. Recognize bases with lower case names.          *
c             Allow backbones splines based on atoms other than P      *
c     Ver 2.0 Added analysis of ligand positions +               10/11 *
c             Tot rise/twist + Optional frame O/P + New cda format     *
c     Ver 2.1 Added analysis of ion positions                    12/11 *
c     Ver 2.2 Corrected problem with line=.t.                    03/12 *
c     Ver 2.3 Corrected problem with circ=.t. for ions           03/12 *
c     Ver 2.4 Added ion identity for rmsf                        04/12 *
c     Ver 2.5 Optional solvent I/P and topology from pdb file    10/13 *
c     Ver 2.6 New ion location algorithm (speed and accuracy)    03/14 *
c     Ver 2.7 Analyze multiple segments of an MD trajectory      05/15 *
c     Ver 2.8nc Amber NetCDF, Gromacs XTC trajectory input       11/15 *
c     Ver 2.9nc Analyze curvature (radius and register)          04/16 *
c     Ver 3.0nc Simplified namelist and PDB input                10/16 *
c***********************************************************************
      include 'curves_data.inc'
#ifdef NETCDF
      include 'netcdf.inc'
#endif
#ifdef XTC
      integer*4 xtc_handle, xtc_status, xtc_natom
      real*4 xtc_cort, xtc_box
      dimension xtc_cort(3,n1),xtc_box(3)
#endif
      character*1 base,type,na,nt,nke
      character*4 ban,mnam,munit,snam,sunit,nat,nback,
     & inam,name
      character*20 ilnam,string*9
      character*80 lini,cline
      character*128 file,ftop,lis,lib,ibld,sol,back
      character ibnam*8,nbo*6,dir(-1:1)*5,ext*3,eto*3
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 lpa,start,traj,box,there,skip
      integer*4 inx(n1p),intdat(1000),ihold(4),nhold(4)
      dimension cort(2000,3),ibrack(200,2)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/bisect/bsx(3)
      common/cha/file,ftop,lis,lib,ibld,sol,back
      common/dat/wback,wbase,rvfac,isym,itst,itnd,itdel,itbkt,
     1 naxlim,circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),upm(n3,4,3),nback(4)
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
      common/str/nbo(40),nat(40,8),ish(40,8),ntor,nke(40)
      data dir/'3''-5''',' ','5''-3'''/
      traj=.false.
c----------------------------------------------------------namelist data
      file=' '
      ftop=' '
      lis=' '
      lib=' '
      ibld=' '
      back='P'
      sol=' '
 
      wback=2.9
      wbase=3.5
      rvfac=7.5
 
      isym=1
      itst=0
      itnd=0
      itdel=1
      itbkt=0
      naxlim=0
 
      circ=.false.
      line=.false.
      zaxe=.false.
      fit=.true.
      test=.false.
      ions=.false.
      refo=.false.
      frames=.false.
      axfrm=.false.
 
      nspl=3
      itcount=0
#ifdef XTC
      xtc_handle = -1
      xtc_status = -1
#endif
c---------------------------------------------------------input run data
      call nml
         if(itbkt.gt.0) then
         do i=1,itbkt
         read(5,*) ibrack(i,1),ibrack(i,2)
         enddo
         endif
      if(zaxe.and.line) then
      write(6,600) '  ---- Choose between ZAXE and LINE ----'
600   format(/,2x,a,/)
      stop
      endif
      if(zaxe.and.circ) then
      write(6,600) '  ---- choose between ZAXE and CIRC ----'
      stop
      endif
      if(line.and.circ) then
      write(6,600) '  ---- choose between LINE and CIRC ----'
      stop
      endif
      kfi=index(lis,' ')-1
      open(unit=3,file=lis(:kfi)//'.cda',status='new')
      if(ions) open(unit=7,file=lis(:kfi)//'.cdi',status='new')
      if(frames) open(unit=8,file=lis(:kfi)//'.fra',status='new')
         ib=index(back,' ')-1
         if(ib.gt.0) then
         if(index(back,'/').eq.0) then
         do i=1,4
         nback(i)=back(:ib)
         enddo
            else
            i=1
            jl=1
            do j=1,ib+1
            if(back(j:j).eq.'/'.or.j.eq.ib+1) then
            nback(i)=back(jl:j-1)
            i=i+1
            jl=j+1
            endif
            enddo
            endif
         endif
c--------------------------------------------optional old base reference
      if(refo) then
      th1=132.193d0
      th2=-54.512d0
      dis=4.5033d0
         else
         th1=141.47d0
         th2=-54.41d0
         dis=4.7024d0
      endif
c-------------------------------------------------parse nucleotide input
         do k=1,4
         do i=1,n3
         ni(i,k)=0
         enddo
         enddo
      read(5,*) nst,(idr(j),j=1,4)
      if(nst.eq.0.or.nst.gt.4) then
      write(6,600) ' ---- Choose 1 -> 4 strands ----'
      stop
      endif
      nlev=0
      do k=1,nst
      idr(k)=sign(1,idr(k))
      read(5,5) lini
5     format(a)
      do i=1,80
      ich=ichar(lini(i:i))
      if((ich.lt.48.and.ich.ne.13.and.ich.ne.32.and.ich.ne.45)
     1 .or.ich.gt.58) then
      write(6,600) '---- Strand input should be a number or '':'' ----'
      stop
      endif
      enddo
      m=0
      n=0
      start=.false.
      jump=0
      jmul=0
6     m=m+1
         if(lini(m:m).eq.':') then
         lini(m:m)=' '
         jump=-1
         endif
            if(lini(m:m).eq.'*') then
            lini(m:m)=' '
            jmul=-1
            endif
      if(.not.start.and.lini(m:m).ne.' ') then
      start=.true.
      ms=m
      endif
      if(start.and.lini(m:m).eq.' ') then
      start=.false.
      read(lini(ms:m-1),*) numb
 
      if(jmul.eq.-1) then
      if(ni(n,k).ne.0) then
      write(6,600) '  ---- * only after 0 ----'
      stop
      endif
      n=n+numb-1
      if(n.gt.n3) then
      write(6,600) '  ---- Too many levels ----'
      stop
      endif
      jmul=0
      if(m.lt.80) goto 6
      endif
 
      n=n+1
      if(n.gt.n3) then
      write(6,600) '  ---- Too many levels ----'
      stop
      endif
      ni(n,k)=numb
      if(jump.eq.-1) jump=n+1
         if(jump.eq.n) then
         n=n-1
         if(numb.le.0.or.ni(n,k).le.0) then
         write(6,600) '  ---- '':'' only between positive numbers ----'
         stop
         endif
         jd=sign(1,numb-ni(n,k))
         js=ni(n,k)+jd
         do j=js,numb,jd
         n=n+1
         if(n.gt.n3) then
         write(6,600) '  ---- Too many levels ----'
         stop
         endif
         ni(n,k)=j
         enddo
         jump=0
         endif
      endif
 
      if(m.lt.80) goto 6
      if(n.gt.nlev) nlev=n
      enddo
      if(isym.le.0.or.isym.gt.nlev) then
      write(6,600) ' ---- Invalid isym value ----'
      stop
      endif
c-----------------------check and if necessary permute strands or levels
      do i=1,nlev
      start=.false.
      do k=1,nst
      if(ni(i,k).gt.0) start=.true.
      enddo
      if(.not.start) then
      write(6,600) '---- Must have at least one base per level ----'
      stop
      endif
      enddo
 
      if(idr(1).lt.0) then
      istr=0
      do k=2,nst
      if(idr(k).gt.0.and.istr.eq.0) istr=k
      enddo
      if(istr.gt.0) then !.....cyclically permute
      write(6,21)
21    format(2x,'... cyclic permutation to make first strand 5''->3''')
      do k=1,nst
      ihold(k)=idr(k)
      enddo
      do i=1,nlev
      do k=1,nst
      nhold(k)=ni(i,k)
      enddo
      do k=1,nst
      kup=k+istr-1
      if(kup.gt.nst) kup=kup-nst
      ni(i,k)=nhold(kup)
      if(i.eq.1) idr(k)=ihold(kup)
      enddo
      enddo
         else !.....invert levels
         write(6,22)
22       format(2x,'... levels inverted to make first strand 5''->3''')
         do k=1,nst
         idr(k)=-idr(k)
         enddo
         do k=1,nst
         do i=1,int(nlev/2)
         inv=nlev+1-i
         nkeep=ni(i,k)
         ni(i,k)=ni(inv,k)
         ni(inv,k)=nkeep
         enddo
         enddo
         endif
      endif
c--------------------------------------------reference data for analysis
      call setup(ions)
c------------------------------------------------------------------input
      kfi=index(file,' ')-1
4     inquire(file=file(:kfi),exist=there)
      if(.not.there) then
      write(6,600) '  ---- Input file missing ----'
      stop
      endif
      open(unit=1,file=file(:kfi),status='old')
      ext=file(index(file,'.',.true.)+1:)
      if(ext.eq.'trj'.or.ext.eq.'nc'.or.
     1   ext.eq.'xtc'.or.ext.eq.'pgp') then
         if(itst.gt.0.and.itnd.eq.0) itnd=itst
         if(itnd.gt.0.and.itst.eq.0) itst=1
           if(itst.eq.0.and.itnd.eq.0) then
           itst=1
           itnd=10000000
           endif
             if(itst.ne.itnd) then
             traj=.true.
             test=.false.
             if(axfrm) write(6,*)
     1       '  ---- No axis frames with trajectories ----'
             axfrm=.false.
             endif
         endif
c----------------------------------------------------------read topology
         if(ext.eq.'trj'.or.ext.eq.'nc'.or.ext.eq.'xtc') then
         kft=index(ftop,' ')-1
         inquire(file=ftop(:kft),exist=there)
         if(.not.there) then
         write(6,600) '  ---- Top file missing ----'
         stop
         endif
         open(unit=2,file=ftop(:kft),status='old')
         eto=ftop(kft-2:)
         call intop(eto,kto,inx)
         close(2)
         endif
         if(ext.eq.'trj') then
         nat3=3*kto
         nlin=nat3/10
         if(nlin*10.lt.nat3) nlin=nlin+1
         read(1,*)
         do i=1,nlin
         read(1,*)
         enddo
         read(1,5) cline
         box=.false.
            if(cline(77:77).ne.'.') then
            box=.true.
            nlin=nlin+1
            endif
         close(1)
         if(box) write(6,5) '  ... TRJ includes box dimensions'
         open(unit=1,file=file(:kfi),status='old')
         read(1,*)
#ifdef NETCDF
         elseif(ext.eq.'nc') then
         ncstatus = nf_open(file(:kfi), nf_nowrite, ncid)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
c--- check consistency: kto and itnd
         ncstatus = nf_inq_dimid(ncid, 'atom', nccid)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
         ncstatus = nf_inq_dimlen(ncid, nccid, nckto)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
            if(kto.gt.nckto) then ! accept partial topology
            write(6,600) '  ---- Top / NetCDF molecule mismatch ----'
            stop
            endif
         ncstatus = nf_inq_dimid(ncid, 'frame', nccid)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
         ncstatus = nf_inq_dimlen(ncid, nccid, ncitnd)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
            if(itnd.gt.ncitnd) then ! consistent with trj behaviour
            itnd = ncitnd
            endif
c--- get coordinates
         ncstatus = nf_inq_varid(ncid, 'coordinates', nccid)
         if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
#endif
#ifdef XTC
         elseif(ext.eq.'xtc') then
         call f77_molfile_init
         call f77_molfile_open_read(
     1     xtc_handle, xtc_natom, file(:kfi), 'auto')
         if (xtc_handle.lt.0) call xtcerror(
     1    'Cannot open trajectory file. Fatal.     ', 1)
c--- check consistency: kto
         if(kto.ne.xtc_natom) then ! only accept exact match
            write(6,600) '  ---- Top / XTC molecule mismatch ----'
            stop
            endif
#endif
        endif
c-----------------------------------------------------if traj begin loop
      itread=0

      itact=1
      itj=itst
400   if(ext.eq.'trj') then
         if(itj.gt.itact) then
         print *,' SKIPPING ',itj-itact,' itj= ',itj,' itact= ',itact
         nsk=(itj-itact)*nlin
         do i=1,nsk
         read(1,*)
         enddo
         itact=itj
         endif
      ilow=1
      ihig=min(2000,kto)
14    read(1,15,end=70) (cort(i,1),cort(i,2),cort(i,3),i=1,ihig-ilow+1)
15    format(10f8.3)
      do i=ilow,ihig
      j=inx(i)
         if(j.gt.0) then
         do k=1,3
         corm(j,k)=cort(i-ilow+1,k)
         enddo
         endif
      enddo     
      ilow=ilow+2000
      ihig=min(ihig+2000,kto)
      if(ilow.le.kto) goto 14

      if(box) read(1,*)
      itread=itread+1
      elseif(ext.eq.'nc') then
#ifdef NETCDF
      ilow=1
      ihig=min(2000,kto)
17    ncstatus = nf_get_varm_double(ncid, nccid, 
     1        (/ 1, ilow, itj /),
     1        (/ 3, ihig-ilow+1, 1 /),
     1        (/ 1,  1,  1 /),
     1        (/ 2000, 1, 6000 /), cort)
      if(ncstatus.ne.nf_noerr) call ncerror(ncstatus, 1)
      do i=ilow,ihig
      j=inx(i)
         if(j.gt.0) then
         do k=1,3
         corm(j,k)=cort(i-ilow+1,k)
         enddo
         endif
      enddo     
      ilow=ilow+2000
      ihig=min(ihig+2000,kto)
      if(ilow.le.kto) goto 17
      itread=itread+1
#else
      write(6,600) 'To read Amber netCDF files, please compile Curves+ w&
     &ith netCDF support (see manual).'
      stop
#endif
      elseif(ext.eq.'xtc') then
#ifdef XTC
      nsk=itj-itact
         do i=1,nsk
         xtc_status=0
         call f77_molfile_read_next( ! skip read
     1     xtc_handle, xtc_natom, corm, xtc_box, xtc_status)
         if (xtc_status.eq.0) call xtcerror(
     1    'Cannot skip next snapshot. Fatal.       ', 1)
         enddo
         itact=itj
      xtc_status=1
 18   call f77_molfile_read_next( ! real read
     1 xtc_handle, xtc_natom, xtc_cort, xtc_box, xtc_status)
      if (xtc_status.eq.0) call xtcerror(
     1    'Fatal error reading from XTC file.      ', 1)
      if (xtc_status.eq.2) goto 112 ! EOF
      do i=1,xtc_natom
      j=inx(i)
         if(j.gt.0) then
         do k=1,3
         corm(j,k)=xtc_cort(k,i)
         enddo
         endif
      enddo     
      itread=itread+1
#else
      write(6,600) 'To read Gromacs XTC files, please compile Curves+ wi&
     &th XTC support (see manual).'
      stop
#endif
         else
         call input(ext)
         if(kam.eq.0) goto 112
         itread=itread+1
         endif
      call locate(itj,itcount,skip)
      if(skip) goto 111
      if(itj.gt.itst.or.itcount.gt.0) goto 16
c--------------------------------------------------------------find ions
      if(ions) then
      ki=0
      do i=1,5
      kion(i)=0
      enddo
      do j=1,nion
      klis(j)=0
      enddo
      do i=1,kam
      name=mnam(i)
         do j=1,nion
         if(name.eq.inam(j).or.inam(j).eq.'*') then
         ki=ki+1
         ichg=ilib(j)
         ilis(ki,1)=i
         ilis(ki,2)=j
            if(ichg.eq.1) then
            kion(2)=kion(2)+1
            else if(ichg.eq.2) then
            kion(3)=kion(3)+1
            else if(ichg.eq.-1) then
            kion(4)=kion(4)+1
            else if(ichg.eq.-2) then
            kion(5)=kion(5)+1
            endif
         endif
         enddo
      enddo
      if(ki.gt.n3*37) stop '  ---- Problem kion > intdat size ----'
      kion(1)=ki
      endif
c------------------------------------------------------input and summary
      write(6,20) nst,kam,kcen
20    format(/2x,'Strands = ',i4,' Atoms = ',i5,' Units = ',i5)
      write(6,30) nlev
30    format(/2x,'Combined strands have ',i4,' levels ...'/)
      do k=1,nst
         num=0
         do i=1,nlev
         if(ni(i,k).ne.0) num=num+1
         enddo
      imin=1
      imax=min(50,nlev)
      write(6,40) k,num,dir(idr(k)),(na(i,k),i=imin,imax)
40    format(2x,'Strand ',i2,' has ',i3,' bases (',a5,'): ',50a1)
38       if(nlev.gt.imax) then
         imin=imax+1
         imax=min(imin+49,nlev)
         write(6,42) (na(i,k),i=imin,imax)
42       format(35x,50a1)
         goto 38
         endif
      enddo
      if(circ) write(6,60)
60    format(/2x,'Ring closure requested')
c---------------------------------------------------------excluded bases
      start=.true.
      do k=1,nst
          do i=1,nlev
          if(ni(i,k).lt.0) then
            if(start) then
            write(6,*)
            start=.false.
            endif
          write(6,62) na(i,k),i,k
62        format(2x,'Base ',a,i3,' in strand ',i1,
     1    ' excluded from axis calculation')
          endif
          enddo
      enddo
      call flush(6)
c-------------------------------------------------------------parameters
16       do i=1,nlev
         do j=1,n5
         bpdat(i,j)=300.
         enddo
         enddo

c        if(mod(itj,10000).eq.0) then
c         write(6,63) itj
c63       format(/2x,'Read ',i7,' snapshots')
c         call flush(6)
c         endif

      call params(totbend)
c==========================================================O/P flat file
         k=0
         jtot=21+nst*(ntor-1)+(nst-1)*4
         do i=1,nlev
         do j=1,jtot
         k=k+1
         intdat(k)=nint(bpdat(i,j)*100)
            if(k.eq.1000) then
            write(3,72) (intdat(m),m=1,1000)
            k=0
            endif
         enddo
         enddo
         k=k+1
         intdat(k)=nint(totbend*100)
         write(3,72) (intdat(i),i=1,k)
72       format(200i6)
c-------------------------------------------------------------------ions
         if(kion(1).gt.0) then
            if(itj.eq.itst) then
            write(7,73) nlev,nst,nion
            write(7,74) (na(i,1),i=1,nlev)
            write(7,73) (ilib(i),i=1,nion)
            write(7,75) (inam(i),i=1,nion)
            write(7,73)
73          format(200i7)
74          format(100a1)
75          format(40a4)
            endif
         write(7,73) kisa
         k=0
         do i=1,kisa
         do j=1,3
         k=k+1
         intdat(k)=nint(pari(i,j)*1000)
            if(k.eq.1000) then
            write(7,73) (intdat(m),m=1,1000)
            k=0
            endif
         enddo
         k=k+1
         intdat(k)=iunit(i)
            if(k.eq.1000) then
            write(7,73) (intdat(m),m=1,1000)
            k=0
            endif
         enddo
         if(k.gt.0) write(7,73) (intdat(i),i=1,k)
         endif
c-----------------------------------------------------------------------
111   itact=itact+1
      itj=itj+itdel
      if(itj.le.itnd) goto 400
         if(itcount.lt.itbkt) then
         itcount=itcount+1         
         itst=ibrack(itcount,1)
         itnd=ibrack(itcount,2)
         write(6,116) itcount,itst,itnd
116      format(/2x,'Bracket ',i3,' ITST= ',i7,' ITND= ',i7)
         itj=itst
         goto 400
         endif
112   close(3)
      if(ions) close(7)
      if(frames) close(8)
70    if(itread.gt.1) write(6,79) ext,itread
#ifdef NETCDF
      ncstatus = nf_close(ncid)
#endif
#ifdef XTC
      if(xtc_handle.ge.0) call f77_molfile_close_read(
     1 xtc_handle, xtc_status)
#endif
79    format(/2x,'... ',a3,' loop read ',i8,' snapshots'/)
c------------------------------------------------------------axis frames
      if(axfrm) then
      kfi=index(lis,' ')-1
      open(unit=8,file=lis(:kfi)//'.afr',status='new')
      write(8,'(i4)') nlev
      write(8,114) (((uvw(i,j,k),k=1,3),j=1,4),i=1,nlev)
114   format(12f12.7)
         if(.not.circ) then
         do j=20,21
         bpdat(1,j)=0.
         bpdat(nlev-2,j)=0.
         bpdat(nlev-1,j)=0.
         enddo
         endif
         write(8,115) (bpdat(i,20),bpdat(i,21),i=1,nlev)
115      format(2f12.5)
      close(8)
      endif
c--------------------------------------------------------------ion build
      if(ibld.ne.' ') then
      kfi=index(ibld,' ')-1
      open(unit=10,file=ibld(:kfi)//'.cdi',status='old')
      read(10,73) nlev,nst,nion
      read(10,74) (na(i,1),i=1,nlev)
      read(10,73) (ilib(i),i=1,nion)
      read(10,75) (inam(i),i=1,nion)
      read(10,73) (klis(i),i=1,nion)
      read(10,73) kisa
      read(10,73) (intdat(i),i=1,kisa*4)
         k=0
         do i=1,kisa
         do j=1,3
         k=k+1
         pari(i,j)=float(intdat(k))/1000
         enddo
         k=k+1
         ilis(i,2)=mod(intdat(k),100)
         ilis(i,1)=(intdat(k)-ilis(i,2))/100
         enddo
      close(10)
      call ionbld
      call pdbout('_C',0)
      endif
c-----------------------------------------------------------------------
      end
