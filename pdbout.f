      subroutine pdbout(qual,icon)
      include 'curves_data.inc'
      character*1 na,nt
      character*2 qual
      character*4 snam,sunit
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      dimension igx(n1)
c------------------------------------------------resort grv for duplexes
          igrv=0
          do i=1,kas
          igx(i)=i
             if(qual.eq.'_b'.and.nunis(i).eq.nst+1) then
             igrv=i
             goto 12
             endif
          enddo
12    k=igrv-1
      if(qual.eq.'_b'.and.nst.ge.2.and.igrv.gt.0) then
      do i=nst+1,2*nst
         do j=igrv,kas,2
         if(nunis(j).eq.i) then
         k=k+1
         igx(k)=j
         k=k+1
         igx(k)=j+1
         endif
         enddo
       enddo
      endif
c--------------------------------------------------------------------pdb
14    kfi=index(lis,' ')-1
      open(unit=2,file=lis(:kfi)//qual//'.pdb',status='new')
      write(2,5) 'HEADER    '//lis(:kfi)//' from Curves'
5     format(a)
      do i=1,kas
      ix=igx(i)
      write(2,10) 'HETATM',i,snam(ix),sunit(ix)(:3),sunit(ix)(4:4),
     1 nunis(ix),cors(ix,1),cors(ix,2),cors(ix,3)
10    format(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
      enddo
c----------------------------------------------------------------connect
      if(icon.gt.0) then
      do il=1,khs-1
      iu=il+1
      mil=mats(il)
      miu=mats(iu)
      if(mil.ne.80000.and.miu.ne.80000) write(2,15) mil,miu
15    format('CONECT',2i5)
      enddo
      endif
      close(2)
      return
      end
