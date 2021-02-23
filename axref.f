      subroutine axref
      include 'curves_data.inc'
      character*1 na,nt
      character*4 snam,sunit,inam
      character*128 file,ftop,lis,lib,lig,ibld,sol,back
      logical*2 circ,line,zaxe,fit,test,ions,refo,lpa,traj,axfrm,
     1 frames,wat
      dimension f(4,3),pl(3),pu(3),ul(3),uu(3),dya(3),
     1 so(3),sn(3),sto(3)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,wat
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/ion/pari(n1,3),ilis(n1,2),klis(40),ilib(40),kion(5),
     1 kisa,nion,nspl,inam(40)
      common/lam/cors(n1,3),snam(n1),sunit(n1),nunis(n1),
     1 mats(3*n1),kas,khs,kces
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/sto/bpdat(n3,n5),traj
      if(test) then
c--------------------------------------------------------step screw axes
      r=0.5d0
      m=0
      lm=0
      do i=1,nlev
      do k=1,npl(i)
      do j=1,3
      cors(m+1,j)=upl(i,k,j+3)-r*upl(i,k,j)
      cors(m+2,j)=upl(i,k,j+3)+r*upl(i,k,j)
      enddo
         do j=1,2
         m=m+1
         lm=lm+1
         mats(lm)=m
         snam(m)='H'
         nunis(m)=i
         sunit(m)='Scrw'
         enddo
         lm=lm+1
         mats(lm)=80000
      enddo
      enddo
      kas=m
      if(.not.traj.or.(itst.eq.itnd)) call pdbout('_s',2)
c---------------------------------------------------------averaged axes
      r=2.0d0
      m=0
      lm=0
      do i=1,nlev
      do j=1,3
      cors(m+1,j)=uvw(i,4,j)-r*uvw(i,3,j)
      cors(m+2,j)=uvw(i,4,j)+r*uvw(i,3,j)
      enddo
         do j=1,2
         m=m+1
         lm=lm+1
         mats(lm)=m
         snam(m)='H'
         nunis(m)=i
         sunit(m)='Aver'
         enddo
         lm=lm+1
         mats(lm)=80000
      enddo
      kas=m
      khs=lm-1
      if(.not.traj.or.(itst.eq.itnd)) call pdbout('_v',2)
      endif
c--------------------------------------------------------------draw axis
      call intaxe
      return
      end
