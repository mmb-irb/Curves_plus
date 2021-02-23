      subroutine axis
      include 'curves_data.inc'
      character*1 na,nt
      character*4 nback
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,lpa
      dimension r1(4,3),r2(4,3)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
         do i=1,nlev
         npl(i)=0
         enddo
c--------------------------------------use upm for aver p point and dyad
      iup=nlev
      if(circ) iup=nlev+isym
      do i=1+isym,iup !...level
      iu=i
      il=iu-isym
      if(i.gt.nlev) iu=i-nlev
         do k=1,nst
         lpa(il,k)=.false.
         enddo
      do k=1,nst !...strand
      if(ni(il,k).gt.0.and.ni(iu,k).gt.0) then
      lpa(il,k)=.true.
      do m=1,4
      do j=1,3
      r1(m,j)=ref(il,k,m,j)
      r2(m,j)=ref(iu,k,m,j)
      enddo
      enddo
      call screw(r1,r2,-1) !...screw axis from system il to iu
      dx=upm(il,4,1)-usc(4)
      dy=upm(il,4,2)-usc(5)
      dz=upm(il,4,3)-usc(6)
      dotl=dx*usc(1)+dy*usc(2)+dz*usc(3)
         dx=upm(iu,4,1)-usc(4)
         dy=upm(iu,4,2)-usc(5)
         dz=upm(iu,4,3)-usc(6)
         dotu=dx*usc(1)+dy*usc(2)+dz*usc(3)
      iln=npl(il)+1
      iun=npl(iu)+1
      npl(il)=iln
      npl(iu)=iun
 
crl   copy screw axis to helical axis vector upl 1->3 for lower and upper level
crl   copy axis ref point to upl 4->6, so base ref to point is perp. to helical
 
      do j=1,3
      upl(il,iln,j)=usc(j)
      upl(il,iln,j+3)=usc(j+3)+dotl*usc(j)
      upl(iu,iun,j)=usc(j)
      upl(iu,iun,j+3)=usc(j+3)+dotu*usc(j)
      enddo
c      if(test) write(6,20) il,iu,k,(upl(il,iln,j),j=1,6),var(6)
c20    format(2x,i3,'-',i3,' / ',i1,' U: ',3f8.2,' P: ',3f8.2,' T ',f8.2)
      endif
      enddo
      enddo
      return
      end
