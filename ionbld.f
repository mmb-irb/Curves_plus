      subroutine ionbld
      include 'curves_data.inc'
      character*1 na,nt,base,type,first
      character*4 mnam,munit,snam,sunit,name,ban,nback,inam
      character*8 ibnam
      character*20 ilnam
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 nind,lpa,traj,start,need(20)
      dimension r1(4,3),r2(4,3),t(3,3)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
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
c------------------------------------------------------------locate ions
      npt=(nlev-1)*nspl+1
      do k=1,kisa
      in=ilis(k,2)
      posn=pari(k,1)
      rmin=pari(k,2)
      ang=pari(k,3)
      imin=1+nint((posn-1)*nspl)
      rx=uint(imin,3,1)
      ry=uint(imin,3,2)
      rz=uint(imin,3,3)
      ca=cos(ang)
      sa=sin(ang)
      xx=rmin*uint(imin,2,1)
      yy=rmin*uint(imin,2,2)
      zz=rmin*uint(imin,2,3)
      cors(k,1)=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     & (rx*rz*(1-ca)+ry*sa)*zz+uint(imin,4,1)
      cors(k,2)=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     & (ry*rz*(1-ca)-rx*sa)*zz+uint(imin,4,2)
      cors(k,3)=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     & (rz*rz+(1-rz*rz)*ca)*zz+uint(imin,4,3)
      snam(k)=inam(in)
      sunit(k)='ION'
      nunis(k)=1
      enddo
      kas=kisa
      khs=0
      return
      end
