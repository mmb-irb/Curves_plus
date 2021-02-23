      subroutine ionpar
      include 'curves_data.inc'
      character*1 na,nt,base,type,first
      character*4 mnam,munit,snam,sunit,name,ban,nback,inam
      character*8 ibnam
      character*20 ilnam
      logical*2 sneg(n3),need(20),circ,line,zaxe,fit,test,ions,refo,
     1 axfrm,frames,nind,lpa,traj,start
      parameter (shell=30.d0)
      dimension rloc(n3*50),r1(4,3),r2(4,3),t(3,3),iloc(n3*50)
      common/bas/ibnam(20),bref(20,15,3),th1,th2,dis,ibref(20),
     1 ban(20,15),nbas,base(20),type(20)
      common/bisect/bsx(3)
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
      kaso=kas
      if(.not.traj) then
      write(6,64) (kion(i),i=1,5)
64    format(/2x,'(I) ',i3,' ions input (',i3,'+  ',i3,'++ ',
     1 i3,'- ',i3,'-- )',
     1 //2x,'Located',6x,'I-Axe',3x,'I-Dis',3x,'I-Ang'/)
      endif
c------------------------------------------------------------locate ions
      ks=0
      do k=1,kion(1)
      in=ilis(k,1)
      x0=corm(in,1)
      y0=corm(in,2)
      z0=corm(in,3)
      rmin=1.d5
      do i=1,nlev
      sneg(i)=.false.
      dx=x0-uvw(i,4,1)
      dy=y0-uvw(i,4,2)
      dz=z0-uvw(i,4,3)
      dot=dx*uvw(i,3,1)+dy*uvw(i,3,2)+dz*uvw(i,3,3)
      if(dot.lt.0) sneg(i)=.true.
      enddo
c--------------------------------------find perpendicular vector to axis
         iminim=0
         rglo=1.d4
         iglo=-1
         do il=1,nlev-1
         if(sneg(il).neqv.sneg(il+1)) then
         imin=-1
         ilow=1+(il-1)*nspl
         ihig=il*nspl+1
         iminim=iminim+1
         bsx(1)=x0
         bsx(2)=y0
         bsx(3)=z0
         call bisection(ilow,ihig,imin)
         dx=bsx(1)-uint(imin,4,1)
         dy=bsx(2)-uint(imin,4,2)
         dz=bsx(3)-uint(imin,4,3)
         rmin=sqrt(dx*dx+dy*dy+dz*dz)
            if(imin.ge.0.and.rmin.lt.rglo) then
            rglo=rmin
            iglo=imin
            endif
         endif
         enddo
         if(rglo.gt.shell) goto 21
c------------------------------------------------------------located ion
c     if(imin.eq.1.or.imin.eq.npt) goto 21
      ks=ks+1
      jt=ilis(k,2)
      iunit(ks)=nunit(ilis(k,1))*100+ilis(k,2)
      klis(jt)=klis(jt)+1
         vx=(x0-uint(iglo,4,1))/rglo
         vy=(y0-uint(iglo,4,2))/rglo
         vz=(z0-uint(iglo,4,3))/rglo
      posn=1+float(iglo-1)/nspl
      pari(ks,1)=posn
      pari(ks,2)=rglo
      px=uint(iglo,2,1)
      py=uint(iglo,2,2)
      pz=uint(iglo,2,3)
      dot=vx*px+vy*py+vz*pz
      pari(ks,3)=acos(dot)
      dx=py*vz-pz*vy
      dy=pz*vx-px*vz
      dz=px*vy-py*vx
      dot=dx*uint(iglo,3,1)+dy*uint(iglo,3,2)+dz*uint(iglo,3,3)
      if(dot.lt.0) pari(ks,3)=-pari(ks,3)
*-------------------------------------------------add vector to axis O/P
      if(.not.traj) then
      kas=kas+1
      cors(kas,1)=x0
      cors(kas,2)=y0
      cors(kas,3)=z0
      snam(kas)=mnam(in)
      sunit(kas)='ION'
      nunis(kas)=nunit(in)
         mats(khs+1)=80000
         mats(khs+2)=iglo
         mats(khs+3)=kas
         khs=khs+3
      endif
c--------------------------------------------------------calculate angle
      if(.not.traj.or.(itst.eq.itnd)) write(6,23) ks,mnam(in),
     & nunit(in),pari(ks,1),pari(ks,2),pari(ks,3)*crd
23    format(i3,') ',a4,i3,2f8.2,f8.1)
21    enddo
      kisa=ks
      return
      end
