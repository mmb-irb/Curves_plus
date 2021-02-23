      subroutine ligpar
      include 'curves_data.inc'
      character*1 na,nt,base,type,first
      character*4 mnam,munit,snam,sunit,name,ban,bln,
     1 liga,nback,inam
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
      common/lgd/ilnam(20),paxe(20),blef(20,40,3),parl(20,9),
     1 ilref(20),iluni(20,2),ilfrm(20,-3:20),bln(20,40),nlig,liga(20)
      common/mac/corm(n1,3),mnam(n1),munit(n1),nunit(n1),
     1 iunit(n1),ncen(0:n2),kam,kcen
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
      common/sto/bpdat(n3,n5),traj
      if(.not.traj) write(6,4)
4     format(/2x,'(L) Ligands ',5x,'Lshif',3x,'Lslid',3x,'Lrise',3x,
     1 'Ltilt',3x,'Lroll',3x,'Ltwis',3x,'L-Axe',3x,'L-Dis',3x,'L-Ang'/)
      kaso=kas
      do kl=1,klig
      x0=rel(kl,4,1)
      y0=rel(kl,4,2)
      z0=rel(kl,4,3)
      if(plig(kl).gt.0) then
      ilev=int(plig(kl))
      frac=plig(kl)-ilev
      imin=1+(ilev-1)*nspl+nint(nspl*frac)
      else
c--------------------------------------find perpendicular vector to axis
         amin=1.
         do i=1,kaso
         dx=x0-uint(i,4,1)
         dy=y0-uint(i,4,2)
         dz=z0-uint(i,4,3)
         r=sqrt(dx*dx+dy*dy+dz*dz)
         dot=(dx*uint(i,3,1)+dy*uint(i,3,2)+dz*uint(i,3,3))/r
         if(abs(dot).lt.amin) then
         rmin=r
         amin=dot
         imin=i
         endif
         enddo
         if(abs(amin).gt.acrit) goto 50
      posn=1+float(imin-1)/nspl
      ilev=int(posn)
      frac=posn-ilev
      endif
      do i=1,4
      do j=1,3
      r1(i,j)=rel(kl,i,j)
      r2(i,j)=uint(imin,i,j)
      enddo
      enddo
*-------------------------------------------------add vector to axis O/P
      if(test) then
      do j=1,3
      cors(kas+1,j)=uint(imin,4,j)
      cors(kas+2,j)=rel(kl,4,j)
      enddo
      mats(khs+1)=80000
      mats(khs+2)=kas+1
      mats(khs+3)=kas+2
      do j=1,2
      kas=kas+1
      snam(kas)='H'
      sunit(kas)='LIG'
      nunis(kas)=kl
      enddo
      khs=khs+3
      endif
c---------------------------------------------------calculate parameters
      call screw(r1,r2,0)
      ux=uint(imin,3,1)
      uy=uint(imin,3,2)
      uz=uint(imin,3,3)
      px=uint(imin,2,1)
      py=uint(imin,2,2)
      pz=uint(imin,2,3)
      dx=x0-uint(imin,4,1)
      dy=y0-uint(imin,4,2)
      dz=z0-uint(imin,4,3)
      dot=dx*ux+dy*uy+dz*uz
      dx=dx-dot*ux
      dy=dy-dot*uy
      dz=dz-dot*uz
      r=sqrt(dx*dx+dy*dy+dz*dz)
      vx=dx/r
      vy=dy/r
      vz=dz/r
      dot=vx*px+vy*py+vz*pz
      angle=acos(dot)
      dx=py*vz-pz*vy
      dy=pz*vx-px*vz
      dz=px*vy-py*vx
      dot=dx*ux+dy*uy+dz*uz
      if(dot.lt.0) angle=-angle
*-------------------------------------------------------store parameters
        do j=1,3
        parl(kl,j)=var(j)
        enddo
        do j=4,6
        parl(kl,j)=var(j)*cdr
        enddo
        parl(kl,7)=ilev+frac
        parl(kl,8)=dist
        parl(kl,9)=angle
      is=ilig(kl)
      if(.not.traj) write(6,22) kl,munit(is),nunit(is),
     1 (var(j),j=1,6),ilev+frac,dist,angle*crd
22    format(2x,i3,') ',a3,i4,2x,3f8.2,3f8.1,2f8.2,f8.1)
50    enddo
      return
      end
