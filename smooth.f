      subroutine smooth
      include 'curves_data.inc'
      parameter (iw=4)          !... smoothing width
      character*1 na,nt
      character*4 nback
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfr,frames,lpa
      dimension w(-5:5),upa(6)
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      common/geo/ref(n3,4,5,3),rel(n6,4,3),upm(n3,4,3),plig(n6),
     1 ilig(n6),klig,nback(4)
      common/hel/upl(n3,0:8,6),uvw(n3,4,3),npl(n3),lpa(n3,4)
      common/mat/ni(n3,4),nu(n3,4),idr(4),nst,nlev,na(n3,4),nt(n3,4)
 
crl   UVW(i,m,j) contains the axis reference system for each level i
crl   m=1 pseudodyad (pd), m=2 perp. to pd, m=3 helical axis, m=4 ref. point
 
      if(test) write(6,5)
5     format(/2x,'Axis calculation ...'/)
c-------------------------------------------------------------------zaxe
      if(zaxe) then
      do i=1,nlev
         do m=1,4
         do j=1,3
         uvw(i,m,j)=0.d0
         enddo
         enddo
      uvw(i,3,3)=1.d0
      uvw(i,4,3)=upm(i,4,3)
      enddo
      if(test) write(6,10) (uvw(1,3,j),j=1,3),(uvw(1,4,j),j=1,3)
 10   format(2x,'Zaxe  u1: ',3f7.3,' p1: ',3f7.3)
c-------------------------------------------------------------------line
      else if(line) then
         if(circ) then
         write(6,11)
11       format(/2x,'Option Line=.t. has no sense for circular DNA')
         stop
         endif
      uax=0.d0
      uay=0.d0
      uaz=0.d0
      pax=0.d0
      pay=0.d0
      paz=0.d0
         ic=0
         do i=1,nlev
         do k=1,npl(i)
         ic=ic+1
         uax=uax+upl(i,k,1)
         uay=uay+upl(i,k,2)
         uaz=uaz+upl(i,k,3)
         pax=pax+upl(i,k,4)
         pay=pay+upl(i,k,5)
         paz=paz+upl(i,k,6)
         enddo
         enddo
      r=sqrt(uax*uax+uay*uay+uaz*uaz)
      uax=uax/r
      uay=uay/r
      uaz=uaz/r
      pax=pax/ic
      pay=pay/ic
      paz=paz/ic
         do i=1,nlev
         uvw(i,3,1)=uax
         uvw(i,3,2)=uay
         uvw(i,3,3)=uaz
         dx=upm(i,4,1)-pax
         dy=upm(i,4,2)-pay
         dz=upm(i,4,3)-paz
         dot=dx*uax+dy*uay+dz*uaz
         uvw(i,4,1)=pax+dot*uax
         uvw(i,4,2)=pay+dot*uay
         uvw(i,4,3)=paz+dot*uaz
         enddo
      if(test) write(6,15) (uvw(1,3,j),j=1,3),(uvw(1,4,j),j=1,3)
 15   format(2x,'Line  u1: ',3f7.3,' p1: ',3f7.3)
      else
c-------------------------------------------------------curvilinear axis
            w(0)=1.d0
            do i=1,iw
c           w(i)=cos(i*pi/(2.d0*(iw+1.d0))) !... alt. trig. smo. fun.
            w(i)=1.-float(i)**2/(iw+1.d0)**2
            w(-i)=w(i)
            enddo
c----------------------------------------------------------------average
         do i=1,nlev
         uax=0.d0
         uay=0.d0
         uaz=0.d0
         pax=0.d0
         pay=0.d0
         paz=0.d0
         do k=1,npl(i)
         uax=uax+upl(i,k,1)
         uay=uay+upl(i,k,2)
         uaz=uaz+upl(i,k,3)
         pax=pax+upl(i,k,4)
         pay=pay+upl(i,k,5)
         paz=paz+upl(i,k,6)
         enddo
         r=sqrt(uax*uax+uay*uay+uaz*uaz)
         uax=uax/r
         uay=uay/r
         uaz=uaz/r
         upl(i,0,1)=uax
         upl(i,0,2)=uay
         upl(i,0,3)=uaz
         pax=pax/npl(i)
         pay=pay/npl(i)
         paz=paz/npl(i)
         dx=upm(i,4,1)-pax
         dy=upm(i,4,2)-pay
         dz=upm(i,4,3)-paz
         dot=dx*uax+dy*uay+dz*uaz
         upl(i,0,4)=pax+dot*uax
         upl(i,0,5)=pay+dot*uay
         upl(i,0,6)=paz+dot*uaz
         enddo
c-----------------------------------------------------------------smooth
         do i=1,nlev
         x0=upm(i,4,1)
         y0=upm(i,4,2)
         z0=upm(i,4,3)
         dx0=upm(i,1,1)
         dy0=upm(i,1,2)
         dz0=upm(i,1,3)
         wt=0.
         ic=0.
         uax=0.d0
         uay=0.d0
         uaz=0.d0
         pax=0.d0
         pay=0.d0
         paz=0.d0
         do id=-iw,iw
         is=i+id
            if(is.lt.1) then
            if(.not.circ) goto 20
            is=is+nlev
            endif
            if(is.gt.nlev) then
            if(.not.circ) goto 20
            is=is-nlev
            endif
         ic=ic+1
         ux=upl(is,0,1)
         uy=upl(is,0,2)
         uz=upl(is,0,3)
         px=upl(is,0,4)
         py=upl(is,0,5)
         pz=upl(is,0,6)
         uax=uax+ux*w(id)
         uay=uay+uy*w(id)
         uaz=uaz+uz*w(id)
         dot=(x0-px)*ux+(y0-py)*uy+(z0-pz)*uz
         pax=pax+(px+dot*ux)*w(id)
         pay=pay+(py+dot*uy)*w(id)
         paz=paz+(pz+dot*uz)*w(id)
         wt=wt+w(id)
20       continue
         enddo
         r=sqrt(uax*uax+uay*uay+uaz*uaz)
         uax=uax/r
         uay=uay/r
         uaz=uaz/r
         uvw(i,3,1)=uax
         uvw(i,3,2)=uay
         uvw(i,3,3)=uaz
         pax=pax/wt
         pay=pay/wt
         paz=paz/wt
         dot=(x0-pax)*uax+(y0-pay)*uay+(z0-paz)*uaz
         uvw(i,4,1)=pax+dot*uax
         uvw(i,4,2)=pay+dot*uay
         uvw(i,4,3)=paz+dot*uaz
         if(test) write(6,25) i,(uvw(i,3,j),j=1,3),(uvw(i,4,j),j=1,3)
 25      format(2x,'Smth ',i3,' U: ',3f8.2,' P: ',3f8.2)
         enddo
      endif
      return
      end
