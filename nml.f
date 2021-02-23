      subroutine nml
      include 'curves_data.inc'
      parameter (n_real=2,n_int=5,n_log=9,n_cha=8,n_tot=24)
      character vp(n_cha)*32,input(n_tot)*10,date*25,day*9,lj*1
      character*128 file,ftop,lis,lib,lig,ibld,sol,back,
     1 lini,vc(n_cha)
      logical*2 circ,line,zaxe,fit,test,ions,refo,axfrm,frames,
     1 vo(n_log),iflag(n_tot),first,last,start
      integer*4 vi(n_int),nmls(n_tot)
      dimension vr(n_real)
      common/cha/file,ftop,lis,lib,lig,ibld,sol,back
      common/dat/wback,wbase,isym,itst,itnd,itdel,naxlim,
     1 circ,line,zaxe,fit,test,ions,refo,axfrm,frames
      equivalence (vc(1),file),(vr(1),wback),(vi(1),isym),(vo(1),circ)
      data input/'wback','wbase','isym','itst','itnd','itdel',
     1 'naxlim','circ','line','zaxe','fit','test','ions','refo',
     1 'axfrm','frames','file','ftop','lis','lib','lig','ibld',
     1 'sol','back'/
      ninr=n_int+n_real
      nlog=n_log+ninr
         do i=1,n_cha-1
         vc(i)=' '
         enddo
         do i=1,n_tot
         iflag(i)=.false.
         nmls(i)=index(input(i),' ')-1
         enddo
      first=.true.
      last=.false.
 10   read(5,5) lini
 5    format(a)
      im=index(lini,'&')
      if(im.gt.0) then
      if(.not.first) last=.true.
      if(first.and.index(lini(im+1:),'&').ne.0) last=.true.
      endif
      do k=1,128
      if(lini(k:k).eq.'=') then
      kl=k
      start=.true.
      do j=k-1,1,-1
      lj=lini(j:j)
      if(start.and.lj.ne.' ') then
      start=.false.
      jh=j
      else if(.not.start.and.(lj.eq.' '.or.lj.eq.',')) then
      jl=j+1
      goto 15
      endif
      enddo
      goto 50
 15      do i=1,n_tot
         if(lini(jl:jh).eq.input(i)) then
         iflag(i)=.true.
 17         kl=kl+1
            if(lini(kl:kl).eq.' ') goto 17
            do j=kl,128
            lj=lini(j:j)
            if(lj.eq.' '.or.lj.eq.','.or.lj.eq.'&') then
            kh=j-1
            goto 19
            endif
            enddo
 19      if(i.le.n_real) then
         read(lini(kl:kh),*,err=50) vr(i)
         goto 25
         else if(i.le.ninr) then
         read(lini(kl:kh),*,err=50) vi(i-n_real)
         goto 25
         else if(i.le.nlog) then
         read(lini(kl:kh),*,err=50) vo(i-ninr)
         goto 25
         else
         if(lini(kl:kl).eq.'''') kl=kl+1
         if(lini(kh:kh).eq.'''') kh=kh-1
         if(kh.ge.kl) read(lini(kl:kh),5,err=50) vc(i-nlog)(1:kh-kl+1)
         goto 25
         endif
         endif
         enddo
         goto 50
      endif
 25   enddo
      first=.false.
      if(.not.last) goto 10
c-----------------------------------------------------------------output
         if(vi(2).gt.1.and.vi(3).eq.1) vi(3)=vi(2)
         if(vc(1).ne.' ') then
         kfi=index(vc(1),' ')-1
         do i=kfi,1,-1
         if(vc(1)(i:i).eq.'.') goto 4
         enddo
         vc(1)=vc(1)(:kfi)//'.pdb'
4        kfi=index(vc(3),' ')-1
         open(unit=6,file=vc(3)(:kfi)//'.lis',status='new')
         endif
      call getdate(date)
      day=date(9:10)//date(4:8)//date(23:)
      write(6,200) day
200   format(
     1/5x,'**************************************',15x,'**************',
     1/5x,'**** CURVES+ Version 2.6  04/2014 ****',15x,'*  ',a9,   ' *',
     1/5x,'**************************************',15x,'**************',
     1 //)
      do i=1,n_tot
      nm=nmls(i)
      if(iflag(i)) then
      do j=1,nm
      ic=ichar(input(i)(j:j))
      input(i)(j:j)=char(ic-32)
      enddo
      endif
      enddo
         do i=1,n_cha
         is=1
         ie=0
         do j=128,1,-1
         if(i.lt.n_cha.and.is.eq.1.and.vc(i)(j:j).eq.'/') is=j+1
         if(ie.eq.0.and.vc(i)(j:j).eq.' ') ie=j-1
         enddo
         vp(i)=vc(i)(is:ie)
         enddo
      write(6,8) (input(nlog+j),vp(j),j=1,n_cha)
8     format(2x,a5,': ',a32,:,2x,a5,': ',a32)
      write(6,*)
      write(6,12) (input(j),vr(j),j=1,n_real)
12    format(4(:,2x,a6,':',f6.2))
      write(6,*)
      write(6,14) (input(n_real+j),vi(j),j=1,n_int)
14    format(5(:,2x,a6,':',i6))
      write(6,*)
      write(6,16) (input(ninr+j),vo(j),j=1,n_log)
16    format(5(:,2x,a6,':',l6))
 
      write(6,*)
      return
50    write(6,55) lini(jl:jh)
55    format(/2x,'---- Error in namelist input for ',a,' ----'/)
      stop
      end
