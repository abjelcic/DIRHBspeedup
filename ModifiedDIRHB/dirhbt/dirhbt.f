c======================================================================c
      subroutine base(lpr)
c======================================================================c
c
c     determines the basis for Dirac equation in cartesian coordinates
c
c     NB        number of parity-blocks
c     IA(ib,1): begin of the large components of block b is IA(b,1)+1
c     IA(ib,2): begin of the small components of block b is IA(b,2)+1
c     ID(ib,1): dimension large components of block b
c     ID(ib,2): dimension small components of block b
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      dimension ipar(2),isim(2)
c
      common /basnnn/ n0f,n0b
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /bconf / no,nxyzb(nobx,3)
      common /gfviv / iv(0:igfv)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /icalc/ ivpair
      common /numsep/ npm
c
      data ipar/1,2/,isim/1,1/
c
      if (lpr) then
      write(l6,*) '****** BEGIN BASE **********************************'
      endif
c
c======================================================================c
c     Oscillator-Base for Fermions:
c======================================================================c
c----------------------------------------------------------------------c
c     Construction of the different Parity blocks.
c----------------------------------------------------------------------c
      nb  = 2
      nxm = 0
      nym = 0
      nzm = 0
      nfx0 = 0
      ngx0 = 0
c
      il = 0
      ik = 0
      mv = 0
c
c     loop over parity-simplex blocks
      do ib = 1,nb
	     ipf = ipar(ib)
         isi = isim(ib)
	     mu(ib) = 2
         write(tb(ib),'(2a1)') tp(isi),tp(ipf)
c
c        loop over large and small components
         do ifg = 1,2
	        ia(ib,ifg) = il
            if (ifg.eq.1) then
               ip = ipf
            else
               ip = 3 - ipf
            endif
c
c           loop over major quantum number nn
            do nn = ip-1,n0f+ifg-1,2
c           loop over quantum number nx
               do ix = 0,nn
c           loop over quantum number ny
                  do iy = 0,nn
c           loop over quantum number nz
                     do iz = 0,nn
                        if (ix+iy+iz.eq.nn) then
                           il = il + 1
	                       if (il.gt.ntx) stop
     &                        'in BASE:ntx too small'
                           nxyz(il,1) = ix
                           nxyz(il,2) = iy
                           nxyz(il,3) = iz
                           ns(il)     = isi
	                       np(il)     = ifg
c                           if (nn.lt.10) then
                              write(tt(il),100) ix,iy,iz!,isi,tp(ip)
c                           else
c                              write(tt(il),101) ix,iy,iz,isi
c                           endif
                           nxm = max0(nxm,ix)
                           nym = max0(nym,iy)
                           nzm = max0(nzm,iz)
                        endif
   10                enddo   ! iz
                  enddo   ! iy
               enddo   ! ix
            enddo   ! nn
            id(ib,ifg) = il - ia(ib,ifg)
	     enddo   ! ifg
	     nf = id(ib,1)
	     do n2=1,nf
	        do n1=n2,nf
	           mv = mv + 1
	        enddo
	     enddo

	     nfx0 = max(nfx0,id(ib,1))
	     ngx0 = max(ngx0,id(ib,2))
	     ik = ik + max0(id(ib,1),id(ib,2))
c
      enddo   ! ib
      nt = il
      nk = ik

      nsep = 0
      do ix=0,npm,2
         do iy=0,npm-ix,2
            do iz=0,npm-ix-iy,2
               nsep = nsep + 1
            enddo
         enddo
      enddo

      if (nxm.gt.nmax)   stop 'in BASE: nxx too small '
      if (nym.gt.nmax)   stop 'in BASE: nyx too small '
      if (nzm.gt.nmax)   stop 'in BASE: nzx too small '
      if (nfx0.gt.nfx)   stop 'in BASE: nfx too small '
      if (ngx0.gt.ngx)   stop 'in BASE: ngx too small '
      if (nk.gt.nkx)     stop 'in BASE: nkx too small '
      if (nsep.gt.nsepx) stop 'in BASE: nsepx too small '
      if (mv.gt.mvfx)    stop 'in BASE: mvfx too small '
c
c
c
c
      if (lpr) then
         write(l6,*) ' '
         write(l6,*)   ' Maximal values:             needed     given'
         write(l6,103) ' Number of blocks: nb  = ',nb,nbx
         write(l6,103) ' Number of levels  nt  = ',nt,ntx
         write(l6,103) ' Number of levels  nk  = ',nk,nkx
         write(l6,103) ' Maximal nx:       nxm = ',nxm,nmax
         write(l6,103) ' Maximal ny:       nym = ',nym,nmax
         write(l6,103) ' Maximal nz:       nzm = ',nzm,nmax
         write(l6,103) ' Maximal nf            = ',nfx0,nfx
         write(l6,103) ' Maximal ng            = ',ngx0,ngx
         write(l6,103) ' Number of two-body pairs  = ',mv,mvx
      endif
c
c----------------------------------------------------------------------c
c     Printout
c----------------------------------------------------------------------c
      if (lpr) then
         do ib = 1,nb
	        nf  = id(ib,1)
	        ng  = id(ib,2)
	        nh  = nf + ng
	        mf  = ia(ib,1)
	        mg  = ia(ib,2)
            write(l6,'(/,a,6i4)') tb(ib),nf,ng,mf,mg
            do i = mf+1,mf+nh
               write(l6,102) i,'   nx = ',nxyz(i,1),
     &                         '   ny = ',nxyz(i,2),
     &                         '   nz = ',nxyz(i,3),
     &                         '  sim = ',ns(i),' ',tt(i)
	           if (i.eq.mf+nf) write(l6,'(3x,61(1h-))')
            enddo  ! i
         enddo   ! ib
      endif
c
c
c
c======================================================================c
c     Oscillator Base for Bosons:
c======================================================================c
      if (mod(n0b,2).ne.0) stop ' in BASE: n0b must be even'
c
c---- boson basis even in x,y,z (densities)
      call boson(n0b,no,nxyzb,lpr)

      nobm  = no
      nmax0 = max(nxm,nym,nzm,n0b)
c
c
c  100 format(1h[,3i1,1h],i2,a1)
  100 format(1h[,3i2,1h])
  101 format(4i2)
  102 format(i4,a,i2,a,i2,a,i2,a,i2,a,a)
  103 format(a,2i10)
c
      if (lpr) then
      write(l6,*) '****** END BASE ************************************'
      endif

      return
c-end-BASE
      end
c======================================================================c

      subroutine boson(nnmax,nbos,nxyzb,lpr)

c======================================================================c
c
c     determines the boson basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      dimension nxyzb(nobx,3)
c
      common /basnnn/ n0f,n0b
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      il   = 0
      nxbm = 0
      nybm = 0
      nzbm = 0
c
      do ix = 0,n0b,2
      do iy = 0,n0b,2
      do iz = 0,n0b,2
c
         if (ix+iy+iz.le.n0b) then
            il = il + 1
            if (il.gt.nobx) stop ' in BASE: nobx too small'
            nxyzb(il,1) = ix
            nxyzb(il,2) = iy
            nxyzb(il,3) = iz

            nxbm = max0(nxbm,ix)
            nybm = max0(nybm,iy)
            nzbm = max0(nzbm,iz)
            if (nxbm.gt.nmax) stop 'in BASE: nxbx too small '
            if (nybm.gt.nmax) stop 'in BASE: nybx too small '
            if (nzbm.gt.nmax) stop 'in BASE: nzbx too small '
         endif
      enddo   ! iz
      enddo   ! iy
      enddo   ! ix
      nbos  = il
c
      if (lpr) then
         write(l6,'(//,a,/)') ' Boson base: '
         do i = 1,nbos
            write(l6,102) i,nxyzb(i,1),nxyzb(i,2),nxyzb(i,3)
         enddo
         write(l6,*) ' '
         write(l6,*)   'Maximal values:             needed     given'
         write(l6,103) ' Number of levels no   = ',nbos,nobx
      endif
c
c
  102 format(i5,'[',i2,i2,i2,']')
  103 format(a,2i10)
c
      return
c-end-BOSON
      end

c======================================================================c

      subroutine broyden(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character tb*5                                           ! blokap
c
      ! broyden iteration sizes
      parameter (nn = 2*MVFX+2*MVGX+2*MVFX+2+2)
      parameter (mm = 7)
c
      common /baspar/ hom,hb0,b0
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /physco/ hbc,alphi,r0
      common /deldel/ de(nhhx,nb2x)
      common /hamham/ hh(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /cstr2/ q0c,q2c,c0,c2,alaq0,alaq2
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /broyde/ bbeta(mm,mm),df(nn,mm),dv(nn,mm),
     &                bwork(mm),curv(nn),bw0,ibwork(mm)
      common /broyde1/ vin(nn)
      common /broyde2/ ibroyd
      dimension vou(nn)
      data bmix /0.7d0/
c
      if (lpr) then
      write(l6,*) '****** BEGIN BROYDEN *************************'
      endif
      xmi = xmix
      si  = zero
      do i=1,nn
         vou(i)=zero
      enddo
c
      if (ii.eq.0) then
         do i=1,nn
            vin(i)=zero
         enddo
         ipos=0
         do it = 1,itx
            do ib=1,nb
               nf=id(ib,1)
               ng=id(ib,2)
               nh=nf+ng
               m = ib+(it-1)*nbx
               do i2 = 1,nf
                  do i1 = i2,nf
                     ipos=ipos+1
                     vin(ipos) = hh(i1+(i2-1)*nh,m)
                  enddo
               enddo
               do i2 = 1,ng
                  do i1 = i2,ng
                     ipos=ipos+1
                     vin(ipos) = hh(nf+i1+(nf+i2-1)*nh,m)
                  enddo
               enddo
               do i2 = 1,nf
                  do i1 = i2,nf
                     ipos=ipos+1
                     vin(ipos) = de(i1+(i2-1)*nh,m)
                  enddo
               enddo
            enddo  !ib
            ipos=ipos+1
            vin(ipos)=ala(it)
         enddo  !it
         ipos=ipos+1
         vin(ipos)=alaq0
         ipos=ipos+1
         vin(ipos)=alaq2
         return
      endif
      ipos=0
      do it = 1,itx
         do ib=1,nb
            nf=id(ib,1)
            ng=id(ib,2)
            nh=nf+ng
            m = ib+(it-1)*nbx
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  vou(ipos) = hh(i1+(i2-1)*nh,m) -vin(ipos)
               enddo
            enddo
            do i2 = 1,ng
               do i1 = i2,ng
                  ipos=ipos+1
                  vou(ipos) = hh(nf+i1+(nf+i2-1)*nh,m) -vin(ipos)
               enddo
            enddo
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  vou(ipos) = de(i1+(i2-1)*nh,m) -vin(ipos)
               enddo
            enddo
         enddo
         ipos=ipos+1
         vou(ipos)=ala(it)-vin(ipos)
      enddo
      ipos=ipos+1
      vou(ipos)=alaq0-vin(ipos)
      ipos=ipos+1
      vou(ipos)=alaq2-vin(ipos)
      ! calculate si
      do i=1,nn
         si = max( si, abs(vou(i)) )
      enddo
      ! broyden's mixing procedure starts here...
      if ( mm .eq. 0 .or. ii.eq. 1 .or. ibroyd.eq.0) then ! linear mixing
         do i = 1, nn
            vin(i) = vin(i) + xmi*vou(i)
         enddo
         ilast = 0
      else                                                ! broyden mixing

         iuse = min( ii-1 - 1, mm )
         ipos = ii - 2 - ( (ii-3)/mm )*mm
         inex = ii - 1 - ( (ii-2)/mm )*mm

         if( ii .eq. 2 ) then
            do j = 1, mm
               do i = 1, nn
                  df(i,j) = zero
                  dv(i,j) = zero
               enddo
            enddo
            bw0 = 0.01d0
         else
            do i = 1, nn
               df(i,ipos) = vou(i) - df(i,ipos)
               dv(i,ipos) = vin(i) - dv(i,ipos)
            enddo

            dnorm = sqrt( dnrm2(nn,df(1,ipos),1)**2.0d0 )
            call dscal( nn, 1.0d0 / dnorm, df(1,ipos), 1 )
            call dscal( nn, 1.0d0 / dnorm, dv(1,ipos), 1 )
         endif

         do i = 1, iuse
            do j = i+1, iuse
               bbeta(i,j) = ddot( nn, df(1,j), 1, df(1,i), 1 )
            enddo
            bbeta(i,i) = 1.0d0 + bw0*bw0
         enddo

         call dsytrf( 'U', iuse, bbeta, mm, ibwork, bwork, mm, info )
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRF V+S '

         call dsytri( 'U', iuse, bbeta, mm, ibwork, bwork, info )
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRI V+S '

         do i = 1, iuse
            do j = i+1, iuse
               bbeta(j,i) = bbeta(i,j)
            enddo
            bwork(i) = ddot( nn, df(1,i), 1, vou, 1 )
         enddo

         do i = 1, nn
            curv(i) = bmix * vou(i)
         enddo
         do i = 1, iuse
            gamma = 0.0d0
            do j = 1, iuse
               gamma = gamma + bbeta(j,i) * bwork(j)
            enddo
            do k = 1, nn
               curv(k) = curv(k) - gamma * ( dv(k,i) + bmix * df(k,i) )
            enddo
         enddo

         call dcopy( nn, vou, 1, df(1,inex), 1 )
         call dcopy( nn, vin, 1, dv(1,inex), 1 )

         curvature = ddot( nn, vou, 1, curv, 1 )
         if( curvature .gt. 0.0d0 ) then
            ilast = 1
            do i = 1, nn
               vin(i) = vin(i) + curv(i)
            enddo
c            write(*,100) 'broyden mixing: mm =', iuse, 'c=',curvature
         else
            ilast = 0
            do i = 1, nn
               vin(i) = vin(i) + xmi*vou(i)
            enddo
c            write(*,100) 'linear  mixing: mm =', iuse, 'c=',curvature
         endif
      endif
 100  format(10x,a,i2,2x,a,f16.8)
      ! broyden's mixing procedure ends here

      ! set the new matrix elements
      ipos=0
      do it = 1,itx
         do ib=1,nb
            nf=id(ib,1)
            ng=id(ib,2)
            nh=nf+ng
            m = ib+(it-1)*nbx
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  hh(i1+(i2-1)*nh,m) = vin(ipos)
                  hh(i2+(i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
             do i2 = 1,ng
               do i1 = i2,ng
                  ipos=ipos+1
                  hh(nf+i1+(nf+i2-1)*nh,m) = vin(ipos)
                  hh(nf+i2+(nf+i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  de(i1+(i2-1)*nh,m) = vin(ipos)
                  de(i2+(i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
         enddo  !ib
         ipos=ipos+1
         ala(it) =vin(ipos)
      enddo  !it
      ipos=ipos+1
      alaq0=vin(ipos)
      ipos=ipos+1
      alaq2=vin(ipos)
c
c
      if (lpr) then
      write(l6,*) '****** END BROYDEN *************************'
      endif
c
      return
c-end-BROYDEN
      end

c=====================================================================c

      subroutine canon(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpr

      character*8 tbb(nhfbx)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tk*8,tx*8                                       ! blolev
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
c
      dimension aa(nhhx),dd(nhhx),v2(nhx),z(nhx),eb(nhx),h(nhx),d(nhx)
      dimension aaka(nhhx),ak(nhx)

      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4),akacan(nkx,4)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /gfvsq / sq(0:igfv)
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /hamham/ hh(nhhx,nb2x)
      common /deldel/ de(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
c
      data ash/100.0/
      if(lpr)
     &write(l6,*) ' ****** BEGIN CANON ********************************'
c
c
c---- BCS structure in the canonical basis:
c
      do it = 1,2
         if (lpr) write(l6,100) tit(it)
  100    format(' single-particle energies and gaps ',1x,
     &          'in the canonical basis: ',a,/,1x,66(1h-))
         if (lpr) write(l6,1000) ' ',' N pi ','  [ nx ny nz]',
     &                           'smax','eecan','vvcan','decan'
1000     format(a2,a6,a13,2a11,2x,2a11)
c
c---- loop of the j-blocks

         klp = 0
         kla = 0
         do ib = 1,nb
	        mul = mu(ib)
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            nhfb= nh + nh
            i0f = ia(ib,1)
            i0g = ia(ib,2)
	        mf   = ib + (it-1)*nbx
            kf  = kd(ib,it)
            kg  = kd(ib,it+2)
	        k0f = ka(ib,it)
	        k0g = ka(ib,it+2)

c------- transformation to the canonical basis
c------- calculation of the generalized density V*VT
            do n2 = 1,nh
               do n1 = 1,nh
                  s = zero
                  saka=zero
                  do k = k0f+1,k0f+kf
                     s=s+fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
                     saka = saka
     &                    + 0.5*(fguv(n1,k,it)*fguv(nh+n2,k,it)+
     &                           fguv(n2,k,it)*fguv(nh+n1,k,it))
                  enddo
                  s1 = zero
                  do k = k0g+1,k0g+kg
                     s1 = s1 + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
                  enddo
                  aa(n1+(n2-1)*nh) = s + ash*s1
                  if(n2 .le. nf .and. n1 .le. nf) then
                     aaka(n1+(n2-1)*nh) = saka
                  else
                     aaka(n1+(n2-1)*nh) = 0.d0
                  endif
               enddo
            enddo   ! n2
            call sdiag(nh,nh,aa,v2,dd,z,1)

            eps=1.0e-6
            call degen(nh,nh,v2,dd,hh(1,mf),eb,eps,aa,z)
c
c        major component of the wave function should be > 0
            do k = 1,nh
	           cmax = zero
	           do i = 1,nh
	              if (abs(dd(i+(k-1)*nh)).gt.abs(cmax))
     &               cmax = dd(i+(k-1)*nh)
               enddo   ! i
	           if (cmax.lt.zero) then
	              do i = 1,nh
		             dd(i+(k-1)*nh) = - dd(i+(k-1)*nh)
	              enddo   ! i
	           endif
	        enddo   ! k
c------- diagonal matrix elements of HH and DE in the canonical basis
            do k = 1,nh
               hk = zero
               dk = zero
               akk=zero
               do n2 = 1,nh
                  h2 = zero
                  d2 = zero
                  akk2=zero
                  do n1 = 1,nh
                     h2 = h2 + dd(n1+(k-1)*nh)*hh(n1+(n2-1)*nh,mf)
                     d2 = d2 + dd(n1+(k-1)*nh)*de(n1+(n2-1)*nh,mf)
                     akk2=akk2 + dd(n1+(k-1)*nh)*aaka(n1+(n2-1)*nh)
                  enddo
                  hk = hk + h2*dd(n2+(k-1)*nh)
                  dk = dk + d2*dd(n2+(k-1)*nh)
                  akk = akk + akk2*dd(n2+(k-1)*nh)
               enddo
               h(k) = hk
               d(k) = dk
               ak(k)=akk
            enddo   ! k

c------- reordering according to the energy h(k)
            call ordx2(nh,h,d,v2,ak,dd)

            do k = 1,nh
               if (v2(k).lt.zero .or. v2(k).gt.2.d0) v2(k) = zero
               if (v2(k).gt.one) v2(k) = one
            enddo ! k

            kacan(ib,it)   = klp
            do k=1,nf
               klp=klp+1
               eecan(klp,it)=h(ng+k)
               decan(klp,it)=d(ng+k)
               vvcan(klp,it)=v2(ng+k)
               akacan(klp,it)=ak(ng+k)
               mulcan(klp,it)=mul
               do i=1,nh
                  fgcan(i,klp,it)=dd(i+(ng+k-1)*nh)
               enddo
            enddo
            kdcan(ib,it) = klp - kacan(ib,it)

            kacan(ib,it+2) = kla
            do k=1,ng
               kla=kla+1
               eecan(kla,it+2)=h(k)
               decan(kla,it+2)=d(k)
               vvcan(kla,it+2)=v2(k)
               akacan(kla,it+2)=ak(k)
               mulcan(kla,it+2)=mul
               do i=1,nh
                  fgcan(i,kla,it+2)=dd(i+(k-1)*nh)
               enddo
            enddo
            kdcan(ib,it+2) = kla - kacan(ib,it+2)
c------- printout for particles
            if (ib.eq.1) e0 = h(ng+1)
            k1 = kacan(ib,it)+1
            k2 = kacan(ib,it)+kdcan(ib,it)
            do k = k1,k2
               e1 = eecan(k,it)
               v1 = vvcan(k,it)
               d1 = decan(k,it)
               ak1=akacan(k,it)
c           search for the main oscillator component
               smax = zero
               do i = 1,nh
                  s = abs(fgcan(i,k,it))
                  if (s.gt.smax) then
                     smax = s
                     imax = i
                  endif
               enddo
               dx = fgcan(imax,k,it)
c
               if (lpr) write(l6,101) k,tp(ib),
     &                    tt(i0f+imax),smax,
     &                    e1,v1,d1
  101 format(i4,'  ',a1,'  ','  ',a8,'  ',6x,f5.2,3f12.4)
           enddo  ! k
c---- end loop over ip-blocks
   10    enddo !ib
         nkcan(it)   = klp
         nkcan(it+2) = kla

      enddo   ! it
c
  102 if (lpr)
     &write(l6,*) ' ****** END CANON **********************************'
c
      return
c-end-CANON
      end

c======================================================================c
c
      subroutine degen(na,n,ea,dd,bb,eb,eps,zz,z)
c
c======================================================================c
c
c     EA    is a set of partially degenerate eigenvalues of some matrix
c     DD    are the corresponding eigenvectors DD.
c     BB    is a matrix, which is diagonalized in the subspaces of
c           degenerate eigenvales EA.
c     EB    contains the eigenvalues of BB in these subspaces.
c     EPS   determines, to which accuracy the eigenvalues EA are
c           considered to be degenerate
c     ZZ,Z  are auxiliary arrays
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      dimension bb(na,n),dd(na,n),ea(n),eb(n)
      dimension zz(na,n),z(n)
c
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
c---- check for degeneracies
      k1 = 0
      do i1 = 1,n
         k1 = k1 + 1
         if (k1.ge.n) goto 20
         do k2 = k1+1,n
            if (ea(k2)-ea(k1).gt.eps) goto 10
         enddo
   10    me = k2 - 1
         ma = k1
c
c----    diagonalize together with bb
         if (me.gt.k1) then
            m0 = ma - 1
            mm = me - m0

            do m1 = ma,me
               do k = 1,n
                  s = zero
                  do i = 1,n
                     s = s + dd(i,m1) * bb(i,k)
                  enddo   ! i
                  z(k) = s
               enddo   ! k
               do m2 = ma,me
                  s = zero
                  do k = 1,n
                     s = s + z(k) * dd(k,m2)
                  enddo   ! k
                  zz(m1-m0,m2-m0) = s
               enddo   ! m2
            enddo   ! m1
            call sdiag(na,mm,zz,eb(ma),zz,z,+1)
            do i = 1,n
               do m = ma,me
                  s = zero
                  do l = ma,me
                     s = s + dd(i,l) * zz(l-m0,m-m0)
                  enddo   ! l
                  z(m) = s
               enddo   ! m
               do m = ma,me
                  dd(i,m) = z(m)
               enddo   ! m
            enddo   ! i
            k1 = me
         endif
      enddo   ! i1
c
   20 return
c-end-DEGEN
      end


c=====================================================================c
      subroutine cmcd
c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
      character nucnam*2

      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /gfvsq / sq(0:igfv)
      common /gfviv / iv(0:igfv)
      common /physco/ hbc,alphi,r0
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4),akacan(nkx,4)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)

      if (icm.ne.2) return

      faccm = -one/two/amu/amas
      ax1 = one/b0/bx
      ax2 = ax1**2
      ay1 = one/b0/by
      ay2 = ay1**2
      az1 = one/b0/bz
      az2 = az1**2

      do it=1,itx
         ecmd(it) = zero
         sumx=zero
         sumy=zero
         sumz=zero
         do ib=1,nb
            nf = id(ib,1)
            do ifg=1,2
               nd = id(ib,ifg)
               i0 = ia(ib,ifg)
               nfg = (ifg-1)*nf
               do i2=1,nd
                  nx2 = nxyz(i0+i2,1)
                  ny2 = nxyz(i0+i2,2)
                  nz2 = nxyz(i0+i2,3)
                  nn2 = nx2+ny2+nz2
                  do i1=1,nd
                     nx1 = nxyz(i0+i1,1)
                     ny1 = nxyz(i0+i1,2)
                     nz1 = nxyz(i0+i1,3)
                     nn1 = nx1+ny1+nz1
                     if(iv(nn1) .eq. iv(nn2)) then
                        k1 = kacan(ib,it)+1
                        k2 = kacan(ib,it)+kdcan(ib,it)
                        rshell = zero
                        do k=k1,k2
                           rshell = rshell + vvcan(k,it)*
     &                       fgcan(nfg+i1,k,it)*fgcan(nfg+i2,k,it)
                        enddo
                        dtx = zero
                        dty = zero
                        dtz = zero
                        if(ny1 .eq. ny2 .and. nz1.eq.nz2) then
                           if(nx1 .eq. nx2-2) dtx = sq(nx2)*sq(nx2-1)
                           if(nx1 .eq. nx2) dtx = -(2*nx2+1)
                           if(nx1 .eq. nx2+2) dtx = sq(nx2+1)*sq(nx2+2)
                           sumx = sumx - dtx*rshell*ax2/two
                        endif
                        if(nx1 .eq. nx2 .and. nz1.eq.nz2) then
                           if(ny1 .eq. ny2-2) dty = sq(ny2)*sq(ny2-1)
                           if(ny1 .eq. ny2)   dty = -(2*ny2+1)
                           if(ny1 .eq. ny2+2) dty = sq(ny2+1)*sq(ny2+2)
                           imy = abs((ny1-ny2)/2)
                           ify = iv(imy)

                           sumy = sumy - dty*ify*rshell*ay2/two
                        endif
                        if(nx1 .eq. nx2 .and. ny1.eq.ny2) then
                           if(nz1 .eq. nz2-2) dtz = sq(nz2)*sq(nz2-1)
                           if(nz1 .eq. nz2)   dtz = -(2*nz2+1)
                           if(nz1 .eq. nz2+2) dtz = sq(nz2+1)*sq(nz2+2)
                           sumz = sumz - dtz*rshell*az2/two
                        endif
                     endif
                  enddo !i1
               enddo !i2
            enddo !ifg
         enddo !ib
         ecmd(it) = two*(sumx+sumy+sumz)*faccm*hbc
      enddo !it
      ecmd(3) = ecmd(1)+ecmd(2)
      return
      end
c=====================================================================c
      subroutine cmcn
c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpx,lpy,lpz
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
      character nucnam*2

      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /gfvsq / sq(0:igfv)
      common /gfviv / iv(0:igfv)
      common /physco/ hbc,alphi,r0
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4),akacan(nkx,4)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)

      if(icm .lt. 2) then
         ecm(1) = -0.75d0*hom
         ecm(2) = -0.75d0*hom
         ecm(3) = -0.75d0*hom
         return
      endif

      faccm = -one/two/amu/amas
      ax1 = one/b0/bx/sq(2)
      ay1 = one/b0/by/sq(2)
      az1 = one/b0/bz/sq(2)

      do it=1,itx
         ecmn(it) = zero
         do ib1=1,nb
            k11 = kacan(ib1,it)+1
            k12 = kacan(ib1,it)+kdcan(ib1,it)
            nf1 = id(ib1,1)
            do ib2=1,nb
               k21 = kacan(ib2,it)+1
               k22 = kacan(ib2,it)+kdcan(ib2,it)
               nf2 = id(ib2,1)
               do k1 = k11,k12
                  vk1 = sqrt(vvcan(k1,it))
                  uk1 = sqrt(abs(one-vk1**2))
                  do k2=k21,k22
                     vk2 = sqrt(vvcan(k2,it))
                     uk2 = sqrt(abs(one-vk2**2))
                     vfac = vk1*vk1*vk2*vk2+vk1*uk1*vk2*uk2

                     termx = zero
                     termy = zero
                     termz = zero

                     do ifg=1,2
                        nd1 = id(ib1,ifg)
                        nd2 = id(ib2,ifg)
                        i01 = ia(ib1,ifg)
                        i02 = ia(ib2,ifg)
                        nfg1 = (ifg-1)*nf1
                        nfg2 = (ifg-1)*nf2
                        do n2=1,nd2
                           nx2 = nxyz(i02+n2,1)
                           ny2 = nxyz(i02+n2,2)
                           nz2 = nxyz(i02+n2,3)
                           nn2 = nx2+ny2+nz2
                           do n1=1,nd1
                              nx1 = nxyz(i01+n1,1)
                              ny1 = nxyz(i01+n1,2)
                              nz1 = nxyz(i01+n1,3)
                              nn1 = nx1+ny1+nz1
                              if(iv(nn1) .ne. iv(nn2)) then
                                 wf1 = fgcan(nfg1+n1,k1,it)
                                 wf2 = fgcan(nfg2+n2,k2,it)
                                 tex = zero
                                 tey = zero
                                 tez = zero
                                 lpx = nx1.eq.nx2
                                 lpy = ny1.eq.ny2
                                 lpz = nz1.eq.nz2

                                 if(lpy .and. lpz) then
                                   if(nx1.eq.nx2-1) tex=iv(nx2)*sq(nx2)
                                   if(nx1.eq.nx2+1) tex=iv(nx1)*sq(nx1)
                                   trfac = iv(ifg+1)
                                   termx = termx + trfac*ax1*
     &                                     iv(ny2+1)*tex*wf1*wf2
                                 endif
                                 if(lpx .and. lpz) then
                                   if(ny1.eq.ny2-1) tey=sq(ny2)
                                   if(ny1.eq.ny2+1) tey=sq(ny1)
                                   trfac = iv(ifg+1)
                                   termy = termy + ay1*
     &                                     trfac*tey*wf1*wf2
                                 endif
                                 if(lpx .and. lpy) then
                                   if(nz1.eq.nz2-1) tez=sq(nz2)
                                   if(nz1.eq.nz2+1) tez=-sq(nz1)
                                   trfac = 1.d0!iv(ifg)
                                   termz = termz + az1*
     &                                     trfac*tez*wf1*wf2
                                 endif
                              endif
                           enddo !n1
                        enddo !n2
                     enddo !ifg
                     f = termx**2+termy**2+termz**2
                     ecmn(it) = ecmn(it)+f*2*vfac
                  enddo
               enddo !k1
            enddo !ib2
         enddo !ib1
         ecmn(it) = -faccm*ecmn(it)*hbc
         ecm(it)  = ecmd(it)+ecmn(it)
      enddo !it
      ecmn(3) = ecmn(1)+ecmn(2)
      ecm(3)  = ecm(1)+ecm(2)
      return
      end

c======================================================================c

      subroutine coulom(lpr)

c======================================================================c
c
c     calculation of the Coulomb field
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /gfvsq / sq(0:igfv)
      common /coulmb/ cou(ngauss),drvp(ngauss)
      common /cou2/ walpot(0:ndcoul,0:ndcoul,nupols,3,0:1)
      common /cou3/ walfac(0:ndcoul,3,0:1)
      common /cou4/ bcoef(0:ndcoul,0:kmax2)
      common /sphpol/ lislam(nupols),lismu(nupols),koeloc(0:4,0:4,1:3),
     &                nreloc(0:4,0:4),kreloc(0:4,0:4,1:8,0:3),
     &                nimloc(0:4,0:4),kimloc(0:4,0:4,1:8,0:3)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      dimension b(3)
      dimension fourx(0:kmax2,0:kmax2)
      dimension foury(0:kmax2)
      dimension four(0:kmax2,0:kmax2,0:kmax2)
      dimension fourc(0:ndcoul,0:ndcoul,0:ndcoul)
      dimension funk(0:ngh,0:ngh,0:ngh)
      dimension work1(0:ndcoul,0:ndcoul),work2(0:ndcoul)
      dimension phi(ngauss)
      dimension rmul(0:4,0:4)

      if (icou.eq.0) return
c
      if (lpr) then
      write(l6,*) '****** BEGIN COULOM ********************************'
      endif
c
c---- Volume term

      call transform(four)

      do ky=0,kmax2
         foury(ky)=0.d0
         do kx=0,kmax2
            fourx(ky,kx) = 0.d0
         enddo
      enddo

      rmul(0,0)=rhomul(0,0)
      rmul(2,0)=rhomul(2,0)
      rmul(2,2)=rhomul(2,2)
      rmul(4,0)=rhomul(4,0)
      rmul(4,2)=rhomul(4,2)
      rmul(4,4)=rhomul(4,4)

      b(1) = one/b0/bx
      b(2) = one/b0/by
      b(3) = one/b0/bz
      facoul = 8/alphi/pi/boucou

      do jx=0,ndcoul
         ipx=mod(jx,2)
         greenx=b(1)*(jx+1)
         greenx = greenx**2
         do kz=0,kmax2
            do ky=0,kmax2
               sum =0.d0
               do kx=ipx,kmax2,2
                  sum=sum+four(kx,ky,kz)*bcoef(jx,kx)
               enddo
               fourx(ky,kz)=sum
            enddo !ky
         enddo !kz
         do jy=0,ndcoul
            ipy=mod(jy,2)
            greeny=b(2)*(jy+1)
            greeny=greeny**2
            do kz=0,kmax2
               sum = 0.d0
               do ky=ipy,kmax2,2
                  sum = sum+fourx(ky,kz)*bcoef(jy,ky)
               enddo
               foury(kz)=sum
            enddo !kz
            do jz=0,ndcoul
               ipz=mod(jz,2)
               greenz=b(3)*(jz+1)
               greenz=greenz**2
               sum=0.d0
               do kz=ipz,kmax2,2
                  sum = sum+foury(kz)*bcoef(jz,kz)
               enddo
               green=greenx+greeny+greenz
               volume=sum/green
c--------------Surface term
               sum =0.d0
               do mulpol=1,nupols
                  lam  = lislam(mulpol)
                  mu   = lismu(mulpol)
                  fmul = rmul(lam,mu)
                  if(lam .ne. 0) fmul=2*fmul

                  facsol   = walpot(jx,jy,mulpol,3,0)*walfac(jz,3,0)
     &                     + walpot(jx,jy,mulpol,3,1)*walfac(jz,3,1)
     &                     + walpot(jz,jx,mulpol,2,0)*walfac(jy,2,0)
     &                     + walpot(jz,jx,mulpol,2,1)*walfac(jy,2,1)
     &                     + walpot(jy,jz,mulpol,1,0)*walfac(jx,1,0)
     &                     + walpot(jy,jz,mulpol,1,1)*walfac(jx,1,1)
                  sum = sum + fmul*facsol
               enddo !mulpol
               surface=sum/green
               fourc(jx,jy,jz) = volume+surface/4/pi
            enddo !jz
         enddo !jy
      enddo !jx
      dx = boucou*bx*b0/sq(2)
      dy = boucou*by*b0/sq(2)
      dz = boucou*bz*b0/sq(2)

      do ihx=0,ngh
         x = xb(ihx)
         do jz=0,ndcoul,2
            do jy=0,ndcoul,2
               sum = 0.d0
               do jx=0,ndcoul,2
                  ax = x*(jx+1)*pi/2/dx
                  fx = cos(ax)
                  sum=sum+fx*fourc(jx,jy,jz)
               enddo
               work1(jy,jz)=sum
            enddo !jy
         enddo !jz
         do ihy=0,ngh
            y=yb(ihy)
            do jz=0,ndcoul,2
               sum = 0.d0
               do jy=0,ndcoul,2
                  ay = y*(jy+1)*pi/2/dy
                  fy = cos(ay)
                  sum= sum + fy*work1(jy,jz)
               enddo !jy
               work2(jz) = sum
            enddo !jz
            do ihz=0,ngh
               z=zb(ihz)
               sum =  0.d0
               do jz=0,ndcoul,2
                  az = z*(jz+1)*pi/2/dz
                  fz = cos(az)
                  sum= sum+fz*work2(jz)
               enddo
               ih=1+ihx + ihy*(ngh+1) + ihz*(ngh+1)*(ngh+1)
               cou(ih)=facoul*sum
            enddo !ihz
         enddo !ihy
      enddo !ihx

      if (lpr) then
      write(l6,*) '****** END COULOM **********************************'
      endif
c
      return
C-end-COULOM
      end
c======================================================================c

      subroutine greecou(lpr)

c======================================================================c
c
C     calculation of the Coulomb-propagator
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr) then
      lx = l6
      l6 = 6
      write(l6,*) '****** BEGIN GREECOU *******************************'
      endif
c
c
      if (icou.eq.0) return

      call calc_fourie()
      call surface_int()
      call calc_b()

      if (lpr) then
      write(l6,*) '****** END GREECOU *********************************'
      l6 = lx
      write(*,*) 'Press enter to continue...'
      read(*,*)
      endif
      return
c-end-GREECOU
      end
c======================================================================
      subroutine calc_b()
c======================================================================
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      common /mathco/ zero,one,two,half,third,pi
      common /cou4/ bcoef(0:ndcoul,0:kmax2)
      common /gfviv / iv(0:igfv)

      do j=0,ndcoul
         omega = (j+1)*pi/2/boucou
         efac = sqrt(2*pi)*exp(-omega*omega/2)
         do k=0,kmax2
            if(mod(k+j,2).eq.0) then
               res = efac*hermit(k,omega)*iv(k/2)
            else
               res = zero
            endif
            bcoef(j,k) = res
         enddo
      enddo

      return
      end
c======================================================================c
      subroutine calc_fourie()
c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      common /mathco/ zero,one,two,half,third,pi
      common /cou1/ fourie(0:ndcoul,1:ndneta)
c
      stepom = pi/(2*boucou)
      stepet = 2*boucou/(ndneta+1)

      do ifouri=0,ndcoul
         xfouri = stepom*(ifouri+1)
         if(mod(ifouri,2).eq.0) then
            do indeta=1,ndneta
               eta = -boucou+stepet*indeta
               fourie(ifouri,indeta) = cos(xfouri*eta)
            enddo
         else
            do indeta=1,ndneta
               eta = -boucou+stepet*indeta
               fourie(ifouri,indeta) = sin(xfouri*eta)
            enddo
         endif
      enddo
      do ifouri=0,ndcoul
         do indeta=1,ndneta
            eta = -boucou+stepet*indeta
         enddo
      enddo
      return
      end
c======================================================================c
      subroutine surface_int()
c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /mathco/ zero,one,two,half,third,pi
      common /gfvsq / sq(0:igfv)
      common /gfviv / iv(0:igfv)
      common /cou1/ fourie(0:ndcoul,1:ndneta)
      common /cou2/ walpot(0:ndcoul,0:ndcoul,nupols,3,0:1)
      common /cou3/ walfac(0:ndcoul,3,0:1)
      common /sphpol/ lislam(nupols),lismu(nupols),koeloc(0:4,0:4,1:3),
     &                nreloc(0:4,0:4),kreloc(0:4,0:4,1:8,0:3),
     &                nimloc(0:4,0:4),kimloc(0:4,0:4,1:8,0:3)
      dimension b(3),poskar(3)

      b(1)=one/(b0*bx)
      b(2)=one/(b0*by)
      b(3)=one/(b0*bz)

      stepet = 2*boucou/(ndneta+1)

      do kartez=1,ndkart
         kartex=kartez-2
         if(kartex .le. 0) kartex=kartex+3
         kartey=kartez-1
         if(kartey .le. 0) kartey=kartey+3
         fnorm = one/4/boucou
         do iwside=0,1
            poskar(kartez) = (2*iwside-1)*boucou/b(kartez)/sq(2)
            do mulpol=1,nupols
               lambda = lislam(mulpol)
               mu = lismu(mulpol)
               do jx=0,ndcoul
                  do jy=0,ndcoul
                     walpot(jx,jy,mulpol,kartez,iwside) = zero
                  enddo
               enddo !jx

               do indetx=1,ndneta
                  ifx = 2*(mod(indetx,2)+1)
                  eta = -boucou+stepet*indetx
                  poskar(kartex)=eta/b(kartex)/sq(2)
                  do indety=1,ndneta
                     ify = 2*(mod(indety,2)+1)
                     eta = -boucou+stepet*indety
                     poskar(kartey)=eta/b(kartey)/sq(2)
                     rehar = solhar(lambda,mu,poskar)
                     p1 = poskar(1)**2
                     p2 = poskar(2)**2
                     p3 = poskar(3)**2
                     diswal = sqrt(p1+p2+p3)
                     rehar = rehar/diswal**lambda
                     potwal=4*pi*rehar/diswal**(lambda+1)/(2*lambda+1)
                     potwal=ifx*ify*potwal
                     do jx=0,ndcoul
                        do jy=0,ndcoul
                           walpot(jx,jy,mulpol,kartez,iwside) =
     &                      walpot(jx,jy,mulpol,kartez,iwside)+
     &                      potwal*fourie(jx,indetx)*fourie(jy,indety)
                        enddo !jy
                     enddo !jx
                  enddo !indety
               enddo !indetx
               do jx=0,ndcoul
                  do jy=0,ndcoul
                     fac = fnorm*stepet*stepet/3/3
                     walpot(jx,jy,mulpol,kartez,iwside)
     &                      =walpot(jx,jy,mulpol,kartez,iwside)*fac
                  enddo !jy
               enddo !jx
            enddo !mulpol
            do jz=0,ndcoul
               factor = b(kartez)*b(kartez)*iv(jz/2)*(jz+1)
               if(iwside .eq. 0 .and. mod(jz,2).eq.1) factor=-factor
               walfac(jz,kartez,iwside)=factor
            enddo
         enddo !iwside
      enddo !kartez
      return
      end
c======================================================================c
      real*8 function solhar(lambda,mu,poskar)
c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      common /mathco/ zero,one,two,half,third,pi
      common /sphpol/ lislam(nupols),lismu(nupols),koeloc(0:4,0:4,1:3),
     &                nreloc(0:4,0:4),kreloc(0:4,0:4,1:8,0:3),
     &                nimloc(0:4,0:4),kimloc(0:4,0:4,1:8,0:3)
      dimension poskar(3)

      nterm=nreloc(lambda,mu)
      sum = zero
      do iterm=1,nterm
         fac = dble(kreloc(lambda,mu,iterm,0))

         do kart=1,3
            x=poskar(kart)
            ipow=kreloc(lambda,mu,iterm,kart)
            fac = fac*x**ipow
         enddo
         sum = sum + fac
      enddo
      C1= dble(koeloc(lambda,mu,1))
      C2=dble(koeloc(lambda,mu,2))
      C3=dble(koeloc(lambda,mu,3))
      C = C1*sqrt(C2/C3)

      solhar = C*sum

      return
      end
c======================================================================c
      real*8 function hermit(k,x)
c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'

      dimension u(0:kmax2),v(0:kmax2)
      common /mathco/ zero,one,two,half,third,pi
      common /gfvsq / sq(0:igfv)
      common /gfviv / iv(0:igfv)
      common /gfvsqi/ sqi(0:igfv)

      w4pi = pi**(-0.25d0)
      s2 = sqi(2)
      x2 = x*x
      u(0)=w4pi
      v(0)=-w4pi*x

      do n=1,k
         s = s2*sqi(n)
         u(n) = s*(x*u(n-1)-v(n-1))
         v(n) = s*(x*v(n-1)+ (dble(2*n)-x2)*u(n-1))
      enddo
      hermit = u(k)

      return
      end
c==========================================================================
      subroutine transform(four)
c=========================================================================
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      dimension four(0:kmax2,0:kmax2,0:kmax2)
      common /mathco/ zero,one,two,half,third,pi
      common /gaussh/  xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /hermite2/ hpol2(0:kmax2,0:ngh),dhpol2(0:kmax2,0:ngh)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
      common /gfvsq / sq(0:igfv)
      dimension f1(0:ngh,0:ngh),f2(0:ngh)

      do kz=0,kmax2
         do ky=0,kmax2
            do kx=0,kmax2
               four(kx,ky,kz)=0.d0
            enddo
         enddo
      enddo

      do kz=0,kmax2,2
         do ihx=0,ngh
            do ihy=0,ngh
               sum=0.d0
               do ihz=0,ngh
                   ih=1+ihx + ihy*(ngh+1) + ihz*(ngh+1)*(ngh+1)
                   sum = sum+hpol2(kz,ihz)*ph(ihz)*rv(ih,2)
               enddo !ihz
               f1(ihy,ihx) = sum
            enddo !ihy
         enddo !ihx

         do ky=0,kmax2,2
            do ihx=0,ngh
               sum=0.d0
               do ihy=0,ngh
                  sum = sum+hpol2(ky,ihy)*ph(ihy)*f1(ihy,ihx)
               enddo
               f2(ihx) = sum
            enddo !ihx
            do kx=0,kmax2,2
               sum=0.d0
               do ihx=0,ngh
                  sum = sum + hpol2(kx,ihx)*ph(ihx)*f2(ihx)
               enddo
               four(kx,ky,kz) = sq(2)**3*sum
            enddo !kx
         enddo !ky
      enddo !kz
      return
      end
c==============================================================================
      real*8 function rhomul(lam,mu)
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      common /mathco/ zero,one,two,half,third,pi
      common /sphpol/ lislam(nupols),lismu(nupols),koeloc(0:4,0:4,1:3),
     &                nreloc(0:4,0:4),kreloc(0:4,0:4,1:8,0:3),
     &                nimloc(0:4,0:4),kimloc(0:4,0:4,1:8,0:3)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
      common /gaucor/ wdcor(ngauss)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)

      sum = 0.d0
      do ihx=0,ngh
         x = xb(ihx)
         do ihy=0,ngh
            y = yb(ihy)
            do ihz=0,ngh
               z = zb(ihz)
               ih=1+ihx + ihy*(ngh+1) + ihz*(ngh+1)*(ngh+1)
               nterm = nreloc(lam,mu)
               do n=1,nterm
                  ipx = kreloc(lam,mu,n,1)
                  ipy = kreloc(lam,mu,n,2)
                  ipz = kreloc(lam,mu,n,3)
                  fac = kreloc(lam,mu,n,0)
                  fx  = x**ipx
                  fy  = y**ipy
                  fz  = z**ipz
                  sum = sum + fac*fx*fy*fz*wdcor(ih)*rv(ih,2)
               enddo !n
            enddo  !ihz
         enddo  !ihy
      enddo !ihx
      c1 = dble(koeloc(lam,mu,1))
      c2 = dble(koeloc(lam,mu,2))
      c3 = dble(koeloc(lam,mu,3))
      c  = c1*sqrt(c2/c3)
      rhomul = c*sum

      return
      end

c======================================================================c

      subroutine cstrpot(lpr)

c======================================================================c
c
c     CALCULATION OF THE CONSTRAINING POTENT
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      common /baspar/ hom,hb0,b0
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /cstr1/ betac,gammac,cqad,icstr
      common /cstr2/ q0c,q2c,c0,c2,alaq0,alaq2
      common /cstr3/ vc(ngauss)
      common /cstr4/ calcxx,calcyy,calczz
      common /cstr5/ fac0,fac2
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr) then
      write(l6,*) '****** BEGIN CSTRPOT *******************************'
      endif
c
      fac0=zero
      fac2=zero
      if(icstr.eq.0) then
         do ih=1,ngauss
            vc(ih)=zero
         enddo
         return
      endif

      calc0 = 2*calczz-calcxx-calcyy
      calc2 = calcxx-calcyy
      alaq0=alaq0+c0*(q0c-calc0)
      alaq2=alaq2+c2*(q2c-calc2)

      q0cst = q0c+alaq0/c0
      q2cst = q2c+alaq2/c2

      fac0 = -c0*(q0cst-calc0)
      fac2 = -c2*(q2cst-calc2)
c
      do iz = 0,ngh
         zz = zb(iz)**2
      do iy = 0,ngh
         yy = yb(iy)**2
      do ix = 0,ngh
         xx = xb(ix)**2
         f1 = fac0*(2*zz-xx-yy)
         f2 = fac2*(xx-yy)
         ih=1+ix+iy*(ngh+1)+iz*(ngh+1)*(ngh+1)
         vc(ih) = f1+f2
      enddo   ! ix
      enddo   ! iy
      enddo   ! iz

      if (lpr) then
      write(l6,*) '****** END CSTRPOT *********************************'
      endif
c
      return
c-end-CSTRPOT
      end

c======================================================================c

      subroutine default()

c======================================================================c
c
c     Default data
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
      common /baspar/ hom,hb0,b0
      common /couplf/ ff(ngauss,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /fermi / ala(2),tz(2)
      common /pair  / del(2),spk(2),spk0(2)
      common /cstr2/ q0c,q2c,c0,c2,alaq0,alaq2
      common /initia/ vin,rin,ain,inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /broyde2/ ibroyd
c
c======================================================================c
c---- signs, factorials and mathematical constants
c-----------------------------------------------------------------------
      call gfv()
c======================================================================c
c---- Neutrons and Protons
c-----------------------------------------------------------------------
c     number of isospins
c     itx = 1   only neutrons
c           2   neutrons and protons
c
      itx = 2
c
c======================================================================c



c======================================================================c
c     center off mass correction:
c-----------------------------------------------------------------------
c
c     a) the kinetic energy in the variation can be
c           <T>  = hbc**2/(2*amu)*Delta   (no correction)
c        or <T>' = hbc**2/(2*amu)*(1-1/A)*Delta
c                  the diagonal part of ecm=<P**2>/2Am is included
c     b) the total energy is
c        <T> + <V> - <P**2>/2Am
c
c     icm: 0   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = 3/4*hom
c          1   variation of <T>'
c              hb0 = hb0*(1-1/A)          ecm = 3/4*hom
c          2   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = <P**2>/2M
c
      icm = 0
c
c======================================================================c
c     Coulomb-Field
c-----------------------------------------------------------------------
c     icou: Coulomb-field:  0   not at all
c                           1   only direct term
c                           2   plus exchange
      icou = 1
c     write(6,*) 'without Coulomb: icou =',icou
c     read(*,*)
c
c======================================================================c


c======================================================================c
c     point-coupling and meson-coupling
c----------------------------------------------------------------------c
c     ipc = 0:        meson-coupling
c           1:        point-coupling
      ipc = 0
c
c
c     density-dependence of the coupling constants:
c                           0  none
c                           1  linear density dependence
c                           2  form of Typel and Wolter
      idd=0
c
c======================================================================c
c---- Fermi level
      do it = 1,2
	     ala(it)  = -7.d0
	     spk0(it) =  0.d0
      enddo   ! it
      alaq0 = 0.d0
      alaq2 = 0.d0
c======================================================================c
c     iteration
c----------------------------------------------------------------------c
c
      maxi = 200               ! maximal number of iteration
      si   = one             ! actual error in the main iteration
      epsi = 1.d-6       ! accuracy for the main iteration
      iaut  = 1              ! automatic change of xmix: 0 (no) 1 (yes)
      inxt  = 200            ! number of iterations until next question
      xmix  = 0.5d0          ! starting xmix
      xmax  = 0.7d0          ! maximal value for xmix
      xmix0 = xmix           ! basic xmix, where it returns to
      ibroyd= 1
c======================================================================c
c
c
c---- parameters of the initial potentials
c     inin = 0: fields read, 1: saxon-wood,
      inin  = 1

c     oscillator length b0 (is calcuated for b0 <= 0)
      b0 = -2.320
c======================================================================c
c---- tapes
c----------------------------------------------------------------------c
      l6   = 10
      lin  = 3
      lou  = 6
      lwin = 1
      lwou = 2
      lplo = 11
      laka = 12
c======================================================================c
c---- preparation of density dependence
c----------------------------------------------------------------------c
      do m = 1,4
         lmes(m) = 0
         ggmes(m) = zero
         do i = 1,ngauss
            ff(i,m,1) = one
	    ff(i,m,2) = zero
         enddo   ! ngauss
      enddo   ! m
c
      return
c-end-DEFAULT
      end
c======================================================================c
      blockdata block1
c======================================================================c
      implicit real*8 (a-h,o-z)

      include 'dirhb.par'
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      common /physco/ hbc,alphi,r0
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)

c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutrons','Protons '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
c---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/

c-end-BLOCK1
      end
c======================================================================c
      blockdata block2
c======================================================================c
      implicit real*8 (a-h,o-z)

      include 'dirhb.par'
      common /sphpol/ lislam(nupols),lismu(nupols),koeloc(0:4,0:4,1:3),
     &                nreloc(0:4,0:4),kreloc(0:4,0:4,1:8,0:3),
     &                nimloc(0:4,0:4),kimloc(0:4,0:4,1:8,0:3)
c
c======================================================================c
c---- parameters for the spherical harmonics
c-----------------------------------------------------------------------
c    Y00 = C                   C=c1*sqrt(c2/c3)/sqrt(4*pi)
c    C=1/sqrt(4*pi)
c    lislam(i) : lambda for the i-th spherical harmonic
c    lismu(i)  : mu for the i-th spherical harmonic
c    koeloc(lambda,mu,k): k=1, c1
c                         k=2, c2
c                         k=3, c3
c    nreloc(lambda,mu): number of real terms
c    kreloc(lambda,mu,i,K=0) : coefficient of the term
c    kreloc(lambda,mu,i,K=1) : exponent of the x part
c    kreloc(lambda,mu,i,K=2) : exponent of the y part
c    kreloc(lambda,mu,i,K=3) : exponent of the z part
c-----------------------------------------------------------------------
      data lislam(1) /0/, lismu(1) /0/
      data (koeloc(0,0,K),K=1,3) /1,1,1/
      data nreloc(0,0) /1/
      data (kreloc(0,0,1,K),K=0,3) /1,0,0,0/
      data nimloc(0,0) /0/
c-----------------------------------------------------------------------
c    Y20 = C*(2z^2-x^2-y^2)
c    C=sqrt(5/4)/sqrt(4*pi)
c-----------------------------------------------------------------------
      data lislam(2) /2/, lismu(2) /0/
      data (koeloc(2,0,K),K=1,3) /1,5,4/
      data nreloc(2,0) /3/
      data (kreloc(2,0,1,K),K=0,3) /2,0,0,2/
      data (kreloc(2,0,2,K),K=0,3) /-1,2,0,0/
      data (kreloc(2,0,3,K),K=0,3) /-1,0,2,0/
      data nimloc(2,0) /0/
c-----------------------------------------------------------------------
c    Y22 = C*(x^2-2*I*x*y-y^2)
c    C=sqrt(15/8)/sqrt(4*pi)
c-----------------------------------------------------------------------
      data lislam(3) /2/, lismu(3) /2/
      data (koeloc(2,2,K),K=1,3) /1,15,8/
      data nreloc(2,2) /2/
      data (kreloc(2,2,1,K),K=0,3) / 1,2,0,0/
      data (kreloc(2,2,2,K),K=0,3) /-1,0,2,0/
      data nimloc(2,2) /1/
      data (kimloc(2,2,1,K),K=0,3) /-2,1,1,0/
c-----------------------------------------------------------------------
c    Y40 = C*(8*z^4+3*x^4+3*y^4-24*x^2*z^2-24*y^2*z^2+6*x^2*y^2)
c    C=3/8/sqrt(4*pi)
c-----------------------------------------------------------------------
      data lislam(4) /4/, lismu(4) /0/
      data (koeloc(4,0,K),K=1,3) /1,9,64/
      data nreloc(4,0) /6/
      data (kreloc(4,0,1,K),K=0,3) / 8,0,0,4/
      data (kreloc(4,0,2,K),K=0,3) / 3,4,0,0/
      data (kreloc(4,0,3,K),K=0,3) / 3,0,4,0/
      data (kreloc(4,0,4,K),K=0,3) /-24,2,0,2/
      data (kreloc(4,0,5,K),K=0,3) /-24,0,2,2/
      data (kreloc(4,0,6,K),K=0,3) / 6,2,2,0/
      data nimloc(4,0) /0/
c-----------------------------------------------------------------------
c    Y42 = C*(-x^4+y^4+6z^2x^2-6z^2y^2)+C*i*(12xyz^2-2x^3y-2xy^3)
c    C=sqrt(90/64)/sqrt(4*pi)
c-----------------------------------------------------------------------
      data lislam(5) /4/, lismu(5) /2/
      data (koeloc(4,2,K),K=1,3) /1,90,64/
      data nreloc(4,2) /4/
      data (kreloc(4,2,1,K),K=0,3) /-1,4,0,0/
      data (kreloc(4,2,2,K),K=0,3) / 1,0,4,0/
      data (kreloc(4,2,3,K),K=0,3) / 6,2,0,2/
      data (kreloc(4,2,4,K),K=0,3) /-6,0,2,2/
      data nimloc(4,2) /3/
      data (kimloc(4,2,1,K),K=0,3) / 12,1,1,2/
      data (kimloc(4,2,2,K),K=0,3) / -2,3,1,0/
      data (kimloc(4,2,3,K),K=0,3) / -2,1,3,0/
c-----------------------------------------------------------------------
c    Y44 = C*(x^4-6x^2y^2+y^4)+C*i*(4x^3y-4xy^3)
c    C=sqrt(630/256)/sqrt(4*pi)
c-----------------------------------------------------------------------
      data lislam(6) /4/, lismu(6) /4/
      data (koeloc(4,4,K),K=1,3) /1,630,256/
      data nreloc(4,4) /3/
      data (kreloc(4,4,1,K),K=0,3) / 1,4,0,0/
      data (kreloc(4,4,2,K),K=0,3) /-6,2,2,0/
      data (kreloc(4,4,3,K),K=0,3) / 1,0,4,0/
      data nimloc(4,4) /2/
      data (kimloc(4,4,1,K),K=0,3) / 4,3,1,0/
      data (kimloc(4,4,2,K),K=0,3) /-4,1,3,0/
c-end-BLOCK2
      end

c======================================================================c

      subroutine delta(it,lpr)

c======================================================================c
c
c     calculates the pairing field
c     IT  = 1 neutrons
c           2 protons
c
c     units:    fields and Hamiltonian in MeV
c               eigenvalues in MeV
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      character tb*5                                            !blokap
      character tt*8                                            ! bloqua
      logical lpr,lpx,lpy,lpz

      dimension work(nsepx)
      common /baspar/ hom,hb0,b0
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
      common /vvvsep/ vnn(0:n0fx,0:n0fx,0:npmax/2,3)
      common /deldel/ de(nhhx,nb2x)
      common /tmrpar/ gl(2),gal
      common /gfviv / iv(0:igfv)
      common /numsep/ npm
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut

c
      if (lpr) then
      write(l6,*) '****** BEGIN DELTA***************************'
      endif

      i12 = 0
      do isep=1,nsepx
         work(isep) = 0.d0
      enddo
      do ib=1,nb
         i0 = ia(ib,1)
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf+ng
         m  = ib+(it-1)*nbx
         do i1=1,nh
            do i2=1,nh
               de(i2+(i1-1)*nh,m)=0.d0
            enddo
         enddo
         do i2=1,nf
             nx2=nxyz(i0+i2,1)
             ny2=nxyz(i0+i2,2)
             nz2=nxyz(i0+i2,3)
             do i1=i2,nf
                nx1=nxyz(i0+i1,1)
                ny1=nxyz(i0+i1,2)
                nz1=nxyz(i0+i1,3)

                nx = nx1+nx2
                ny = ny1+ny2
                nz = nz1+nz2

                i12 = i12 + 1
                if (iv(nx).eq.1.and.iv(ny).eq.1.and.iv(nz).eq.1) then
                isep = 0
                do ix = 0,npm,2
                   vx = vnn(nx1,nx2,ix/2,1)
                   do iy = 0,npm-ix,2
                      vy = vnn(ny1,ny2,iy/2,2)
                      do iz = 0,npm-ix-iy,2
                         isep = isep + 1
                         vz = vnn(nz1,nz2,iz/2,3)
                         ifac=iv(abs(ny1-ny2)/2)
                         work(isep) = work(isep)
     &                      +ifac*vx*vy*vz*aka(i12,it)
                      enddo !iz
                   enddo !iy
                enddo !ix
                endif
             enddo !i1
          enddo !i2
      enddo !ib

      do ib=1,nb
         i0 = ia(ib,1)
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf+ng
         m  = ib+(it-1)*nbx
         do i4=1,nf
             nx4=nxyz(i0+i4,1)
             ny4=nxyz(i0+i4,2)
             nz4=nxyz(i0+i4,3)
             do i3=i4,nf
                nx3=nxyz(i0+i3,1)
                ny3=nxyz(i0+i3,2)
                nz3=nxyz(i0+i3,3)
                nx=nx3+nx4
                ny=ny3+ny4
                nz=nz3+nz4
                lpx = iv(nx).eq.1
                lpy = iv(ny).eq.1
                lpz = iv(nz).eq.1
                if(lpx.and.lpy.and.lpz) then
                   isep = 0
                   s = 0.d0
                   do ix = 0,npm,2
                      vx = vnn(nx3,nx4,ix/2,1)
                      do iy = 0,npm-ix,2
                         vy = vnn(ny3,ny4,iy/2,2)
                         do iz = 0,npm-ix-iy,2
                            isep = isep + 1
                            vz = vnn(nz3,nz4,iz/2,3)
                            ifac=iv(abs(ny3-ny4)/2)
                            s = s + ifac*vx*vy*vz*work(isep)
                         enddo !iz
                      enddo !iy
                   enddo !ix
                   s = -0.5d0*s*gl(it)
              else
                s=0.d0
              endif
                de(i3+(i4-1)*nh,m) = s
                de(i4+(i3-1)*nh,m) = s
             enddo !i3
          enddo !i4
      enddo !ib
c
      if (lpr) then
      write(l6,*) '****** END DELTA ******************************'
      endif
c
      return
C-end-DELTA
      end

c=====================================================================c
      subroutine densit(lpr)
c=====================================================================c
c     calculates the densities in r-space at Gauss-meshpoints
c---------1---------2---------3---------4---------5---------6---------7-
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      dimension rshell(2)
      dimension frsn(0:kmax2,0:kmax2,0:kmax2)
      dimension frsp(0:kmax2,0:kmax2,0:kmax2)
      dimension frvn(0:kmax2,0:kmax2,0:kmax2)
      dimension frvp(0:kmax2,0:kmax2,0:kmax2)
      dimension tmp(0:ngh,0:ngh,0:ngh)
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /gaucor/ wdcor(ngauss)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gfviv / iv(0:igfv)
      common /herexp/  C00(0:kmax2,0:kmax,0:kmax)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
c
      if (lpr) then
      write(l6,*) '****** BEGIN DENSIT ********************************'
      endif
c
      do i = 1,ngauss
         do it = 1,2
            rs(i,it)   = zero
            rv(i,it)   = zero
            drs(i,it)  = zero
            drv(i,it)  = zero
            rshell(it) = zero
         enddo
      enddo
      do kx=0,kmax2
         do ky=0,kmax2
            do kz=0,kmax2
               frsn(kx,ky,kz) = 0.d0
               frsp(kx,ky,kz) = 0.d0
               frvn(kx,ky,kz) = 0.d0
               frvp(kx,ky,kz) = 0.d0
            enddo
         enddo
      enddo
c
c---- loop over simplex-parity-blocks
      do ib = 1,nb
c
c------- loop over contributions from large and small components
         do ifg = 1,2
            nd  = id(ib,ifg)
            nh  = id(ib,1)+id(ib,2)
            i0  = ia(ib,ifg)
            nfg = (ifg-1)*id(ib,1)
            ivv = iv(ifg)
c
c---------- loop of oscillator basis states n2 and n1
            do n2 =  1,nd
               nx2 = nxyz(i0+n2,1)
               ny2 = nxyz(i0+n2,2)
               nz2 = nxyz(i0+n2,3)
            do n1 = n2,nd
               nx1 = nxyz(i0+n1,1)
               ny1 = nxyz(i0+n1,2)
               nz1 = nxyz(i0+n1,3)
               nx = nx1+nx2
               ny = ny1+ny2
               nz = nz1+nz2
c
               if (iv(nx).eq.1.and.iv(ny).eq.1.and.iv(nz).eq.1) then
                  i12 = (2 - n2/n1)*iv(iabs(ny1-ny2)/2)
c---------- loop over neutrons and protons
                  do it = 1,itx
                     m=ib+(it-1)*nbx
                     rshell(it) = i12*rosh(nfg+n1+(nfg+n2-1)*nh,m)
                  enddo   ! it
                  do kz=0,nz,2
                     fz = C00(kz,nz1,nz2)
                     do ky=0,ny,2
                        fy = C00(ky,ny1,ny2)
                        do kx=0,nx,2
                           fx = C00(kx,nx1,nx2)
                           frvn(kx,ky,kz) = frvn(kx,ky,kz)
     &                                     +rshell(1)*fx*fy*fz
                           frvp(kx,ky,kz) = frvp(kx,ky,kz)
     &                                     +rshell(2)*fx*fy*fz
                           frsn(kx,ky,kz) = frsn(kx,ky,kz)
     &                                     -ivv*rshell(1)*fx*fy*fz
                           frsp(kx,ky,kz) = frsp(kx,ky,kz)
     &                                     -ivv*rshell(2)*fx*fy*fz
                        enddo !kx
                     enddo !ky
                  enddo !kz
               endif

c---------- end loop of oscillator basis states n2 and n1
            enddo   ! n1
            enddo   ! n2
c
c------- end loop over large and small components
         enddo   ! ifg
c
c---- end loop over blocks ib
      enddo   ! ib
      call invden(frvn,1)
      call invden(frvp,2)
      call invden(frsn,3)
      call invden(frsp,4)
      call invlap(frvn,1)
      call invlap(frvp,2)
      call invlap(frsn,3)
      call invlap(frsp,4)
c---- check, whether integral over dro vanishes
      s1 = zero
      s2 = zero
      do i = 1,ngauss
         s1 = s1 + drv(i,1)*wdcor(i)
         s2 = s2 + drv(i,2)*wdcor(i)
      enddo
      if (lpr) write(l6,*) ' integral over dro',s1,s2
c
c---- normalization and renormalization to particle number
      do it = 1,itx
         s  = zero
         do i = 1,ngauss
            s  =  s + rv(i,it)*wdcor(i)
         enddo
         if (lpr) write(l6,'(a,i3,2f15.8)') '  Integral over rv:',it,s

      enddo  ! it
c
      do i = 1,ngauss
	    ro(i,1)  = rs(i,1) + rs(i,2)
	    ro(i,2)  = rv(i,1) + rv(i,2)
	    ro(i,3)  =-rs(i,1) + rs(i,2)
	    ro(i,4)  =-rv(i,1) + rv(i,2)
	    dro(i,1) = drs(i,1) + drs(i,2)
	    dro(i,2) = drv(i,1) + drv(i,2)
	    dro(i,3) =-drs(i,1) + drs(i,2)
	    dro(i,4) =-drv(i,1) + drv(i,2)
      enddo   ! i

      if (lpr) then
      write(l6,*) '****** END DENSIT **********************************'
      endif
c
      return
C-end-DENSIT
      end
c=============================================================================
c
      subroutine invden(four,iden)
c
c=============================================================================
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      common /hermite/ hpol(0:kmax2,0:ngh),dhpol(0:kmax2,0:ngh),
     &                 ddhpol(0:kmax2,0:ngh)
      common /baspar/ hom,hb0,b0
      common /basdef/ beta0,gamma0,bx,by,bz
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
      common /mathco/ zero,one,two,half,third,pi

      dimension four(0:kmax2,0:kmax2,0:kmax2)
      dimension tmp1(0:kmax2,0:kmax2)
      dimension tmp2(0:kmax2)
      do ihx=0,ngh
         xx = xh(ihx)*xh(ihx)
         fx = exp(-xx)/bx/b0
         do kz=0,kmax2,2
            do ky=0,kmax2,2
               sum=0.d0
               do kx=0,kmax2,2
                  sum = sum + four(kx,ky,kz)*hpol(kx,ihx)
               enddo
               tmp1(ky,kz) = sum
            enddo !ky
         enddo !kz
         do ihy=0,ngh
            yy = xh(ihy)*xh(ihy)
            fy = exp(-yy)/by/b0
            do kz=0,kmax2,2
               sum = zero
               do ky=0,kmax2,2
                  sum = sum+tmp1(ky,kz)*hpol(ky,ihy)
               enddo
               tmp2(kz)=sum
            enddo !kz
            do ihz=0,ngh
               zz = xh(ihz)*xh(ihz)
               fz = exp(-zz)/bz/b0
               sum=0.d0
               do kz=0,kmax2,2
                  sum=sum+tmp2(kz)*hpol(kz,ihz)
               enddo
               i = 1+ihx+ihy*(ngh+1)+ihz*(ngh+1)*(ngh+1)
               if (iden.eq.1) then
                  rv(i,1)=sum*fx*fy*fz
               else if(iden .eq. 2) then
                  rv(i,2)=sum*fx*fy*fz
               else if(iden .eq. 3) then
                  rs(i,1)=sum*fx*fy*fz
               else
                  rs(i,2)=sum*fx*fy*fz
               endif
            enddo !ihz
         enddo !ihy
      enddo !ihx


      return
      end
c=============================================================================
c
      subroutine invlap(four,ilap)
c
c=============================================================================
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
      common /hermite/ hpol(0:kmax2,0:ngh),dhpol(0:kmax2,0:ngh),
     &                 ddhpol(0:kmax2,0:ngh)
      common /baspar/ hom,hb0,b0
      common /basdef/ beta0,gamma0,bx,by,bz
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /rhorho  / rs(ngauss,2),rv(ngauss,2),
     &                  drs(ngauss,2),drv(ngauss,2)
      common /mathco/ zero,one,two,half,third,pi
      dimension four(0:kmax2,0:kmax2,0:kmax2)
      dimension tmp1x(0:kmax2,0:kmax2)
      dimension tmp2x(0:kmax2)
      dimension tmp1y(0:kmax2,0:kmax2)
      dimension tmp2y(0:kmax2)
      dimension tmp1z(0:kmax2,0:kmax2)
      dimension tmp2z(0:kmax2)

      ax = 1.d0/(b0*bx)**2
      ay = 1.d0/(b0*by)**2
      az = 1.d0/(b0*bz)**2

      do ihx=0,ngh
         xx = xh(ihx)*xh(ihx)
         fx = exp(-xx)/bx/b0
         do kz=0,kmax2,2
            do ky=0,kmax2,2
               sumx=0.d0
               sumy=0.d0
               sumz=0.d0
               do kx=0,kmax2,2
                  sumx = sumx + four(kx,ky,kz)*ddhpol(kx,ihx)
                  sumy = sumy + four(kx,ky,kz)*hpol(kx,ihx)
                  sumz = sumz + four(kx,ky,kz)*hpol(kx,ihx)
               enddo
               tmp1x(ky,kz) = sumx
               tmp1y(ky,kz) = sumy
               tmp1z(ky,kz) = sumz
            enddo !ky
         enddo !kz
         do ihy=0,ngh
            yy = xh(ihy)*xh(ihy)
            fy = exp(-yy)/by/b0
            do kz=0,kmax2,2
               sumx = 0.d0
               sumy = 0.d0
               sumz = 0.d0
               do ky=0,kmax2,2
                  sumx = sumx+tmp1x(ky,kz)*hpol(ky,ihy)
                  sumy = sumy+tmp1y(ky,kz)*ddhpol(ky,ihy)
                  sumz = sumz+tmp1z(ky,kz)*hpol(ky,ihy)
               enddo
               tmp2x(kz)=sumx
               tmp2y(kz)=sumy
               tmp2z(kz)=sumz
            enddo !kz
            do ihz=0,ngh
               zz = xh(ihz)*xh(ihz)
               fz = exp(-zz)/bz/b0
               sumx=0.d0
               sumy=0.d0
               sumz=0.d0
               do kz=0,kmax2,2
                  sumx = sumx+tmp2x(kz)*hpol(kz,ihz)
                  sumy = sumy+tmp2y(kz)*hpol(kz,ihz)
                  sumz = sumz+tmp2z(kz)*ddhpol(kz,ihz)
               enddo
               i = 1+ihx+ihy*(ngh+1)+ihz*(ngh+1)*(ngh+1)
               if (ilap.eq.1) then
                  drv(i,1)=(sumx*ax+sumy*ay+sumz*az)*fx*fy*fz
               else if(ilap .eq. 2) then
                  drv(i,2)=(sumx*ax+sumy*ay+sumz*az)*fx*fy*fz
               else if(ilap .eq. 3) then
                  drs(i,1)=(sumx*ax+sumy*ay+sumz*az)*fx*fy*fz
               else
                  drs(i,2)=(sumx*ax+sumy*ay+sumz*az)*fx*fy*fz
               endif
            enddo !ihz
         enddo !ihy
      enddo !ihx


      return
      end

c=====================================================================c
      subroutine denssh(it,lpr)
c=====================================================================c
c     calculates the densities in h.o. space
c---------1---------2---------3---------4---------5---------6---------7-
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
      character tb*5                                            ! blokap
c
c
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
      common /pair  / del(2),spk(2),spk0(2)
c
      if (lpr) then
      write(l6,*) '****** BEGIN DENSSH ********************************'
      endif
c
      spk(it) = zero
c
c---- loop over simplex-parity-blocks
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf+ng
         k1 = ka(ib,it) + 1
         ke = ka(ib,it) + kd(ib,it)
         k1a = ka(ib,it+2) + 1
         kea = ka(ib,it+2) + kd(ib,it+2)
         m=ib+(it-1)*nbx
c
c---------- loop of oscillator basis states n2 and n1
         do n2 =  1,nh
            do n1 = n2,nh
c---------- loop over neutrons and protons
               s = zero
               do k=k1,ke
                  s=s+fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
               enddo    ! k
               do k=k1a,kea
                  s=s+fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
               enddo    ! k
               rosh(n1+(n2-1)*nh,m) = 2*s
               rosh(n2+(n1-1)*nh,m) = 2*s
            enddo
         enddo
         do n2=1,nf
            do n1=n2,nf
               i12 = (2 - n2/n1)
               il=il+1
               sk = zero
               do k=k1,ke
                  sk = sk - fguv(nh+n1,k,it)* fguv(n2,k,it)
               enddo
               do k=k1a,kea
                  sk = sk - fguv(nh+n1,k,it+2)*fguv(n2,k,it+2)
               enddo
               aka(il,it) = -2*i12*sk
               if(n1 .eq. n2) spk(it) = spk(it)+aka(il,it)
            enddo
         enddo
      enddo   ! ib
      spk(it)=half*spk(it)
      if (lpr) then
      write(l6,*) '****** END DENSSH **********************************'
      endif
c
      return
C-end-DENSSH
      end

c======================================================================c
      subroutine dinout(is,lpr)
c======================================================================c
c
c     IS = 1 : for ININK = 0  reads pairing field Delta from tape
c                  ININK = 1  calculates pairing field Delta
c     IS = 2 : writes pairing field Delta to tape
c
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /initia/ vin,rin,ain,inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      data laka/12/
c
      if(inink.eq.1 .and. is.eq.1) return
      if (lpr)
     &write(l6,*) ' ****** BEGIN DINOUT ****************************'
c
c
c==== initialize or read pairing potential
      if (is.eq.1) then
         do ib = 1,nb2x
            do i = 1,nhhx
	           de(i,ib) = zero
	        enddo   ! i
         enddo   ! ib
c
c---- reads pairing potential
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1)
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           do n2 = 1,nf
 	              read(laka) (de(n1+(n2-1)*nh,m),n1=n2,nf)
	              do n1 = n2,nf
                     de(n2+(n1-1)*nh,m) = de(n1+(n2-1)*nh,m)
                  enddo   ! n1
               enddo   ! n2
            enddo   ! ib
         enddo   ! it
         close(laka)
c
      if (lpr) write(l6,*) ' Pairing field has been read'
c
c==== writing of the pairing potential
      elseif (is.eq.2) then
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1)
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           do n2 =  1,nf
 	              write(laka) (de(n1+(n2-1)*nh,m),n1=n2,nf)
               enddo   ! n2
            enddo   ! ib
         enddo   ! it
         close(laka)
      else
         stop 'in DINOUT: is wrong'
      endif   ! is
c
      if (lpr) then
         write(l6,*) ' ****** END DINOUT ****************************'
      endif
      return
c-end-DINOUT
      end

c======================================================================c

      subroutine dirhb(it,lpr)

c======================================================================c
c
c     solves the Dirac-Equation in spherical oscillator basis
c     IT  = 1 neutrons
c           2 protons
c
c     units:    fields and Hamiltonian in MeV
c               eigenvalues in MeV
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'

      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
      character tk*8                                            ! blolev
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
      character bbb*1
c
      dimension hb(nhfbqx),e(nhfbx),ez(nhfbx)

      common /hamham/ hh(nhhx,nb2x)
      common /deldel/ de(nhhx,nb2x)
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /fermi / ala(2),tz(2)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,npr(2),nucnam
      common /physco/ hqc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      data epsl/1.d-8/,maxl/150/

c
c
      if (lpr) then
      write(l6,*) '****** BEGIN DIRHFB ***************************'
      endif
c
c
c
c======================================================================c
c     with pairing: Dirac-Hartree-Bogoliubov-Equation
c======================================================================c
c
c---- loop over the different k-blocks


      dl    = 100.d0
      xh    = ala(it) + dl
      xl    = ala(it) - dl
      al    = ala(it)

      do lit=1,maxl
         snold=sn
         sn  = zero
         alx = al
         klp = 0
         kla = 0
         do ib = 1,nb
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            nhb = nh+nh
            m   = ib + (it-1)*nbx

            do i2=1,nh
               do i1=1,nh
                  hla = hh(i1+(i2-1)*nh,m)
                  dla = de(i1+(i2-1)*nh,m)
                  hb(i1+(i2-1)*nhb)       =  hla
                  hb(nh+i1+(nh+i2-1)*nhb) = -hla
                  hb(nh+i1+(i2-1)*nhb)    =  dla
                  hb(nh+i2+(i1-1)*nhb)    =  dla
               enddo
               hb(i2+(i2-1)*nhb)       = hb(i2+(i2-1)*nhb) - al
               hb(nh+i2+(nh+i2-1)*nhb) = hb(nh+i2+(nh+i2-1)*nhb)+al
            enddo

c
c------- Diagonalisation of the Dirac equation:
            call sdiag(nhb,nhb,hb,e,hb,ez,+1)
c
c------- store eigenvalues and wave functions
c------- particles
            ka(ib,it) = klp
            do k = 1,nf
	       klp = klp + 1
	       ibk(klp,it) = ib
	       write(tk(klp,it),100) k,tb(ib)
 100        format(i2,a5)
               equ(klp,it) = e(nh+k)
               do i = 1,nhb
                  fguv(i,klp,it) = hb(i+(nh+k-1)*nhb)
               enddo
               v2 = zero
               do i = 1,nh
                  v2 = v2 + fguv(nh+i,klp,it)**2
               enddo
               sn=sn+2*v2
            enddo   ! nf
 	    kd(ib,it) = klp - ka(ib,it)
c
c------- antiparticles
            ka(ib,it+2) = kla
            do k = 1,ng
               kla = kla + 1
               ibk(kla,it+2) = ib
               write(tk(kla,it+2),100) ng+1-k,tb(ib)
               equ(kla,it+2) = e(ng-k+1)
               do i = 1,nhb
                  fguv(i,kla,it+2) = hb(i+(ng-k)*nhb)
               enddo
               v2 = zero
               do i = 1,nh
                  v2 = v2 + fguv(nh+i,kla,it+2)**2
               enddo
               sn = sn+2*v2
            enddo   ! k
            kd(ib,it+2) = kla - ka(ib,it+2)
c
         enddo   ! ib
         nk(it)   = klp
         nk(it+2) = kla

         dn = sn-tz(it)
         if (lit.gt.1) dd = (sn - snold)/(al - alold)
         alold = al
         if (abs(dn).lt.epsl) goto 30
         if (dn.lt.zero) then
            xl = al
         else
            xh = al
         endif
         if (lit.eq.1) then
            if(dabs(dn).le.0.1d0) then
               al=al-dn
            else
               al = al - 0.1d0*sign(one,dn)
            endif
         else
c
c           secant method
            if (dabs(dd).lt.1.d-20) dd = 1.d-20
            al    = al - dn/dd
            if (al.lt.xl.or.al.gt.xh) then
c
c              bisection
               al = half*(xl+xh)
               bbb = 'B'
            endif
         endif
c         write(6,*) sn
         if (abs(al-alold).lt.epsl) goto 30
c
c         if (lpr .and. lit.gt.10) then
c            write(l6,102) lit,'. L-Iteration: ',bbb,alold,dn,al
c            write(6,102)  lit,'. L-Iteration: ',bbb,alold,dn,al
c           read(*,*)
  102       format(i4,a,a1,3f13.8)
            bbb = ' '
c         endif
c
c---- end of lambda-loop
      enddo
      write(l6,'(a,i4,a)')
     &     ' Lambda-Iteration interupted after',lit-1,' steps'
      stop
   30 if (lpr) then
         write(l6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
         write(6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
  101    format(i4,a,i4,3f13.8)
      endif
      ala(it) =      al

c
c
      if (lpr) then
      write(l6,*) '****** END DIRHFB ******************************'
      endif
c
      return
C-end-DIRHFB
      end

c======================================================================c

      subroutine aprint(is,it,ns,ma,n1,n2,a,t1,t2,text)

c======================================================================c
C
C     IS = 1    Full matrix
C          2    Lower diagonal matrix
c          3    specially stored symmetric matrix
C
C     IT = 1    numbers for rows and columns
C          2    text for rows and numbers for columns
C          3    text for rows and columns
C
C     NS = 1     FORMAT  8F8.4      80 Coulums
C     NS = 2     FORMAT  8f8.2      80 Coulums
C     NS = 3     FORMAT 17F4.1      80 Coulums
C     NS = 4     FORMAT 30F4.1     120 Coulums
C     NS = 5     FORMAT  5F12.8     80 Coulums
C     NS = 6     FORMAT  5F12.4     80 Coulums
C     NS = 7     FORMAT  4E13.6     80 Coulums
C     NS = 8     FORMAT  8E15.8    130 Coulums
c
c----------------------------------------------------------------------c
      implicit double precision (a-h,o-z)
C
      character*8 t1(n1),t2(n2)
      character text*(*)
C
      dimension a(ma*n2)
C
      character*30 fmt1,fmt2
      character*20 fti,ftt,fmt(8),fmti(8),fmtt(8)
      dimension nsp(8)
c
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      data nsp/8,8,17,30,5,5,4,8/
      data fmt /'8f8.4)',            '8F8.2)',
     &          '17f4.1)',           '30f4.1)',
     &          '5f12.8)',           '5f12.4)',
     &          '4e13.6)',           '8e15.8)'/
      data fmti/'(11x,8(i4,4x))',    '(11x,8(i4,4x))',
     &          '(11x,17(1x,i2,1x))','(11x,30(1x,i2,1x))',
     &          '(11x,6(i4,8x))',    '(11x,10(i4,8x))',
     &          '(11x,5(i4,9x))',    '(11x,8(i4,11x))'/
      data fmtt/'(11x,8a8)',         '(11x,8a8)',
     &          '(11x,17a4)',        '(11x,30a4)',
     &          '(11x,5(a8,2x))',    '(11x,5(a8,4x))',
     &          '(11x,4(a8,5x))',    '(11x,8(a8,7x))'/
C
      fmt1   = '(4x,i3,4x,' // fmt(ns)
      fmt2   = '(1x,a8,2x' // fmt(ns)
      fti    = fmti(ns)
      ftt    = fmtt(ns)
      nspalt = nsp(ns)

C
      if (ma.eq.0.or.n1.eq.0.or.n2.eq.0) return
      write(l6,'(//,3x,a)') text
C
      ka = 1
      ke = nspalt
      nteil = n2/nspalt
      if (nteil*nspalt.ne.n2) nteil = nteil + 1
C
      do  10  nt = 1,nteil
      if (n2.gt.nspalt)  write(L6,100)  nt
  100 format(//, 10x,'Part',i5,' of the Matrix',/)
      if(nt.eq.nteil) ke = n2
      if (it.lt.3) then
        write(L6,fti) (k,k=ka,ke)
      else
        write(L6,ftt) (t2(k),k=ka,ke)
      endif
C
      do 20  i=1,n1
         kee=ke
         if (is.ge.2.and.ke.gt.i) kee=i
         if (ka.gt.kee) goto 20
         if (is.eq.3) then
            if (it.eq.1) then
               write(l6,fmt1) i,(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            else
               write(l6,fmt2) t1(i),(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            endif
         else
            if (it.eq.1) then
               write(l6,fmt1) i,(a(i+(k-1)*ma),k=ka,kee)
            else
               write(l6,fmt2) t1(i),(a(i+(k-1)*ma),k=ka,kee)
            endif
         endif
   20 continue
c
      ka=ka+nspalt
      ke=ke+nspalt
   10 continue
C
      return
C-end-APRINT
      end
C=======================================================================

      subroutine gfv

C=======================================================================
C
C     Calculates sign, sqrt, factorials, etc. of integers and half int.
C
c     iv(n)  =  (-1)**n
c     sq(n)  =  sqrt(n)
c     sqi(n) =  1/sqrt(n)
c     sqh(n) =  sqrt(n+1/2)
c     shi(n) =  1/sqrt(n+1/2)
c     fak(n) =  n!
c     fad(n) =  (2*n+1)!!
c     fdi(n) =  1/(2*n+1)!!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
c     wfd(n) =  sqrt((2*n+1)!!)
c     gm2(n) =  gamma(n+1/2)
c     gmi(n) =  1/gamma(n+1/2)
c     wg(n)  =  sqrt(gamma(n+1/2))
c     wgi(n) =  1/sqrt(gamma(n+1/2))
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
c     include 'dirhb.par'
      parameter (igfv = 100)
c
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /gfvfad/ fad(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvfdi/ fdi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwfd/ wfd(0:igfv)
      common /gfvgm2/ gm2(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
c
      zero  = 0.d0
      one   = 1.d0
      two   = 2.d0
      half  = one/two
      third = one/3.d0
      pi    = 4*atan(one)
c
      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  one
      fad(0) =  one
      fi(0)  =  one
      fdi(0) =  one
      wf(0)  =  one
      wfi(0) =  one
      wfd(0)=  one
c     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
      do i = 1,igfv
         iv(i)  = -iv(i-1)
         sq(i)  = dsqrt(dfloat(i))
         sqi(i) = one/sq(i)
         sqh(i) = sqrt(i+half)
         shi(i) = one/sqh(i)
         fak(i) = i*fak(i-1)
	 fad(i) = (2*i+1)*fad(i-1)
         fi(i)  = one/fak(i)
         fdi(i) = one/fad(i)
         wf(i)  = sq(i)*wf(i-1)
         wfi(i) = one/wf(i)
	 wfd(i) = sqrt(fad(i))
         gm2(i) = (i-half)*gm2(i-1)
         gmi(i) = one/gm2(i)
         wg(i)  = sqh(i-1)*wg(i-1)
         wgi(i) = one/wg(i)
      enddo
c
c     write(6,*) ' ****** END GFV *************************************'
      return
c-end-GFV
      end
c======================================================================c

      function itestc()

c======================================================================c
c
C    yields 1, if interrupted
c           2, if convergence
c
c    the iteration is determined by the parameter XMIX
c    it can be fixed, or automatically adjusted according to
c
c     IAUT = 0 fixed value for XMIX
c          = 1 automatic adjustment of XMIX
c
c     INXT = 0 with out inputs from the console
c          > 0 after INXT iterations INX and XMIX are read
c
c     INX  = 0 immediate interruption of the iteration
c          > 0 further INX steps with fixed XMIX, which is read
c          < 0 further ABS(INX) steps with automatic change of XMIX
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      integer itestc
c
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /broyde2/ ibroyd
c
      if (ii.ge.2.and.si.lt.epsi) then
         itestc = 2
         return
      endif
c
c     if you want to change xmix at the console, please remove the
c     next line !
      if (ibroyd.eq.0) then
         inxt = maxi
c     change of XMIX by reading from the console:
         if (inxt.eq.ii) then
            write(6,*)
     &     'next stop? (0 right now, >0 fixed xmix, <0 autom. xmix)'
            read(*,*)  inx
c        for running the program without console
            inx=-30

            if (inx.eq.0) then
               itestc = 1
               return
            endif
            if (inx.lt.0) then
               iaut = 1
            endif
            if (inx.gt.0) then
               iaut = 0
            endif
            inxt = ii+iabs(inx)
            write(6,*) 'new value for xmix?'
            read(*,*) xmix
            write(6,*) inxt,xmix
            xmix0 = xmix
         endif
c
c     automatic change of XMIX:
         if ((si.lt.siold).and.iaut.eq.1) THEN
            xmix = xmix * 1.0
            if (xmix.gt.xmax) xmix = xmax
         else
            xmix = xmix0
         endif
      endif
      siold  = si
      itestc = 0
c
      return
c-end-ITESTC
      end
C=======================================================================

      subroutine lingd(ma,mx,n,m,a,x,d,ifl)

C=======================================================================
C
C     solves the system of linear equations A*X = B
C     at the beginning the matrix B is stored in X
C     during the calculation it will be overwritten
C     D is the determinant of A
C
C-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
C
c     solves the system A*X=B, where B is at the beginning on X
c     it will be overwritten lateron, d is the determinant
C
      dimension  a(ma,m),x(mx,m)
C
      data tollim/1.d-10/,one/1.d0/,zero/0.d0/
C
      ifl=1
      p=zero
      do 10 i=1,n
         q=zero
         do 20 j=1,n
 20      q=q+ abs(a(i,j))
         if (q.gt.p)   p=q
 10   continue
      tol=tollim*p
      d=one
      do 30 k=1,n
         p=zero
         do 40 j=k,n
            q = abs(a(j,k))
            if (q.lt.p) goto 40
            p=q
            i=j
 40      continue
         if (p.gt.tol) goto 70
         write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
  200    format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e10.4,
     1 ' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1)
         ifl=-1
         return
   70    cp=one/a(i,k)
         if (i.eq.k) goto 90
         d=-d
         do 81 l=1,m
            cq=x(i,l)
            x(i,l)=x(k,l)
   81       x(k,l)=cq
         do 80 l=k,n
            cq=a(i,l)
            a(i,l)=a(k,l)
   80       a(k,l)=cq
   90       d=d*a(k,k)
            if (k.eq.n) goto 1
            k1=k+1
            do 120 i=k1,n
               cq=a(i,k)*cp
               do 106 l=1,m
  106             x(i,l)=x(i,l)-cq*x(k,l)
               do 120 l=k1,n
  120             a(i,l)=a(i,l)-cq*a(k,l)
   30 continue
    1 do 126 l=1,m
  126    x(n,l)=x(n,l)*cp
         if (n.eq.1) return
         n1=n-1
         do 140 k=1,n1
            cp=one/a(n-k,n-k)
            do 140 l=1,m
               cq=x(n-k,l)
               do 141 i=1,k
  141             cq=cq-a(n-k,n+1-i)*x(n+1-i,l)
  140          x(n-k,l)=cq*cp
c
c
      return
c-end-LINGD
      end
C=======================================================================

      subroutine nucleus(is,npro,te)

C=======================================================================
C
C     is = 1 determines the symbol for a given proton number npro
c          2 determines the proton number for a given symbol te
c
C-----------------------------------------------------------------------
C
      PARAMETER (MAXZ=140)
C
      CHARACTER TE*2,T*(2*MAXZ+2)
C
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMr10111213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
c
c ... Rf is called also as Ku (kurchatovium)
c ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
c ... Ha is called also as Db (Dubnium)
c
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong'
         te = t(2*npro+1:2*npro+2)
         return
      else
c
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
	       if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
c
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
c
      stop
C-END-NUCLEUS
      END
c
c=======================================================================

      subroutine ord(n,e)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
            endif
         endif
   10 continue
c
      return
c-end-ORD
      end

c=======================================================================

      subroutine ordi(n,e,mu)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),mu(n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               mk    = mu(k)
               mu(k) = mu(i)
               mu(i) = mk
            endif
         endif
   10 continue
c
      return
c-end-ORDI
      end
c=======================================================================

      subroutine ordii(n,e,ix,mx)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),ix(n),mx(n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               mk    = mx(k)
               mx(k) = mx(i)
               mx(i) = mk
               ik    = ix(k)
               ix(k) = ix(i)
               ix(i) = ik
            endif
         endif
   10 continue
c
      return
c-end-ORDI
      end
c=======================================================================

      subroutine ordx(n,e,a1,a2,bb)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),a1(n),a2(n),bb(n,n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               x     = a1(k)
               a1(k) = a1(i)
               a1(i) = x
               x     = a2(k)
               a2(k) = a2(i)
               a2(i) = x
               do j = 1,n
                  x       = bb(j,k)
                  bb(j,k) = bb(j,i)
                  bb(j,i) = x
               enddo
            endif
         endif
   10 continue
c
      return
c-end-ORD3
      end
c=======================================================================

      subroutine ordx2(n,e,a1,a2,a3,bb)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),a1(n),a2(n),a3(n),bb(n,n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               x     = a1(k)
               a1(k) = a1(i)
               a1(i) = x
               x     = a2(k)
               a2(k) = a2(i)
               a2(i) = x
               x     = a3(k)
               a3(k) = a3(i)
               a3(i) = x
               do j = 1,n
                  x       = bb(j,k)
                  bb(j,k) = bb(j,i)
                  bb(j,i) = x
               enddo
            endif
         endif
   10 continue
c
      return
c-end-ORD3
      end
c======================================================================c

      subroutine osc(n,l,x,rnl)

c======================================================================c
c
c     calculates the radial functions for the spherical oscillator
c
c     the wave function R_nl(r) of the spherical oscillator are:
c
c     phi(r,Omega) = b^(-3/2) * R_nl(r) * Y_ljm(Omega)
c
c     R_nl(r) = N_nl * r**l * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
c
c     N_nl    = sqrt(2 * (n-1)!/(n+l-1/2)!)     and    x=r/b
c
c     R_nl is normalized in such way that the norm integral reads
c
c     \int dr r**2 R_nl(r)^2 = 1
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c     include 'dirhb.par'
      parameter (igfv = 100)
c
      dimension rnl(n)
c
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      xx = x*x
      if (l.eq.0) then
	 xl = one
      else
	 xl = x**l
      endif
      rnl(1) = sq(2)*wgi(l+1)*exp(-half*xx)*xl
      rnl(2) = rnl(1)*(l+1.5d0-xx)*shi(l+1)
      do i = 3,n
         rnl(i) = ((2*i+l-2.5d0-xx)*rnl(i-1) -
     &             sq(i-2)*sqh(i-2+l)*rnl(i-2))*sqi(i-1)*shi(i-1+l)
      enddo
c
      return
c-end-OSC
      end
C=======================================================================

      subroutine sdiag(nmax,n,a,d,x,e,is)

C=======================================================================
C
C     A   matrix to be diagonalized
C     D   eigenvalues
C     X   eigenvectors
C     E   auxiliary field
C     IS = 1  eigenvalues are ordered and major component of X is positiv
C          0  eigenvalues are not ordered
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)
C
      data tol,eps/1.e-32,1.e-10/
C
      if (n.eq.1) then
         d(1)=a(1,1)
         x(1,1)=1.
         return
      endif
c
      do 10 i=1,n
      do 10 j=1,i
   10    x(i,j)=a(i,j)
c
ccc   householder-reduktion
      i=n
   15 if (i-2) 200,20,20
   20 l=i-2
      f=x(i,i-1)
      g=f
      h=0
      if (l) 31,31,32
   32 do 30 k=1,l
   30 h=h+x(i,k)*x(i,k)
   31 s=h+f*f
      if (s-tol) 33,34,34
   33 h=0
      goto 100
   34 if (h) 100,100,40
   40 l=l+1
      g= dsqrt(s)
      if (f.ge.0.) g=-g
      h=s-f*g
      hi=1.d0/h
      x(i,i-1)=f-g
      f=0.0
      if (l) 51,51,52
   52 do 50 j=1,l
      x(j,i)=x(i,j)*hi
      s=0.0
      do 55 k=1,j
   55 s=s+x(j,k)*x(i,k)
      j1=j+1
      if (l-j1) 57,58,58
   58 do 59 k=j1,l
   59 s=s+x(k,j)*x(i,k)
   57 e(j)=s*hi
   50 f=f+s*x(j,i)
   51 f=f*hi*.5d0
c
      if (l) 100,100,62
   62 do 60 j=1,l
      s=x(i,j)
      e(j)=e(j)-f*s
      p=e(j)
      do 65 k=1,j
   65 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
   60 continue
  100 continue
      d(i)=h
      e(i-1)=g
      i=i-1
      goto 15
c
ccc   Bereitstellen der Transformationmatrix
  200 d(1)=0.0
      e(n)=0.0
      b=0.0
      f=0.0
      do 210 i=1,n
      l=i-1
      if (d(i).eq.0.) goto 221
      if (l) 221,221,222
  222 do 220 j=1,l
      s=0.0
      do 225 k=1,l
  225 s=s+x(i,k)*x(k,j)
      do 226 k=1,l
  226 x(k,j)=x(k,j)-s*x(k,i)
  220 continue
  221 d(i)=x(i,i)
      x(i,i)=1
      if (l) 210,210,232
  232 do 230 j=1,l
      x(i,j)=0.0
  230 x(j,i)=0.0
  210 continue
c
ccc   Diagonalisieren der Tri-Diagonal-Matrix
      DO 300 L=1,N
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h
c
ccc   Test fuer Splitting
      do 310 j=l,n
      if ( abs(e(j)).le.b) goto 320
  310 continue
c
ccc   test fuer konvergenz
  320 if (j.eq.l) goto 300
  340 p=(d(l+1)-d(l))/(2*e(l))
      r= dsqrt(p*p+1.d0)
      pr=p+r
      if (p.lt.0.) pr=p-r
      h=d(l)-e(l)/pr
      do 350 i=l,n
  350 d(i)=d(i)-h
      f=f+h
c
ccc   QR-transformation
      p=d(j)
      c=1.d0
      s=0.0
      i=j
  360 i=i-1
      if (i.lt.l) goto 362
      g=c*e(i)
      h=c*p
      if ( abs(p)- abs(e(i))) 363,364,364
  364 c=e(i)/p
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.d0/r
      goto 365
  363 c=p/e(i)
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r
      s=1.d0/r
      c=c/r
  365 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 368 k=1,n
         h=x(k,i+1)
         x(k,i+1)=x(k,i)*s+h*c
  368    x(k,i)=x(k,i)*c-h*s
      goto 360
  362 e(l)=s*p
      d(l)=c*p
      if ( abs(e(l)).gt.b) goto 340
c
ccc   konvergenz
  300 d(l)=d(l)+f
c
      if (is.eq.0) return
ccc   ordnen der eigenwerte
      do 400 i=1,n
      k=i
      p=d(i)
      j1=i+1
      if (j1-n) 401,401,400
  401 do 410 j=j1,n
      if (d(j).ge.p) goto 410
      k=j
      p=d(j)
  410 continue
  420 if (k.eq.i) goto 400
      d(k)=d(i)
      d(i)=p
      do 425 j=1,n
      p=x(j,i)
      x(j,i)=x(j,k)
  425 x(j,k)=p
  400 continue
c
c     signum
      do k = 1,n
         s = 0.0d0
         do i = 1,n
            h = abs(x(i,k))
            if (h.gt.s) then
               s  = h
               im = i
            endif
         enddo   ! i
         if (x(im,k).lt.0.0d0) then
            do i = 1,n
               x(i,k) = - x(i,k)
       	    enddo
         endif
      enddo   ! k
c
      return
c-end-SDIAG
      end
c=======================================================================
c
      SUBROUTINE gauher(x,w,n)
c     Calculates gaussian mesh points
c
c=======================================================================
c
      INTEGER n,MAXIT
      DOUBLE PRECISION w(n),x(n)
      DOUBLE PRECISION EPS,PIM4
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(dble(2*n+1))-1.85575d0*(2*n+1)**(-.16667d0)
        else if(i.eq.2)then
          z=z-1.14d0*n**.426d0/z
        else if (i.eq.3)then
          z=1.86d0*z-.86d0*x(1)
        else if (i.eq.4)then
          z=1.91d0*z-.91d0*x(2)
        else
          z=2.d0*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=sqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        stop 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
      return
      END
c=======================================================================
c
      function spur(n,a)
c
c=======================================================================
c
c     calculates the trace of the matrix A
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension a(n*n)
c
      s = 0.d0
      do i = 1,n
         s = s + a(i+(i-1)*n)
      enddo
      spur = s
c
      return
c-end-SPUR
      end
c=====================================================================c

      double precision function trabt(ma,mb,n1,n2,aa,bb)

c======================================================================c
c
c     calculaties the trace( A * BT )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      dimension aa(ma,n2),bb(mb,n2)
c
      s = 0.d0
      do i = 1,n1
      do k = 1,n2
         s = s + aa(i,k)*bb(i,k)
      enddo
      enddo
      trabt = s
c
      return
c-end-TRABT
      end
c=====================================================================c

c======================================================================c
c
      subroutine expect(is,lpr)
c
c======================================================================c
c
c     calculates expectation values
c
c     is = 1    short version, during the iteration
c     is = 2    long version, at the end
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character nucnam*2                                 ! common nucnuc
c
      dimension xn(3),xs(3),r2(3),q0(3),q2(3),hh(3)
      dimension bet2(3),gam(3),het4(3)
      dimension ekt(3),epart(3),ept(3)
      dimension emes(4)
      dimension rmom(9,3),a0(3),a2(3)
      dimension npa(3)
c
      common /baspar/ hom,hb0,b0
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /erwar / ea,rms,betg,gamg
      common /cstr4/ calcxx,calcyy,calczz
      common /cstr5/ fac0,fac2
      common /fermi / ala(2),tz(2)
      common /gaucor/ wdcor(ngauss)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /pair  / del(2),spk(2),spk0(2)
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
c
      if (lpr) then
      write(l6,*) '****** BEGIN EXPECT ********************************'
      endif
c
c
c
c======================================================================c
c---- particle number, radii, deformations
c======================================================================c
      npa(3)=0
      do it=1,2
         npa(it)=npr(it)
         npa(3)=npa(3)+npa(it)
      enddo
c
      do it = 1,3
         do i=1,9
            rmom(i,it)=zero
         enddo
      enddo
c
      do iz = 0,ngh
         z  = zb(iz)
         zz = z**2
      do iy = 0,ngh
         y  = yb(iy)
         yy = y**2
      do ix = 0,ngh
         x  = xb(ix)
         xx = x**2
         i  = 1+ix + iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
c
c------- root mean square radius
         rrp= xx + yy
         rr = rrp + zz
c
c------- quadrupole moment
         rq = 3*zz - rr
c
c------- hexadecupole moment
         rh = 8*zz**2 - 24*zz*rrp + 3*rrp**2

c
         do it = 1,itx
            rv     = ro(i,2) + iv(it)*ro(i,4)
	    rv     = half*rv*wdcor(i)
            rs     = ro(i,1) + iv(it)*ro(i,3)
	    rs     = half*rs*wdcor(i)
            rmom(1,it) = rmom(1,it) + rv
            rmom(2,it) = rmom(2,it) + rs
            rmom(3,it) = rmom(3,it) + rv*rr
            rmom(4,it) = rmom(4,it) + rv*rq
            rmom(5,it) = rmom(5,it) + rv*(xx-yy)
            rmom(6,it) = rmom(6,it) + rv*rh
            rmom(7,it) = rmom(7,it) + rv*xx
            rmom(8,it) = rmom(8,it) + rv*yy
            rmom(9,it) = rmom(9,it) + rv*zz
         enddo   ! it
      enddo   ! ix
      enddo   ! iy
      enddo   ! iz

      if(itx.eq.1) then
        do i=1,9
           rmom(i,2) = rmom(i,1)
        enddo
      endif

      rmom(1,3) = rmom(1,1)+rmom(1,2)
      rmom(2,3) = rmom(2,1)+rmom(2,2)
      rmom(3,3) = rmom(3,1)+rmom(3,2)
      rmom(4,3) = rmom(4,1)+rmom(4,2)
      rmom(5,3) = rmom(5,1)+rmom(5,2)
      rmom(6,3) = rmom(6,1)+rmom(6,2)
      rmom(7,3) = rmom(7,1)+rmom(7,2)
      rmom(8,3) = rmom(8,1)+rmom(8,2)
      rmom(9,3) = rmom(9,1)+rmom(9,2)

      fac20 = sqrt(5.d0/16/pi)
      fac22 = sqrt(15.d0/32/pi)
      fac40 = sqrt(9.d0/256/pi)
      f4    = sqrt(pi/16)
      r0    = 1.2*amas**third
      do it=1,itx
         xn(it) = rmom(1,it)
         xs(it) = rmom(2,it)
         r2(it) = rmom(3,it)
         r2(it) = sqrt(r2(it)/xn(it))
         q0(it) = fac20*rmom(4,it)
         q2(it) = fac22*rmom(5,it)
         hh(it) = fac40*rmom(6,it)
      enddo
c
      xn(3) = xn(1) + xn(2)
      xs(3) = xs(1) + xs(2)
      r2(3) = sqrt((npr(1)*r2(1)**2+npr(2)*r2(2)**2)/amas)
      rc    = sqrt(r2(2)**2 + 0.64)
      rms   = r2(3)

      q0(3) = q0(1) + q0(2)
      q2(3) = q2(1) + q2(2)
      hh(3) = hh(1) + hh(2)

      do it=1,3
         a0(it) = q0(it)*4*pi/(3*xn(it)*r0**2)
         a2(it) = q2(it)*4*pi/(3*xn(it)*r0**2)
         if(abs(a0(it)) .gt. 1.d-8) then
            gam(it) = atan(sq(2)*a2(it)/a0(it))
         else
            gam(it) = 0.d0
         endif
         bet2(it) = a0(it)/cos(gam(it))
         het4(it) = f4*rmom(6,it)/(xn(it)*r0**4)
         gam(it) = gam(it)*180.d0/pi
      enddo

      betg  = bet2(3)
      gamg  = gam(3)

      calcxx = rmom(7,3)
      calcyy = rmom(8,3)
      calczz = rmom(9,3)

      ecstr=fac0*rmom(4,3)+fac2*rmom(5,3)

c
c======================================================================c
c---- single particle energies, kinetic energies and pairing energies
c======================================================================c

      do it = 1,itx
c------- kinetic energy
         ekt(it) = ekin(it)
c
c------- particle energy
         epart(it) = epar(it)
c
c------- pairing energy
         ept(it) = epair(it)
c
      enddo   ! it
      if (itx.eq.1) then
	     ekt(2)   = ekt(1)
	     epart(2) = epart(1)
	     ept(2)   = ept(1)
      endif
      ekt(3)   = ekt(1) + ekt(2)
      epart(3) = epart(1) + epart(2)
      ept(3)   = ept(1) + ept(2)
c======================================================================c
c---- field energies
c======================================================================c
c
      call efield(emes,er,ecou)
      escsc  = emes(1)
      escve  = emes(2)
      evesc  = emes(3)
      eveve  = emes(4)
c
c======================================================================c
c---- Total energy
c======================================================================c
      etot0 = ekt(3) + escsc + escve + eveve
     &               + ecou  + ept(3)
      etot1 = epart(3) -ecstr
     &      - escsc - escve - eveve - ecou - er + ept(3)
      etot  = etot0 + ecm(3)
      etest = etot1-etot0
      ea    = etot/amas
c
c
c======================================================================c
c---- printout
c======================================================================c
      if (.not.lpr) return
      write(l6,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
c
c     particle number
      write(l6,200) ' particle number .....',xn
      write(l6,200) ' trace scalar density ',xs
c
c     Lambda
      write(l6,200) ' lambda ..............',ala
c     trace of kappa
      write(l6,200) ' spk..................',spk
c
c     rms-Radius
      write(l6,200) ' rms-Radius ..........',r2
c
c     charge-Radius
      write(l6,201) ' charge-Radius, R0....',rc
c
c     quadrupole-moment
      write(l6,200) ' quadrupole moment ...',(rmom(4,i),i=1,3)
      write(l6,200) '    <Q20> ............',(q0(i),i=1,3)
      write(l6,200) ' beta ................',(bet2(i),i=1,3)
      write(l6,200) ' gamma................',(gam(i),i=1,3)
      write(l6,200) ' a20 .................',(a0(i),i=1,3)
      write(l6,200) ' a22 .................',(a2(i),i=1,3)

c     hexadecupole moment
      write(l6,200) ' hexadecupole moment .',(rmom(6,i),i=1,3)
      write(l6,200) '    <Q40> ............',(hh(i),i=1,3)
      write(l6,200) ' heta ................',(het4(i),i=1,3)
      write(l6,*) ' '
c
      write(l6,*) ' '
c
c     single-particle energy
      write(l6,200) ' Particle Energy .....',epart
      write(l6,202) ' Selfconsistency Test.',etest
      write(l6,*) ' '
c
c     kinetic energy
      write(l6,200) ' Kinetic Energy ......',ekt
c
c     sigma energy
      write(l6,202) ' E-scsc  .............',emes(1)
c
c     omega energy
      write(l6,202) ' E-scve ..............',emes(2)
c
c     delta energy
      write(l6,202) ' E-vesc ..............',emes(3)
c
c     rho energy
      write(l6,202) ' E-veve ..............',emes(4)
c
c     rearrangement contribution
      write(l6,202) ' Rearrangement contr..',er
c----
c
c     Coulomb energy
      write(l6,202) ' Coulomb direct ......',ecou
c
c     Constrained energy
      write(l6,202) ' Constrained energy ..',ecstr
c
c     pairing-energy
      write(l6,200) ' Pairing Energy .....',ept

c     total energy without center of mass correction
      write(l6,202)    ' Sum without E-cm ....',etot0
c
c     center of mass correction
      if (icm.eq.0) then
         write(l6,202) ' E-cm  3/4*hom .......',ecm(3)
      elseif (icm.eq.1) then
         write(l6,202) ' E-cm  3/4*hom .......',ecm(3)
      elseif (icm.eq.2) then
         write(l6,200) ' DE-cm <P**2>/2M .....',ecmd
         write(l6,200) ' NE-cm <P**2>/2M .....',ecmn
         write(l6,200) ' E-cm <P**2>/2M ......',ecm
      else
         stop 'in ERWAR: icm not properly defined'
      endif

c
c     total energy
      write(l6,202) ' Total Energy ........',etot
c
c     energy per particle
      ea = etot/amas
      write(l6,202) ' E/A .................',ea
c
c
  200 format(a,3f15.6)
  201 format(a,15x,2f15.6)
  202 format(a,30x,2f15.6)
c
c---------------------------
  100 format(a,7x,3f15.6)
  101 format(a,37x,f15.6)
  102 format(a,35x,f15.6)
c
c
      if (lpr) then
      write(l6,*) '****** END EXPECT **********************************'
      endif
c
      return
c-end-EXPECT
      end
c======================================================================c
c
      real*8 function ekin(it)
c
c======================================================================c
c
c     calculates the kinetic energy
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tb*5                                            ! blokap
c
      dimension  h0(nhhx)
      common /baspar/ hom,hb0,b0
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /physco/ hbc,alphi,r0
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /single/ sp(nfgx,nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
c
c
      ek = zero
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
	     nhfb = nh + nh
         m   = ib + (it-1)*nbx
c
c------- construction of the free Dirac-operator H0
         emcc2 = 2*amu*hbc
         f = hbc/b0
         do n2 = 1,nf
            do n1 = 1,nf
               h0(n1+(n2-1)*nh) = zero
            enddo
            do n1 = 1,ng
               h0(nf+n1+(n2-1)*nh) = f*sp(n1+(n2-1)*ng,ib)
               h0(n2+(nf+n1-1)*nh) = h0(nf+n1+(n2-1)*nh)
            enddo
         enddo
         do n2 = nf+1,nh
            do n1 = nf+1,nh
               h0(n1+(n2-1)*nh) = zero
            enddo
            h0(n2+(n2-1)*nh) = - emcc2
         enddo
c------- calculation of Tr(H0 * RO)
         ek = ek + trabt(nh,nh,nh,nh,h0,rosh(1,m))
c
      enddo   ! ib
      ekin = ek
c
      return
c-end-EKIN
      end
c======================================================================c
c
      real*8 function epar(it)
c
c======================================================================c
c
c     calculates the particle energy
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tk*8                                            ! blolev
      character tb*5                                            ! blokap
c
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /hamham/ hh(nhhx,nb2x)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
c
      ep = zero
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
	     nhfb = nh + nh
         m   = ib + (it-1)*nbx
         ep = ep + trabt(nh,nh,nh,nh,hh(1,m),rosh(1,m))
      enddo   ! ib
      epar = ep
c
      return
c-end-EPAR
      end
c======================================================================c
c
      real*8 function epair(it)
c
c======================================================================c
c
c     calculates the pairing energy
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      character tb*5                                            ! blokap
c
c
      common /mathco/ zero,one,two,half,third,pi
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /rhoshe/ rosh(NHHX,NB2X),aka(MVFX,2)
      common /deldel/ de(nhhx,nb2x)

      s = zero
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         m  = ib + (it-1)*NBX
         do n2 =  1,nf
         do n1 = n2,nf
            il = il + 1
            s  = s + de(n1+(n2-1)*nh,m)*aka(il,it)
         enddo
         enddo
      enddo   ! ib
      epair = -s/2

      return
c-end-EPAIR
      end

c======================================================================c
c
      subroutine efield(emes,er,ecou)
c
c======================================================================c
c
c     calculates  field energies
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      dimension emes(4)
c
      common /coulmb/   cou(ngauss),drvp(ngauss)
      common /coupld/ ddmes(4)
      common /couplf/ ff(ngauss,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /fields/ phi(ngauss,4)
      common /gaucor/   wdcor(ngauss)
      common /mathco/   zero,one,two,half,third,pi
      common /optopt/   itx,icm,icou,ipc,inl,idd
      common /physco/   hbc,alphi,r0
      common /tapes /   l6,lin,lou,lwin,lwou,lplo
c
      enl = zero
c
c
c======================================================================c
c---- field energies
c======================================================================c
c
c-----meson-fields
c
      if (ipc.eq.0) then
         er=zero
         do m = 1,4
            s = zero
            do i = 1,ngauss
               s = s + ggmes(m)*ff(i,m,1)*phi(i,m)*ro(i,m)*wdcor(i)
               er=er + ggmes(m)*ff(i,m,2)*phi(i,m)
     &                         *ro(i,m)*ro(i,2)*wdcor(i)
            enddo   ! i
            emes(m) = half*hbc*s
         enddo   ! m
         er=er*hbc
c
c-----point coupling
      elseif (ipc.eq.1) then
         er=zero
         do m = 1,4
            s = zero
            do i = 1,ngauss
               s = s + ggmes(m)*ff(i,m,1)*ro(i,m)**2*wdcor(i)
               er = er+ggmes(m)*ff(i,m,2)*ro(i,m)**2*ro(i,2)*wdcor(i)
            enddo   ! i
c
c           derivative terms
            do i = 1,ngauss
               s = s + ddmes(m)*ro(i,m)*dro(i,m)*wdcor(i)
            enddo   ! i
            emes(m) = half*hbc*s
         enddo   ! m
         er=half*er*hbc
c
      else
         stop 'in EFIELD: ipc not properly defined'
      endif   ! ipc
c
c
c
c======================================================================c
c---- Coulomb energy
c======================================================================c
      ecou  = zero
      if (icou.ne.0) then
         do i = 1,ngauss
            rv = half*(ro(i,2) + ro(i,4))
            ecou  = ecou + cou(i)*rv*wdcor(i)
         enddo   ! i
      endif   ! icou
      ecou  = hbc*ecou/2
      return
c-end-EFIELD
      end

c======================================================================c

      subroutine field(lpr)

c======================================================================c
c
c     calculation of the meson-fields in the oscillator basis
c     the fields are given in (fm^-1)
c
c     meson fields:  sig(i) = phi(i,1)*ggsig/gsig
c                    ome(i) = phi(i,2)*ggome/gome
c                    del(i) = phi(i,3)*ggdel/gdel
c                    rho(i) = phi(i,4)*ggrho/grho
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr,lprs
c
      dimension so(ngauss),ph(ngauss)
      common /bconf /  no,nxyzb(nobx,3)
      common /bosonwf/ psi(nobx,ngauss)
      common /couplf/ ff(ngauss,4,2)
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /fields/ phi(ngauss,4)
      common /iterat/  si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/  zero,one,two,half,third,pi
      common /optopt/  itx,icm,icou,ipc,inl,idd
      common /propag/  gg(nobx,nobx,4)
      common /physco/  hbc,alphi,r0
      common /tapes /  l6,lin,lou,lwin,lwou,lplo
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
c
c
      if (ipc.eq.1) return
c
      if (lpr) then
      write(l6,*) '****** BEGIN FIELD *********************************'
      endif
c

      do imes = 1,4
         do i = 1,ngauss
            so(i) = ff(i,imes,1)*ro(i,imes)
	     enddo
         call gordon(no,psi,gg(1,1,imes),so,phi(1,imes))
      enddo

      if (lpr) then
      write(l6,*) '****** END FIELD ***********************************'
      endif
c
      return
C-end-FIELD
      end
c======================================================================c

      subroutine gordon(no,psi,gg,so,phi)

c======================================================================c
C
C     SOLUTION OF THE KLEIN-GORDON-EQU. BY EXPANSION IN OSCILLATOR
c     imes: number of the meson
c     so:   source
c     phi:  meson field
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      dimension so(ngauss),phi(ngauss)
      dimension gg(nobx,nobx),psi(nobx,ngauss)
      dimension rn(nobx),sn(nobx)
c
      common /gaucor/ wdcor(ngauss)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      do i = 1,ngauss
         so(i) = so(i) * wdcor(i)
      enddo   ! i
      do n = 1,no
         s = zero
         do i = 1,ngauss
            s = s + psi(n,i) * so(i)
         enddo   ! i
         rn(n) = s
      enddo   ! n
c
      do n1 = 1,no
         s = zero
         do n2 = 1,no
            s = s + gg(n1,n2) * rn(n2)
         enddo   ! n2
         sn(n1) = s
      enddo   ! n1
c
      do i = 1,ngauss
         s = zero
         do n = 1,no
            s = s + psi(n,i) * sn(n)
         enddo   ! n
         phi(i) = s
      enddo   ! i
c
      return
c-end-GORDON
      end
c======================================================================c

      subroutine greemes(lpr)

c======================================================================c
c
C     calculation of the meson-propagator
c
c     calculates the meson-propagators GG
c     DD is the Laplace operator in oscillator space
c     GG = m**2/(-DD+m**2)
c
c
c     imes = 1:   sigma
c            2:   omega
c            3:   delta
c            4:   rho
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character*20 text
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /couplg/ ggmes(4),lmes(4)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /propag/  gg(nobx,nobx,4)
      common /bconf / no,nxyzb(nobx,3)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /work  / dd(nobx,nobx),gi(nobx,nobx)

c
      if (ipc.eq.1) return
c
      if (lpr) then
      write(l6,*) '****** BEGIN GREEMES *******************************'
      endif
c
c
c     calculation of DD (Laplacian)
      ax2 = one/(b0*bx)**2
      ay2 = one/(b0*by)**2
      az2 = one/(b0*bz)**2
      fac = one/(b0**3*bx*by*bz)
      do i2 = 1,no
         nx2 = nxyzb(i2,1)
         ny2 = nxyzb(i2,2)
         nz2 = nxyzb(i2,3)
         dd(i2,i2) = -(ax2*(nx2+half)+ay2*(ny2+half)+az2*(nz2+half))
c
         do i1 = 1,i2-1
            nx1 = nxyzb(i1,1)
            ny1 = nxyzb(i1,2)
            nz1 = nxyzb(i1,3)
            t   = zero
c
c           x*x
            if (ny1.eq.ny2.and.nz1.eq.nz2) then
               if (nx2.eq.nx1+2) t = ax2*half*sq(nx1+1)*sq(nx2)
               if (nx2.eq.nx1-2) t = ax2*half*sq(nx1)*sq(nx2+1)
            endif
c
c           y*y
            if (nx1.eq.nx2.and.nz1.eq.nz2) then
               if (ny2.eq.ny1+2) t = ay2*half*sq(ny1+1)*sq(ny2)
               if (ny2.eq.ny1-2) t = ay2*half*sq(ny1)*sq(ny2+1)
            endif
c
c           z*z
            if (nx1.eq.nx2.and.ny1.eq.ny2) then
               if (nz2.eq.nz1+2) t = az2*half*sq(nz1+1)*sq(nz2)
               if (nz2.eq.nz1-2) t = az2*half*sq(nz1)*sq(nz2+1)
            endif
c
            dd(i1,i2) = t
            dd(i2,i1) = t
         enddo   ! i1
      enddo   ! i2
c
      do imes = 1,4
         if (lmes(imes).eq.1) then
            f = one/ames(imes)**2
            do n = 1,no
               do k = 1,no
                  gi(n,k)      = -dd(n,k)*f
                  gg(n,k,imes) = zero
               enddo   ! k
               gi(n,n)      = one + gi(n,n)
               gg(n,n,imes) = fac
            enddo   ! n
            call lingd(nobx,nobx,no,no,gi,gg(1,1,imes),d,ifl)
c
         endif   ! lmes
      enddo   ! imes
c
c
      if (lpr) then
      write(l6,*) '****** END GREEMES *********************************'
      endif
c
      return
c-end-GREEMES
      end

c=====================================================================c

      subroutine forces(lpr)

c=====================================================================c
c
c---- options
c     center of mass:
c        icm = 0   hb0 = hbc**2/(2*amu) ecm = 3/4*hom
c              1   hb0 = hb0*(1-1/A)    ecm = 3/4*hom
c              2   hb0 = hb0            ecm = <P**2>/2M
c
c     model-type   DD     density-dependent meson-coupling
c                  PC     point-coupling
c
c---------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      character parname*10                      ! common partyp
c
      common /dforce/ a_s,a_v,a_ts,a_tv,b_s,b_v,b_ts,b_tv,
     &                c_s,c_v,c_ts,c_tv,d_s,d_v,d_ts,d_tv,dsat
      common /masses/ amu,amsig,amome,amdel,amrho
      common /coupld/ ddsig,ddome,dddel,ddrho
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /tmrpar/ gl(2),gal
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo

      if (lpr) then
      write(l6,*) '****** BEGIN FORCES ********************************'
      endif

c===============================================================
      if (parname.eq.'DD-ME2') then
c---------------------------------------------------------------
         amu    =  939.d0                 ! MeV
         amsig  =  550.1238d0             ! MeV
	     amome  =  783.d0                 ! MeV
	     amrho  =  763.d0                 ! MeV
         gsig   =   10.5396d0
         gome   =   13.0189d0
	     grho   =    3.6836d0
         b_s    =    1.0943d0
	     c_s    =    1.7057d0
         c_v    =    1.4620d0
         a_tv   =    0.5647d0
         b_tv   =    0.d0
         c_tv   =    0.d0
         d_tv   =    0.d0
         d_s  = one/sqrt(3.d0*c_s)
         a_s  = (one+c_s*(one+d_s)**2)/(one+b_s*(one+d_s)**2)
         d_v  = one/sqrt(3.d0*c_v)
         facs =two*a_s*(b_s-c_s)*(one-3.d0*c_s*(one+d_s)**2)/
     &                          (one+c_s*(1+d_s)**2)**3
         faco = (one-3.d0*c_v*(one+d_v)**2)/(one+c_v*(1+d_v)**2)**3
         x    = facs/(two*faco)
         fac1 = x+c_v*(one+c_v*(one+d_v)**2)
         fac2 = one+c_v*(one+d_v)**2-x*(one+d_v)**2
         b_v = fac1/fac2
         a_v =(one+c_v*(one+d_v)**2)/(one+b_v*(one+d_v)**2)

         a_ts   = zero
         b_ts   = zero
         c_ts   = zero
         d_ts   = zero

         dsat   =    0.152d0
	 icm    =  2
	 idd    =  1
	 ipc    =  0

c===============================================================
      else if (parname.eq.'DD-PC1') then
c---------------------------------------------------------------
c        G(x) = a + (b + c*x) * exp(-d*x)
c---------------------------------------------------------------

         dsat   =  0.152d0              ! fm^-3
         amu    =  939.d0               ! MeV

c        scalar-isoscalar
         a_s    = -10.0462d0           ! fm^-2
         b_s    =  -9.1504d0           ! fm^-2
         c_s    =  -6.4273d0           ! fm^-2
         d_s    =  +1.3724d0

c        vector-isoscalar
         a_v    =  +5.9195d0           ! fm^-2
         b_v    =  +8.8637d0           ! fm^-2
         c_v    =   0.00000d0           ! fm^-2
         d_v    =  +0.6584d0

c        scalar-isovector
         a_ts   =   zero
         b_ts   =   zero
         c_ts   =   zero
         d_ts   =   zero

c        vector-isovector
         a_tv   =   0.0000d0            ! fm^-2
         b_tv   =   1.8360d0           ! fm^-2
         c_tv   =   0.0000d0            ! fm^-2
         d_tv   =   0.6403d0
c
c----- derivative terms                 ! MeV^-4
         ddsig  =  -0.8149d0
         ddome  =   zero
         dddel  =   zero
         ddrho  =   zero
c
c----------------------------------------------------
c
c        gg = 1
         ggsig = one
         ggome = one
         ggdel = one
         ggrho = one
         amsig=zero
         amome=zero
         amdel=zero
         amrho=zero

         icm    = 2
         idd    = 2
         ipc    = 1
      else
          stop 'Wrong type of force'
      endif   !  DD-PC1
c===============================================================
c
c
c---- masses in units of fm**(-1)
      amu      = amu/hbc
      amsig    = amsig/hbc
      amome    = amome/hbc
      amdel    = amdel/hbc
      amrho    = amrho/hbc
c
      if (ipc.eq.0) then
         ggsig = -(gsig/(amsig+1.d-10))**2
         ggome = +(gome/(amome+1.d-10))**2
         ggdel = -(gdel/(amdel+1.d-10))**2
         ggrho = +(grho/(amrho+1.d-10))**2
         if (abs(gsig).gt.1.d-5) lmes(1) = 1
         if (abs(gome).gt.1.d-5) lmes(2) = 1
         if (abs(gdel).gt.1.d-5) lmes(3) = 1
         if (abs(grho).gt.1.d-5) lmes(4) = 1
      endif

      gal = 0.415d0
      gl0 = -728.d0
      do it = 1,2
          gl(it) = gl0
      enddo   ! it


c---- printout of force:
      if (lpr) call pripar
c
c
      if (lpr) then
      write(l6,*) '****** END FORCES **********************************'
      endif
c
      return
c-end-forces
      end

c=====================================================================c

      subroutine pripar

c=====================================================================c
c
c     prints parameters of the Lagrangian
c
c---------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
c
      character parname*10                     ! common partyp
c
      common /partyp/ parname
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /coupld/ ddmes(4)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /tmrpar/ gl(2),gal
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      write(l6,'(/,a,a10)') ' NEDF Parameters: ',parname
c
      if (ipc.eq.0) then
	 write(l6,'(a,f8.3)') ' AMU   = ',amu*hbc
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' msig  = ',amsig*hbc,
     &                 '  gsig = ',gsig,'  Gsig = ',ggmes(1)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mome  = ',amome*hbc,
     &                 '  gome = ',gome,'  Gome = ',ggmes(2)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mdel  = ',amdel*hbc,
     &                 '  gdel = ',gdel,'  Gdel = ',ggmes(3)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mrho  = ',amrho*hbc,
     &                 '  grho = ',grho,'  Grho = ',ggmes(4)
c
      elseif (ipc.eq.1) then
         write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                         'vecsca','vecvec'
         write(l6,'(a,4f10.6)') ' GG   = ',ggmes
         write(l6,'(a,4f10.6)') ' DD   = ',ddmes
      else
         stop 'in PRIPAR: ipc not properly defined'
      endif   ! ipc=1
      write(l6,*) ' '
      write(l6,'(a)') ' Density dependence parameters:'
      write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                   'vecsca','vecvec'
      write(l6,'(a,4f10.6)') ' a    = ',(a_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' b    = ',(b_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' c    = ',(c_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' d    = ',(d_m(m),m=1,4)
      write(l6,'(a,f10.6)') ' dsat = ',dsat

      write(l6,*) ' '
      write(l6,'(a)') ' TMR pairing: Tian,Ma,Ring, PRB 676, 44 (2009)'
      write(l6,'(a,f11.5,a)') ' TMR pairing strength    gl =: ',gl(1),
     &              ' [MeV*fm^3]'
      write(l6,'(a,f11.5,a)') ' TMR pairing width        a =: ',
     &                           sqrt(gal),' [fm]'
      return
C-end-PRIPAR
      end

C======================================================================c

      subroutine gamma()

c======================================================================c
c
c     calculats the Dirac-Matrix in the Hartee-equation
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'

      character tb*5                                            ! blokap
c
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(ngauss,2),vms(ngauss,2),vpstot(ngauss,2)
      common /hamham/ hh(nhhx,nb2x)
      common /single/ sp(nfgx,nbx)
c
      emcc2 = 2*amu*hbc
      f = hbc/b0
c
      do it=1,itx
         do ib=1,nb
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            i0f = ia(ib,1)
            i0g = ia(ib,2)
            m   = ib+(it-1)*nbx
c
c---- sigma*gradient
            do i2 = 1,nf
               do i1 = 1,ng
                  hh(nf+i1+(i2-1)*nh,m) = f*sp(i1+(i2-1)*ng,ib)
               enddo
            enddo
c
c---- V+S and V-S
            call pot(i0f,nh,nf,vpstot(1,it),hh(1,m))
            call pot(i0g,nh,ng,vms(1,it),hh(nf+1+nf*nh,m))
c
c---- shift by emc2
            do i = nf+1,nh
               hh(i+(i-1)*nh,m) = hh(i+(i-1)*nh,m) - emcc2
            enddo
c
c---- symmetrize HH
            do i2 = 1,nh
               do i1 = i2+1,nh
                  hh(i2+(i1-1)*nh,m) = hh(i1+(i2-1)*nh,m)
               enddo
            enddo
         enddo
      enddo
c
      return
C-end-GAMMA
      end
c======================================================================c

      subroutine pot(i0,nh,n,v,aa)

c======================================================================c
c
c     calculates the potentials V+S and V-S in an axial oscillator basis
c
c     i0   initial point for the quantum number array
c     n    dimension
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tt*8                                            ! bloqua
c
      dimension aa(nh,nh),v(0:ngh,0:ngh,0:ngh)
      dimension f1(0:ngh,0:ngh,0:kmax2),f2(0:ngh,0:kmax2,0:kmax2)
      dimension four(0:kmax2,0:kmax2,0:kmax2)
c
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /gfviv / iv(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /gaussh/  xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /hermite/ hpol(0:kmax2,0:ngh),dhpol(0:kmax2,0:ngh),
     &                 ddhpol(0:kmax2,0:ngh)
      common /herexp/  C00(0:kmax2,0:kmax,0:kmax)
c
      do kx=0,kmax2,2
         do ihz=0,ngh
            do ihy=0,ngh
               sum=0.d0
               do ihx=0,ngh
                  sum = sum+hpol(kx,ihx)*ph(ihx)*v(ihx,ihy,ihz)
               enddo !ihx
               f1(ihy,ihz,kx) = sum
            enddo !ihy
         enddo !ihz

         do ky=0,kmax2,2
            do ihz=0,ngh
               sum=0.d0
               do ihy=0,ngh
                  sum = sum+hpol(ky,ihy)*ph(ihy)*f1(ihy,ihz,kx)
               enddo
               f2(ihz,ky,kx) = sum
            enddo !ihz
            do kz=0,kmax2,2
               sum=0.d0
               do ihz=0,ngh
                  sum = sum + hpol(kz,ihz)*ph(ihz)*f2(ihz,ky,kx)
               enddo
               four(kz,ky,kx) = sum
            enddo !kz
         enddo !ky
      enddo !kx

      do i2=1,n
         nx2 = nxyz(i0+i2,1)
         ny2 = nxyz(i0+i2,2)
         nz2 = nxyz(i0+i2,3)
         do i1=i2,n
            nx1 = nxyz(i0+i1,1)
            ny1 = nxyz(i0+i1,2)
            nz1 = nxyz(i0+i1,3)
            nx = nx1+nx2
            ny = ny1+ny2
            nz = nz1+nz2
            if (iv(nx).eq.1.and.iv(ny).eq.1.and.iv(nz).eq.1) then
               iny = iabs(ny1-ny2)/2
               my  = iv(iny)
               sumx = 0.d0
               do kx=0,nx,2
                  sumy = 0.d0
                  do ky=0,ny,2
                     sumz = 0.d0
                     do kz=0,nz,2
                        sumz=sumz+C00(kz,nz1,nz2)*four(kz,ky,kx)
                     enddo
                     sumy=sumy+C00(ky,ny1,ny2)*sumz
                  enddo !ky
                  sumx=sumx+C00(kx,nx1,nx2)*sumy
               enddo !kx
               aa(i1,i2) = my*sumx
            else
               aa(i1,i2) = 0.d0
            endif
         enddo !i1
      enddo ! i2
c
      return
c-end-POT
      end

c======================================================================c

      subroutine gaupol(lpr)

c======================================================================c
c
c     BERECHNET DIE  HERMITE-POLYNOME
c     numn...maximaler (n+1)-Wert in der Basis
c     qh(0:nmax,ngh)...Hermitepolynom mal Normierungskonstante
c            qh(n,i)=n(n) * h(n,zeta)
c     qh1(0:nmax,ngh)...speziell definierte erste Ableitung der
c            Hermitepolynome  qh1(n,i)=
c            n(n) * (2*n*h(n-1,zeta)-zeta*h(n,zeta)
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
      logical lpr
c
      dimension h1(0:ngh),h2(0:ngh)
c
      common /bconf /  no,nxyzb(nobx,3)
      common /bosonwf/ psi(nobx,ngauss)
      common /gaussh/  xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /hermite/ hpol(0:kmax2,0:ngh),dhpol(0:kmax2,0:ngh),
     &                 ddhpol(0:kmax2,0:ngh)
      common /hermite2/ hpol2(0:kmax2,0:ngh),dhpol2(0:kmax2,0:ngh)
      common /herexp/  C00(0:kmax2,0:kmax,0:kmax)
      common /gfvsq /  sq(0:igfv)
      common /gfviv /  iv(0:igfv)
      common /mathco/  zero,one,two,half,third,pi
      common /polher/  qh(0:nmax,0:ngh),qh1(0:nmax,0:ngh)
      common /tapes /  l6,lin,lou,lwin,lwou,lplo
c
      data eps/1.e-9/
c
      if (lpr) then
      write(l6,*) '****** BEGIN GAUPOL *******************************'
      endif
c
c
      do i = 0,ngh
         call hermio(nmax,xh(i),qh(0,i),qh1(0,i))
         w = exp(-half*xh(i)**2)
         do n = 0,nmax
            qh(n,i)  = qh(n,i)*w
            qh1(n,i) = qh1(n,i)*w
         enddo   ! n
      enddo   ! i
      do i=0,ngh
         call hermio(kmax2,xh(i),hpol(0,i),dhpol(0,i))
         x  = xh(i)
         x2 = x*x
         ddhpol(0,i) = 2*(2*x2-1)*hpol(0,i)
         ddhpol(1,i) = -4*x*sq(2)*sq(1)*hpol(0,i)+2*(2*x2-1)*hpol(1,i)
         do n=2,kmax2
            ddhpol(n,i) = 2*sq(n)*sq(n-1)*hpol(n-2,i)
     &                   -4*x*sq(2)*sq(n)*hpol(n-1,i)
     &                   +2*(2*x2-1)*hpol(n,i)
         enddo
      enddo
      do i=0,ngh
         x=sq(2)*xh(i)
         call hermio(kmax2,x,hpol2(0,i),dhpol2(0,i))
      enddo
      do n1=0,kmax
         do i=0,ngh
            h1(i)=hpol(n1,i)
         enddo
         do n2=n1,kmax
            do i=0,ngh
               h2(i)=hpol(n2,i)
            enddo
            do k=0,kmax2
               sum=0.d0
               if(iv(k) .eq. iv(n1+n2)) then
                  do i=1,ngh
                     sum = sum + h1(i)*h2(i)*hpol(k,i)*ph(i)
                  enddo
               endif
               C00(k,n1,n2) = sum
               C00(k,n2,n1) = sum
            enddo
         enddo
      enddo
c
c======================================================================c
c     mixed boson wavefunctions
c======================================================================c
c
c---  w.f of boson for densities
      do n = 1,no
         nx = nxyzb(n,1)
         ny = nxyzb(n,2)
         nz = nxyzb(n,3)
         do iz = 0,ngh
         do iy = 0,ngh
         do ix = 0,ngh
            i=1+ix+iy*(ngh+1)+iz*(ngh+1)*(ngh+1)
            psi(n,i) = qh(nx,ix)*qh(ny,iy)*qh(nz,iz)
         enddo   ! ix
         enddo   ! iy
         enddo   ! iz
      enddo   ! no
c
c======================================================================c
c     fermion wavefunctions with sqrt(wh)
c======================================================================c
      do i = 0,ngh
         x = sqrt(wh(i))
         do n = 0,nmax
            qh(n,i)  = qh(n,i)*x
            qh1(n,i) = qh1(n,i)*x
         enddo   ! n
      enddo   ! i
c
c
c     test of orthogonality
c--------------------------
      if (lpr) then
      eps1=zero
      eps2=zero
      write(l6,*) 'Number of gauss-hermite meshpoints ngh=',ngh
      write(l6,*) 'from xh(1) = ',xh(1),' to ',' xh(ngh) = ',xh(ngh)

      write(l6,*) ' Test of Orthogonality for Hermite-Polynomials'
      do n1 = 0,nmax
      do n2 = n1,nmax,2
         s = zero
         do ih = 1,ngh
            s = s + qh(n1,ih)*qh(n2,ih)
         enddo   ! ih
         if (n1.eq.n2) then
            if (abs(s-one).gt.eps1) eps1=abs(s-one)
         else
            if (abs(s).gt.eps2) eps2=abs(s)
         endif
      enddo   ! n2
      enddo   ! n1
      write(l6,*) 'max. dev. in orthogonality: ',eps1
      write(l6,*) 'max. dev. in normalization: ',eps2
      endif   ! lpr
c
      if (lpr) then
      write(l6,*) '****** END GAUPOL *********************************'
      endif
c
      return
c-end-herpol
      end
c======================================================================c

      subroutine hermio(n,x,u,v)

c======================================================================c
c     Calculates normalized oscillator functions of one dimension
c     phi_n_(x)  =  u(n+1) * exp(-x*x/2)
c     and their derivatives
c     d/dx phi_n_(x)  =  v(n+1) * exp(-x*x/2)
c     recursion from below
c
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      dimension u(0:n),v(0:n)
c
      common /gfvsqi/ sqi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
c
      w4pi=pi**(-0.25d0)
c
      u(0) =  w4pi
      v(0) = -w4pi*x
      if (n.eq.0) return
      s2 = sqi(2)
      x2 = x*x
      do i = 1,n
         s = s2*sqi(i)
         u(i) = s*(x*u(i-1) - v(i-1))
         v(i) = s*(x*v(i-1) + (dfloat(i+i)-x2)*u(i-1))
      enddo
c
      return
c-end-hermio
      end

c======================================================================c

       subroutine gaush(lpr)

c======================================================================c
c
c     Gauss-Hermite integration data
c     ------------------------------
c     for integration from minus to plus infinity
c                     or from 0 to infinity
c
c     ph  =  wh * exp(-xh**2)
c
c     whh=ph
c
c     \int_-\infty^+\infty  f(z) exp(-z**2) dz  =   \sum_i f(xh(i)) ph(i)
c    possible alternative
c     \int_-\infty^+\infty  f(z) dz             =   \sum_i f(xh(i)) wh(i)
c
c     fak = 1  --- for integration from 0 to +infinity
c     fak = 2  --- for integration from -infinity to +infinity
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      dimension xhh(2*ngh),whh(2*ngh)
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c

      if (lpr)
     &write(l6,*) ' ****** BEGIN GAUSH *******************************'
c
c
      call gauher(xhh,whh,2*ngh)
      do i = 1,ngh
         xh(i) = xhh(ngh+1-i)
         ph(i) = whh(ngh+1-i)
         wh(i) = ph(i)*dexp(xh(i)*xh(i))
      enddo
      xh(0)=1.d-10
      ph(0)=1.d-10
      wh(0)=1.d-10

      if(lpr) write(l6,100) ngh
C
      if (.not.lpr) return
      write(l6,101)
c
      write(l6,*) ' xh'
      write(l6,105) (xh(i),i=1,ngh)
c
      write(l6,*) ' wh'
      write(l6,105) (wh(i),i=1,ngh)
c
      write(l6,*) ' ph'
      write(l6,105) (ph(i),i=1,ngh)
c
  100 format('  GAUSH:  G-H-Integration  ngh =',i3)
  101 format(1x,36(1h-))
  105 format(3e19.11)
c
      if (lpr)
     &write(l6,*) ' ****** END GAUSH *********************************'
      return
c-end-GAUSH
      end
c
c==========================================================================
c======================================================================c
c
      subroutine gdd(lpr)
c
c----------------------------------------------------------------------c
c
c     calculates  density-dependence of the coupling constants
c
c======================================================================c
      implicit real*8(a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      dimension gf(4)
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /couplf/ ff(ngauss,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo

      fun1(x,a,b,c,d) = a*(1+b*(x+d)**2)/(1+c*(x+d)**2)
      dun1(x,a,b,c,d) = 2*a*(b-c)*(x+d)/(1+c*(x+d)**2)**2

      fun2(x,a,b,c,d) = exp(-a*(x-one))
      dun2(x,a,b,c,d) = -a*exp(-a*(x-one))

      fun3(x,a,b,c,d) = a+(b+c*x)*exp(-d*x)
      dun3(x,a,b,c,d) = c*exp(-d*x)-d*(b+c*x)*exp(-d*x)

c
      if (lpr) then
      write(l6,*) '****** BEGIN GDD ***********************************'
      endif

c
      if (ipc.eq.0) then
         do i = 1,ngauss
            x=ro(i,2)/dsat
            do m = 1,2
               ff(i,m,1) = fun1(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun1(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
            do m = 3,4
               ff(i,m,1) = fun2(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun2(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
         enddo   ! i
      else if(ipc.eq.1) then
         do i = 1,ngauss
            x=ro(i,2)/dsat
            do m = 1,4
               ff(i,m,1) = fun3(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun3(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
         enddo
      endif
c
      if (lpr) then
      write(l6,*) '****** END GDD *************************************'
      endif

      return
c-end-GDD
      end

c======================================================================c

      subroutine inout(is,lpr)

c======================================================================c
c
c     is = 1: reads meson-fields from tape
c          2: writes meson-fields  to tape
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character*2 nucnam,nucnam1
      character*8 parnam1
      character*15 text
      character tp*1,tis*1,tit*8,tl*1                    ! common textex
      character tb*8,txb*5,tx*8                          ! common texblo
      character parname*10                     ! common partyp
      character parpair*10                               ! common partpp
c
      dimension ga1(2),g01(2),del1(2),spk1(2),dec1(2),tz1(2)
      dimension npr1(2)
c
      common /baspar/ hom,hb0,b0
      common /fermi / ala(2),tz(2)
      common /cstr2/ q0c,q2c,c0,c2,alaq0,alaq2
      common /initia/ vin,rin,ain,inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(ngauss,2),vms(ngauss,2),vpstot(ngauss,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      if (is.eq.1.and.inin.ne.0) return
c
      write(l6,*) ' ****** BEGIN INOUT ********************************'
c
c
c
c
c---- reading of meson fields from tape:
c-------------------------------------
      if (is.eq.1) then
         open(lwin,file='dirhb.wel',status='old')
         read(lwin,103) ala,alaq0,alaq2
c
c------- reading of the potentials
         read(lwin,*)
         read(lwin,101) vms
         read(lwin,*)
         read(lwin,101) vpstot
         close(lwin)
         write(l6,*) ' potentials read from tape ','dirhb.wel'
c
c---- writing to tape:
      else
         open(lwou,file='dirhb.wel',status='unknown')
         write(lwou,104) 'Lambda:   ',ala,alaq0,alaq2
c
c------- writing of the potentials
         write(lwou,*)   'VMS:'
         write(lwou,101) vms
         write(lwou,*)   'VPS:'
         write(lwou,101) vpstot
         close(lwou)

         write(l6,*) ' potentials written to tape dirhb.wel'
      endif
c
  100 format(1x,a2,8i4)
  101 format(4e20.12)
  103 format(10x,4f12.6)
  104 format(a,4f12.6)
c
      write(l6,*) ' ****** END INOUT **********************************'
      return
c-end-INOUT
      end

c=====================================================================c

      subroutine plot(lpr)

c=====================================================================c
C
C     prepares plot of densities in coordinate space
C
c---------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
      logical lpr
c
      dimension pn(nox1),funk(0:ngh,0:ngh)
      dimension zp(0:ngh,0:ngh,3)
c
      common /baspar/ hom,hb0,b0
      common /basdef/ beta0,gamma0,bx,by,bz
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /dens/ ro(ngauss,4),dro(ngauss,4)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /physco/ hbc,alphi,r0
      common /gaussh/ xh(0:NGH),wh(0:NGH),ph(0:NGH)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
c
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN PLOT ********************************'
c
c     number of points for the plot
      mxplx  = 80
      mxply  = 80
      mxplz  = 80
c     plot step in (fm)
      stplx = 0.2
      stply = 0.2
      stplz = 0.2

c---- plot for densities:

      open(lplo,file='dirhb-xz.plo',status='unknown')
      iy=0
      do ix=0,ngh
         do iz=0,ngh
            ih=1+ix + iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
            funk(ix,iz)=ro(ih,2)
         enddo
      enddo
      call splin2(xb,zb,NGH,NGH,funk,zp)
      x=-8.d0
      do ix=0,mxplx
         z=-8.d0
         do iz=0,mxplz
            call splint2(abs(x),abs(z),NGH,NGH,xb,zb,funk,zp,rv)
            write(lplo,100) x,z,rv
            z=z+stplz
         enddo
         x=x+stplx
      enddo
      close(lplo)

      open(lplo,file='dirhb-xy.plo',status='unknown')
      iz=0
      do ix=0,ngh
         do iy=0,ngh
            ih=1+ix + iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
            funk(ix,iy)=ro(ih,2)
         enddo
      enddo
      call splin2(xb,yb,NGH,NGH,funk,zp)

      x=-8.d0
      do ix=0,mxplx
         y=-8.d0
         do iy=0,mxply
            call splint2(abs(x),abs(y),NGH,NGH,xb,yb,funk,zp,rv)
            write(lplo,100) x,y,rv
            y=y+stply
         enddo
         x=x+stplx
      enddo
      close(lplo)

      open(lplo,file='dirhb-yz.plo',status='unknown')
      ix=0
      do iy=0,ngh
         do iz=0,ngh
            ih=1+ix + iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
            funk(iy,iz)=ro(ih,2)
         enddo
      enddo
      call splin2(yb,zb,NGH,NGH,funk,zp)

      y=-8.d0
      do iy=0,mxply
         z=-8.d0
         do iz=0,mxplz
            call splint2(abs(y),abs(z),NGH,NGH,yb,zb,funk,zp,rv)
            write(lplo,100) y,z,rv
            z=z+stplz
         enddo
         y=y+stply
      enddo
      close(lplo)

  100       format(2f10.3,f15.6)
      close(lplo)
c
      if (lpr)
     &write(l6,*) ' ****** END PLOT **********************************'
      return
C-end-PLOT
      end

c======================================================================c

      subroutine poten(lpr)

c======================================================================c
c
c     CALCULATION OF THE POTENTIALS AT GAUSS-MESHPOINTS
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      dimension glt(4)
c
      common /baspar/ hom,hb0,b0
      common /cstr3/ vc(ngauss)
      common /cstr1/ betac,gammac,cqad,icstr
      common /coulmb/ cou(ngauss),drvp(ngauss)
      common /coupld/ ddmes(4)
      common /couplf/ ff(ngauss,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /dens  / ro(ngauss,4),dro(ngauss,4)
      common /fields/ phi(ngauss,4)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(ngauss,2),vms(ngauss,2),vpstot(ngauss,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      if (lpr) then
      write(l6,*) '****** BEGIN POTEN *********************************'
      endif
c
      xmi = xmix
c      xmi = one
c
      do i = 1,ngauss
c
c------- meson-fields
         if (ipc.eq.0) then
	    do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*phi(i,m)
	    enddo   ! m
c
c           rearangement field
            if (idd.ge.1) then
	       re = zero
	       do m = 1,4
                  re = re + ggmes(m)*ff(i,m,2)*phi(i,m)*ro(i,m)
	       enddo   ! m
               glt(2) = glt(2) + re
            endif
c
c------- point-coupling models
	 elseif (ipc.eq.1) then
	    do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*ro(i,m)
	    enddo   ! m
c        derivative terms
	    do m = 1,4
               glt(m) = glt(m) + ddmes(m)*dro(i,m)
	    enddo   ! m
c
c           rearangement field
            if (idd.ge.1) then
	       re = zero
	       do m = 1,4
                  re = re + ggmes(m)*ff(i,m,2)*ro(i,m)**2
	       enddo   ! m
	       glt(2)=glt(2)+half*re
            endif
c
         else
	    stop 'in POTEN: ipc not properly defined'
         endif   ! ipc
c
         s1 = glt(1) - glt(3)                 ! neutron scalar
         s2 = glt(1) + glt(3)                 ! proton  scalar
         v1 = glt(2) - glt(4)                 ! neutron vector
         v2 = glt(2) + glt(4) + cou(i)        ! proton  vector
c
c------- constraining potential
         if (icstr.gt.0) then
            v1 = v1 + vc(i)/hbc
            v2 = v2 + vc(i)/hbc
         endif   ! icstr

         vps1 = (v1 + s1)*hbc
         vps(i,1) = vps1
c
         vms1 = (v1 - s1)*hbc
         vms(i,1) = vms1
c
         vps2 = (v2 + s2)*hbc
         vps(i,2) = vps2
c
         vms2 = (v2 - s2)*hbc
         vms(i,2) = vms2

c------- constraining potential
         vpstot(i,1) = vps(i,1)! + vc(i)
         vpstot(i,2) = vps(i,2)! + vc(i)
c
c
      enddo   ! i

      if (lpr) then
      write(l6,*) '****** END POTEN ***********************************'
      endif
c
      return
c-end-POTEN
      end

c======================================================================c

      subroutine prep(lpr)

c======================================================================c
c
c     preparations
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
      character tb*5                                            ! blokap
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /gaucor/ wdcor(ngauss)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /cstr1/ betac,gammac,cqad,icstr
      common /cstr2/ q0c,q2c,c0,c2,alaq0,alaq2
c
c
      if (lpr) then
      write(l6,*) '****** BEGIN PREP **********************************'
      endif
c
c
c---- basis parameters
      hb0 = hbc/(two*amu)
      hom = 41.0*amas**(-third)
      if (icm.eq.1) hb0 = hb0*(one - one/amas)
      write(l6,*) ' '
      write(l6,'(a,f10.6)') ' hom   =',hom
      write(l6,'(a,f10.6)') ' hb0   =',hb0
      if (b0.le.0.0) then
          b0 = sqrt(two*hb0/hom)
      endif
c
      w  = beta0*half*dsqrt(5.d0/(4.d0*pi))
      bx = exp(w*cos((gamma0-120.d0)*pi/180.d0))
      by = exp(w*cos((gamma0+120.d0)*pi/180.d0))
      bz = exp(w*cos(gamma0*pi/180.d0))
c
      write(l6,'(a,f10.6)') ' b0    =',b0
      write(l6,'(a,f10.6)') ' beta0 =',beta0
      write(l6,'(a,f10.6)') ' gamm0 =',gamma0
      write(l6,'(a,f10.6)') ' bx    =',bx
      write(l6,'(a,f10.6)') ' by    =',by
      write(l6,'(a,f10.6)') ' bz    =',bz
c
c---- coordinates in fm
      do ih = 0,ngh
         xb(ih) = xh(ih)*b0*bx
         yb(ih) = xh(ih)*b0*by
         zb(ih) = xh(ih)*b0*bz
c
c---- integration from 0 to infinity
         wh(ih) = 2*wh(ih)
         ph(ih) = 2*ph(ih)
      enddo  ! ngh
      do ix = 0,ngh
      do iy = 0,ngh
      do iz = 0,ngh
         ih= 1+ix+ iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
         wdcor(ih) = b0**3*bx*by*bz*wh(ix)*wh(iy)*wh(iz)
      enddo   ! iz
      enddo   ! iy
      enddo   ! ix
c
c
c---- constraining field
      if (icstr.ge.1) then
         gam=gammac*pi/180.0
         r02=1.2**2*amas**(2.d0/3.d0)
         a0c=betac*cos(gam)
         a2c=a0c*tan(gam)/sqrt(2.d0)
         q0c=3*amas*r02*a0c/4/pi
         q2c=3*amas*r02*a2c/4/pi
         f0=sqrt(16*pi/5.d0)
         f2=sqrt(32*pi/15.d0)
         q0c=q0c*f0
         q2c=q2c*f2
         c0=cqad*hom/(f0*hbc*b0**2)
         c2=cqad*hom/(f2*hbc*b0**2)
      endif
c
   10 if (itx.eq.1) icou = 0
      if(icou.eq.1) write(l6,*) 'With Coulomb force'
      if(icm.eq.2) write(l6,*) 'With microscopic c.m. correction'

c
  100 format(a,4f10.6)
  101 format(a,i4)
  102 format(4f15.8)
c
      if (lpr) then
      write(l6,*) '****** END PREP ************************************'
      endif
c
      return
c-end PREP
      end

c======================================================================c

      subroutine reader(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character*10 inputf,outputf,outputw,outdel   ! ttttt1
      character parname*10                            ! partyp
      character nucnam*2                                        ! nucnuc
c
      CHARACTER*8 date
      CHARACTER*10 time
      CHARACTER*5 zone
      INTEGER*4 VALUES(8)

      common /basdef/ beta0,gamma0,bx,by,bz
      common /basnnn/ n0f,n0b
      common /fermi / ala(2),tz(2)
      common /initia/ vin,rin,ain,inin,inink
      common /nucnuc/ amas,nama,nneu,npro,nucnam
      common /partyp/ parname
      common /start0/ betai,gammai
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /icalc/ ivpair
      common /pair  / del(2),spk(2),spk0(2)
      common /numsep/ npm
      common /cstr1/ betac,gammac,cqad,icstr
c
      inputf  ='dirhb.dat'
c
      if (lin.eq.0) return
      open(lin,file=inputf,status='old')
c
c
c---- Basisparameters:
      read(lin,'(10x,2i5)') n0f,n0b
c
c---- Deformation-Parameters
      read(lin,'(10x,f10.4)') beta0
      read(lin,'(10x,f10.4)') gamma0
c
c---- Initial Deformation--Parameters of the start potential
      read(lin,'(10x,f10.4)') betai
      read(lin,'(10x,f10.4)') gammai
c
c---- Initialization of wavefunctions:
      read(lin,'(10x,2i5)') inin,inink
c
c---- Initial gaps
      read(lin,'(10x,2f10.4)') del(1),del(2)
      npm=2*n0f
c
c---- Nucleus under consideration
      read(lin,'(a2,i4)') nucnam,nama
c
c---- Parameterset of the Lagrangian
      read(lin,*)
      read(lin,'(12x,a10)') parname

c
c---- quadratic constraint
      read(lin,*)
      read(lin,'(10x,i5)') icstr
      if(icstr .ne. 0 .and. icstr .ne. 1 .and. icstr .ne. 2)
     &     stop 'Wrong value for the icstr in READER'
      read(lin,'(10x,f10.4)') betac
      read(lin,'(10x,f10.4)') gammac
      read(lin,'(10x,f10.4)') cqad
c
      call nucleus(2,npro,nucnam)
c
c======================================================================c
      outputf ='dirhb.out'
      outputw ='dirhb.wel'
      outdel  ='dirhb.del'
c======================================================================c
c
c
      if (l6.ne.6) open(l6,file=outputf,status='unknown')
      call date_and_time( date, time, zone, values )

      write(l6,'(a)')
      write(l6,'(a)') '  ******************************************  '
      write(l6,'(a)') '  *           Fortran 77 Code              *  '
      write(l6,'(a)') '  *         Triaxial H.O. basis            *  '
      write(l6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(l6,'(a)') '  *        with density-dependent          *  '
      write(l6,'(a)') '  * meson-exchange or point coupling force *  '
      write(l6,'(a)') '  *          and separable pairing         *  '
      write(l6,'(a)') '  * -------------------------------------- *  '
      write(l6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(l6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(l6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nama,'  N = ',nama-npro,'  Z = ',npro
      write(l6,'(a,16x,a10)') '  *',parname
      write(l6,'(a)') '  ******************************************  '
c
      write(6,'(a)')
      write(6,'(a)') '  ******************************************  '
      write(6,'(a)') '  *           Fortran 77 Code              *  '
      write(6,'(a)') '  *         Triaxial H.O. basis            *  '
      write(6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(6,'(a)') '  *        with density-dependent          *  '
      write(6,'(a)') '  * meson-exchange or point coupling force *  '
      write(6,'(a)') '  *          and separable pairing         *  '
      write(6,'(a)') '  * -------------------------------------- *  '
      write(6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nama,'  N = ',nama-npro,'  Z = ',npro
      write(6,'(a,16x,a10)') '  *',parname
      write(6,'(a)') '  ******************************************  '

c
      if (lpr) then
      write(l6,*) '****** BEGIN READER ********************************'
      endif
c
c---- Basisparameters:
      write(l6,'(a,2i5)') ' Number of oscillator shells : ',n0f,n0b
c
c---- basis deformation
      write(l6,'(a,f10.4,a,f10.4)')
     &              ' Basis deformation           :  beta0 =',beta0
      write(l6,'(a,f10.4,a,f10.4)')
     &              ' Basis deformation           :  gamma0=',gamma0
c---- initial deformation
      write(l6,'(a,f10.4,a,f10.4)')
     &              ' Initial deformation         :  betai =',betai
      write(l6,'(a,f10.4,a,f10.4)')
     &              ' Initial deformation         :  gammi =',gammai
c
c---- Initialization of wavefunctions:
      write(l6,'(a,2i5)')
     &              ' Initial wavefunctions       : ',inin,inink
c---- Initial gaps
      write(l6,'(a,2f10.4)')
     &              ' Initial pairing gaps        : ',del(1),del(2)
c
c---- mixing parameter
      write(l6,'(a,f10.4)')
     &              ' Mixing-Parameter xmix       : ',xmix
c
c---- constraints
      if (icstr.eq.0) then
          write(l6,'(a)') ' Calculation without quadrupole constraint'
      else
        write(l6,'(a,f10.4)') ' Constrained parameter beta    : ',betac
        write(l6,'(a,f10.4)') ' Constrained parameter gamma   : ',gammac
        write(l6,'(a,f10.4)') ' Spring constant               : ',cqad
      endif
c
      close(lin)
c
      nneu = nama - npro
      amas = nama
c
c---- Nucleus under consideration
      write(l6,'(a,a2,i4,a,i3,a,i3)')
     &              ' Nucleus                     : ',
     &    nucnam,nama,'  N = ',nneu,'  Z = ',npro
      tz(1) = nneu
      tz(2) = npro
c
      if (lpr) then
         write(l6,*) '****** END READER ****************************'
      endif
c
      return
c-end-READER
      end

c======================================================================c

      subroutine resu(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tb*5                                            ! blokap
      character tk*8                                            ! blolev
      character tt*8                                            ! bloqua
      character tph*1
c
      dimension ix(ntx),ez(ntx)
c
      common /baspar/ hom,hb0,b0
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /coulmb/ cou(ngauss),drvp(ngauss)
      common /fermi / ala(2),tz(2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
c
      data emaxp/50.d0/,emina/-1400.d0/
c
      if (.not.lpr) return
c
      write(l6,*) ' ****** BEGIN RESU *********************************'
c
c---- center-of-mass correction
c
      call expect(2,.true.)
c
c
c
c---- single particle energies
      do it = 1,itx
         write(l6,100) tit(it)
  100    format(//,' Quasi-particle properties: ',a,/,1x,66(1h-))
         if (lpr) write(l6,1000) ' ',' N pi ','  [ nx ny nz]',
     &                           'smax','eeqp','vv','uu','p/h'
1000     format(a2,a6,a13,2x,a9,1x,a9,1x,a8,2x,a8,2x,a3)
c
         do ib = 1,nb
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
            k1 = ka(ib,it) + 1
            k2 = ka(ib,it) + kd(ib,it)
            i0 = ia(ib,1)
c
            do k = k1,k2
c---------- search for main oscillator component
               smax = zero
               su = zero
               sv = zero
               do n = 1,nh
                  u = fguv(n,k,it)
                  v = fguv(nh+n,k,it)
                  su = su+u**2
                  sv = sv+v**2
               enddo
               if(su .gt. sv) then
                  do n = 1,nf
                     s = abs(fguv(n,k,it))
                     if (s.gt.smax) then
                        smax = s
                        imax = n
                     endif
                  enddo
                  tph = 'p'
               else
                  do n = 1,nf
                     s = abs(fguv(nh+n,k,it))
                     if (s.gt.smax) then
                        smax = s
                        imax = n
                     endif
                  enddo
                  tph = 'h'
               endif

c
c---------- printing
               if (equ(k,it).lt.60.d0)
     &         write(l6,101) k,tp(ib),tt(i0+imax),smax,
     &                    equ(k,it),su,sv,tph
  101       format(i4,' ',a1,5x,a8,8x,f5.2,3f10.3,2x,a1)
   20       enddo  ! k
         enddo   ! ib
      enddo   ! it
      write(l6,*) ' ****** END RESU ***********************************'
c
      return
c-end-RESU
      end

c======================================================================c

      subroutine singd(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions
c     in the cartesian oscillator basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
c
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /single/ sp(nfgx,nbx)
      common /vvvsep/ vnn(0:n0fx,0:n0fx,0:npmax/2,3)
      common /ppp/ oscn(0:npmax/2,3)
      common /tmrpar/ gl(2),gal
      common /numsep/ npm

c
      if (lpr) then
      write(l6,*) '****** BEGIN VPAIR ********************************'
      endif

      call vn0(1)
      call vn0(2)
      call vn0(3)

      do n2=0,n0fx
         do n1=0,n0fx
            do NN=0,npm,2
               rfac = rmosh(n1,n2,NN)
               nn1 = n2+n1-NN
               if(nn1 .ge. 0) then
                  vnn(n1,n2,NN/2,1) = rfac*oscn(nn1/2,1)
                  vnn(n1,n2,NN/2,2) = rfac*oscn(nn1/2,2)
                  vnn(n1,n2,NN/2,3) = rfac*oscn(nn1/2,3)
               endif
            enddo
         enddo
      enddo

c
c
      if (lpr) then
      write(l6,*) '****** END VPAIR **********************************'
      endif
c
      return
c-end-SINGD
      end

c====================================================================================
      subroutine vn0(nkart)
c====================================================================================
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'

      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /tmrpar/ gl(2),gal
      common /gfvfi / fi(0:igfv)
      common /gfviv / iv(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /ppp/ oscn(0:npmax/2,3)

      a = sqrt(gal)
      if(nkart.eq.1) then
         b = bx*b0
      else if(nkart.eq.2) then
         b = by*b0
      else if(nkart.eq.3) then
         b = bz*b0
      else
         stop 'Wrong value for NKART in vn0'
      endif
      f1 = (2*pi)**(-0.25d0)
      f2 = sqrt(b/(a**2+b**2))
      do n=0,npmax,2
         f3 = wf(n)
         f4 = one/2**(n/2)
         f5 = fi(n/2)
         f6 = ((a**2-b**2)/(a**2+b**2))**(n/2)
         fac = f1*f2*f3*f4*f5*f6
         oscn(n/2,nkart) = fac
      enddo

      return
c-end-VN0
      end
c====================================================================================
      real*8 function rmosh(n1,n2,N)
c====================================================================================
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'

      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfviv / iv(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /mathco/ zero,one,two,half,third,pi

      nn = n1+n2-N
      if (nn .lt. 0) then
         rmosh = 0.d0
         return
      endif

      f1 = wf(n1)*wf(n2)*wfi(N)*wfi(nn)
      f2 = one/2**((nn+N)/2)

      m1 = max(0,n1-N)
      m2 = min(n1,nn)
      rmosh = zero
      do m=m1,m2
         f3 = iv(m)
         f4 = fak(N)*fi(N-n1+m)*fi(n1-m)
         f5 = fak(nn)*fi(m)*fi(nn-m)
         rmosh = rmosh + f3*f4*f5
      enddo
      rmosh = iv(nn)*rmosh*f1*f2
      return
      end

c======================================================================c

      subroutine singf(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions
c     in the cartesian oscillator basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
c
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /single/ sp(nfgx,nbx)
c
      if (lpr) then
      lx = l6
      l6 = 6
      write(l6,*) '****** BEGIN SINGF ********************************'
      endif
c
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
c
c        SIGMA*P
c----------------
         call sigp(nf,ng,ib,sp(1,ib),lpr)
      enddo  ! ib
c
      if (lpr) then
      write(l6,*) '****** END SINGF **********************************'
      endif
c
      return
c-end-SINGF
      end
c=====================================================================c

      subroutine sigp(nf,ng,ib,aa,lpr)

c=====================================================================c
c
c     Calculates <n1|i*f*SIGma*P|n2>
c
c           ( 0         )        ip   : Parity
c           ( 0 0       )        |n1> : neg. Space Reflection
c     tt =  ( * * 0     )        |n2> : pos. Space Reflection
c           ( * * 0 0   )
c           ( * * 0 0 0 )
c
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
      dimension aa(ng,nf)
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
      w = one/sq(2)
      fx  = w/bx
      fy  = w/by
      fz  = w/bz
      i0f = ia(ib,1)
      i0g = ia(ib,2)
      do n2 = 1,nf
         nx2 = nxyz(i0f+n2,1)
         ny2 = nxyz(i0f+n2,2)
         nz2 = nxyz(i0f+n2,3)
         do n1 = 1,ng
            nx1 = nxyz(i0g+n1,1)
            ny1 = nxyz(i0g+n1,2)
            nz1 = nxyz(i0g+n1,3)
c
            s = zero
c
            mxy1=iv(nx1+ny1+1)
            my1= iv(ny1)
c
c-----      sigma_x*grad_x
            if (ny2.eq.ny1.and.nz2.eq.nz1) then
               if (nx2.eq.nx1+1) s = my1*fx*sq(nx2)
               if (nx2.eq.nx1-1) s =-my1*fx*sq(nx1)
c
c-----      sigma_y*grad_y
            elseif (nx2.eq.nx1.and.nz2.eq.nz1) then
               if (ny2.eq.ny1+1) s =-my1*fy*sq(ny2)
               if (ny2.eq.ny1-1) s =-my1*fy*sq(ny1)
c
c-----      sigma_z*grad_z
            elseif (nx2.eq.nx1.and.ny2.eq.ny1) then
               if (nz2.eq.nz1+1) s =-mxy1*fz*sq(nz2)
               if (nz2.eq.nz1-1) s = mxy1*fz*sq(nz1)
            endif
c
            aa(n1,n2) = s
         enddo  ! n1
      enddo  ! n2
c
      if (lpr) then
         write(l6,'(/,a,i3)') 'Block ',ib
         call aprint(1,3,1,ng,ng,nf,aa,tt(i0g+1),tt(i0f+1),'Sigma * P')
      endif
c
      return
c-end-SIGP
      end

c=================== BEGIN SPLINE LIB =========================================
c
c****************************************************************
c
c    To calculate the interpolation value at (xx,yy)
c    before this one, you should run splin2 at first
c
c     (xx,yy) :  the gaussian mesh points
c     m,n     :  dimension of r_perp and z
c     x,y,z,zp:  input function, zp from subroutine splin2
c     ef      :  the output function
c
c****************************************************************
c
      subroutine splint2(xx,yy,m,n,x,y,z,zp,ef)
c
      implicit real*8(a-h,o-z)
	  include "dirhb.par"
c
	  dimension x(0:m),y(0:n),z(0:m,0:n),zp(0:m,0:n,3)
c
c---- find the location of (xx,yy) ------------------------
c
      ilo = 0
      ihi = m
1     if ((ihi-ilo).gt.1) then
         k = (ihi+ilo)/2
         if (x(k).gt.xx) then
            ihi = k
         else
            ilo = k
         endif
         goto 1
      endif
c
      jlo = 0
      jhi = n
2     if ((jhi-jlo).gt.1) then
         k = (jhi+jlo)/2
         if (y(k).gt.yy) then
            jhi = k
         else
            jlo = k
         endif
         goto 2
      endif
c
c---- cal. the interpolation value at (xx,yy) -------------------
c
      vmlo   = fu(yy,y(jlo),y(jhi),z(ilo,jlo),z(ilo,jhi),
     *            zp(ilo,jlo,2),zp(ilo,jhi,2))
      vmhi   = fu(yy,y(jlo),y(jhi),z(ihi,jlo),z(ihi,jhi),
     *            zp(ihi,jlo,2),zp(ihi,jhi,2))
      vm2xlo = fu(yy,y(jlo),y(jhi),zp(ilo,jlo,1),zp(ilo,jhi,1),
     *            zp(ilo,jlo,3),zp(ilo,jhi,3))
      vm2xhi = fu(yy,y(jlo),y(jhi),zp(ihi,jlo,1),zp(ihi,jhi,1),
     *            zp(ihi,jlo,3),zp(ihi,jhi,3))
c
      ef = fu(xx,x(ilo),x(ihi),vmlo,vmhi,vm2xlo,vm2xhi)
c
      return
      end
c****************************************************************
c
c    To do the 2D spline interpolation
c
c     ng :  the index for different parameters
c     x,y:  data points
c     m,n:  dimension of \beta and \gamma
c     z  :  the 2D function (the parameters, like V_{coll} ...
c     zp :  The derivative of z
c
c****************************************************************
c
      subroutine splin2(x,y,m,n,z,zp)
c
      implicit real*8(a-h,o-z)
	  include "dirhb.par"
c
	  dimension x(0:m),y(0:n),z(0:m,0:n),zp(0:m,0:n,3)
      dimension temp1(0:m),temp2(0:n)
c
c------- do spline for r_perp direction -------------------
c
         do 30 j=0, n
c
            do 10 i=0, m
               temp1(i)=z(i,j)
   10       continue
c
            call derlim2(x,temp1,m,derivp,derivk)
            call spline(x,temp1,temp2,m,derivp,derivk)
c
            do 20 i=0, m
               zp(i,j,1)=temp2(i)
   20       continue
c
   30    continue
c-------------------------------------
c
c------- do spline for z direction ------------------
c
         do 60 i=0, m
c
            do 40 j=0, n
               temp1(j)=z(i,j)
   40       continue
c
            call derlim2(y,temp1,n,derivp,derivk)
            call spline(y,temp1,temp2,n,derivp,derivk)
c
            do 50 j=0, n
               zp(i,j,2)=temp2(j)
   50       continue
c
   60    continue
c-------------------------------------
c
c------- do spline of d^2f/dg^2 for \gamma direction ----
c
         do 90 i=1, m
c
            do 70 j=1, n
               temp1(j)=zp(i,j,1)
   70       continue
c
            call derlim2(y,temp1,n,derivp,derivk)
            call spline(y,temp1,temp2,n,derivp,derivk)
c
            do 80 j=1, n
               zp(i,j,3)=temp2(j)
   80       continue
c
   90    continue
      return
      end
c
c
c***********************************************************
c
c    To cal. the first-order derivative
c
c    derivp: df/dx at x(1)
c    derivk: df/dx at x(N)
c
c***********************************************************
c
      subroutine derlim2(x,tem,m,derivp,derivk)
c
      implicit real*8(a-h,o-z)
c
	  dimension  x(0:m),tem(0:m)
c
      t0= (tem(1)-tem(0))/(x(1)-x(0))
      t1=(tem(2)-tem(1))/(x(2)-x(1))
      t2=(t1-t0)/(x(2)-x(0))
      derivp=2*t2*x(0) + (t0-t2*(x(0)+x(1)))
c
      t1= (tem(m-2)-tem(m-1))/(x(m-2)-x(m-1))
      t2=(tem(m-1)-tem(m))/(x(m-1)-x(m))
      t3=(t2-t1)/(x(m)-x(m-2))
      derivk=2*t3*x(m) + (t1-t3*(x(m-1)+x(m-2)))
c
      return
      end
c
c
c======================================================================c

      subroutine spline(x,y,y2,n,yp0,ypn)

c======================================================================c
c
c     SPLINE-Routine of "Numerical Recipies" p.88
c
c input:
c     X,Y       tabulated function (0..N)
c     YP0,YPN   first derivatives at point 0 and N
c               of the interpolating spline-function
c               (if larger then 1.e30, natural spline: y''= 0)
c output:
c     Y2        second derivatives of the interpolating function
c               is used as input for function SPLINT
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      dimension x(0:n),y(0:n),y2(0:n),u(0:nmax)
c
      if (nmax.lt.n) stop ' in SPLINE: nmax too small'
      if (yp0.gt.0.999d30) then
         y2(0) = 0.0
         u(0)  = 0.0
      else
         y2(0) = -0.5d0
         u(0)  = (3.d0/(x(1)-x(0))) * ((y(1)-y(0))/(x(1)-x(0))-yp0)
      endif
      do 11 i = 1,n-1
         sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p     = sig*y2(i-1) + 2.d0
         y2(i) = (sig - 1.d0)/p
         u(i)  = (6.d0*( (y(i+1)-y(i))/(x(i+1)-x(i)) -
     &                   (y(i)-y(i-1))/(x(i)-x(i-1)) )/
     &                   (x(i+1)-x(i-1)) - sig*u(i-1))/p
   11    continue
      if (ypn.gt..999d30) then
         qn = 0.0
         un = 0.0
      else
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k = n-1,0,-1
         y2(k) = y2(k)*y2(k+1)+u(k)
   12 continue
c
      return
c-end-SPLINE
      end
c======================================================================c

      subroutine splint(is,xa,ya,y2a,n,x,y,y1,y2)

c======================================================================c
c
c     SPLINT-Routine of "Numerical Recipies" p.89
c
c input:
c     XA,YA     tabulated function (0:N)
c     Y2A       first derivatives (output von SPLINE)
c     X         given value on the abscissa
c
c output:
c   is = 0:  Y  interpolated value  Y(x)
c   is = 1;  y1 in addition interpolated value of the derivativ dY/dx
c   is = 2;  y2 in addition interpolated value of 2nd derivativ d2Y/dx2
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      dimension xa(0:n),ya(0:n),y2a(0:n)
      data sixth/0.1666666666666666667d0/
c
      klo = 0
      khi = n
    1 if (khi-klo.gt.1) then
         k = (khi+klo)/2
         if (xa(k).gt.x) then
            khi = k
         else
            klo = k
         endif
         goto 1
      endif
      h = xa(khi)-xa(klo)
      if (h.eq.0.0) stop ' in SPINT: bad xa input '
      hi = 1.d0/h
      a  = (xa(khi)-x)*hi
      b  = (x-xa(klo))*hi
c
c     value of the function
      y = a*ya(klo)+b*ya(khi)+
     &    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)*sixth
c
c     first derivative
      if (is.lt.1) return
      y1 = hi*(-ya(klo)+ya(khi)) +
     &    (-(3*a**2-1)*y2a(klo)+(3*b**2-1)*y2a(khi))*h*sixth
c
c     second derivative
      if (is.lt.2) return
      y2 = a*y2a(klo) + b*y2a(khi)
c
      return
c-end-SPLINT
      end
c*************************************************************
c
c    To calculate the value and first derivative of
c    of interpolation at v
c
c*************************************************************
c
      double precision function fu(v,alo,ahi,blo,bhi,slo,shi)
      double precision v,alo,ahi,blo,bhi,slo,shi
      double precision h,a,b,sk1,sk2
      h = ahi-alo
      if (h.eq.0.) then
       write(*,*) 'pause in routine fu'
       write(*,*) ' ... bad XA input'
      endif
      a = (ahi-v)/h
      b = (v-alo)/h
      sk1=a*blo+b*bhi
      sk2=((a*a*a-a)*slo+(b*b*b-b)*shi)*(h*h)/6.0
      fu=sk1+sk2
      return
      end
c
c-----------------------------------------------------------
c
      double precision function difu(v,alo,ahi,blo,bhi,slo,shi)
      double precision v,alo,ahi,blo,bhi,slo,shi
      double precision h,a,b,sk1,sk2
      h = ahi-alo
      if (h.eq.0.) then
       write(*,*) 'pause in routine difu'
       write(*,*) ' ... bad XA input'
      endif
      a = (ahi-v)/h
      b = (v-alo)/h
      sk1=(bhi-blo)/h
      sk2=(-(3*a*a-1)*slo+(3*b*b-1)*shi)*h/6.0
      difu=sk1+sk2
      return
      end

c======================================================================c

      subroutine start(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
      logical lpr

      common /initia/ vin,rin,ain,inin,inink

      if(inin.eq.1) call startpot(lpr)
      if(inink.eq.1) call startdel(lpr)

      return
      end

c======================================================================c

      subroutine startpot(lpr)

c======================================================================c
c
c     initializes potentials
c     inin = 0:   reads fields from tape lwin
c            1:   calculates field in saxon-woods form
c
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
c
      character nucnam*2                                 ! common nucnuc
c
      dimension vso(2),r0v(2),av(2),rso(2),aso(2)
      dimension rrv(2),rls(2),vp(2),vls(2)
c
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /coulmb/ cou(ngauss),drvp(ngauss)
      common /fields/ phi(ngauss,4)
      common /gaussb/ xb(0:ngh),yb(0:ngh),zb(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /initia/ vin,rin,ain,inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(ngauss,2),vms(ngauss,2),vpstot(ngauss,2)
      common /start0/ betai,gammai
      common /tapes / l6,lin,lou,lwin,lwou,lplo
c
c=======================================================================
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
      data v0/-71.28/,akv/0.4616/
      data r0v/1.2334,1.2496/,av/0.615,0.6124/
      data vso/11.1175,8.9698/
      data rso/1.1443,1.1401/,aso/0.6476,0.6469/
c=======================================================================
c
      if (lpr) then
      write(l6,*) '****** BEGIN START *********************************'
      endif
c

c---- saxon-woods potential
      if (lpr) then
         write(l6,'(/,a)') ' Initial potentials of Saxon-Woods shape '
         write(l6,100) ' v0     = ',v0
         write(l6,100) ' kappa  = ',akv
         write(l6,100) ' lambda = ',vso
         write(l6,100) ' r0     = ',r0v
         write(l6,100) ' a      = ',av
         write(l6,100) ' r0-so  = ',rso
         write(l6,100) ' a-so   = ',aso
         write(l6,100) ' betai  = ',betai
         write(l6,100) ' gammi  = ',gammai
      endif
c
      betas = betai*half*dsqrt(5.d0/(4.d0*pi))
c
      fac = one + betas*cos(gammai*pi/180.d0)
      fac = fac*(one + betas*cos((gammai+120.d0)*pi/180.d0))
      fac = fac*(one + betas*cos((gammai-120.d0)*pi/180.d0))
      fac = fac**(-third)
c
      do ix = 0,ngh
         xx = xb(ix)**2
      do iy = 0,ngh
         yy = yb(iy)**2
      do iz = 0,ngh
         zz = zb(iz)**2
         i  = 1+ix + iy*(ngh+1) + iz*(ngh+1)*(ngh+1)
         rr = xx + yy + zz
         r  = sqrt(rr)
c
c------- Woods Saxon
         ctet = zz/rr
         cphi = (xx-yy)/rr
         p20  = 3*ctet - one
         p22  = sq(3)*cphi
         p2   = p20*cos(gammai*pi/180)+p22*sin(gammai*pi/180)
         facb = fac*(one + betas*p2)
c
         do it = 1,itx
            ita = 3-it
            rrv(it) = r0v(it)*amas**third
            rls(it) = rso(it)*amas**third
            vp(it)  = v0*(one - akv*(npr(it)-npr(ita))/amas)
            vls(it) = vp(it) * vso(it)
            u  =  vp(it) /(one + exp( (r - rrv(it)*facb) / av(it) ))
            w  = -vls(it)/(one + exp( (r - rls(it)*facb) / aso(it)))
            vps(i,it) = u
            vms(i,it) = w
c
         enddo   ! it
         if (itx.eq.1) then
            vms(i,2) = vms(i,1)
            vps(i,2) = vps(i,1)
         endif   ! itx=1
c
c------- Coulomb
         cou(i) = zero
         if (icou.ne.0) then
            rc = rrv(2)
            if (r.lt.rc) then
               c = half*(3/rc-r*r/(rc**3))
            else
               c = one/r
            endif
            cou(i)   = c*npr(2)/alphi
            vps(i,2) = vps(i,2) + hbc*cou(i)
            vms(i,2) = vms(i,2) + hbc*cou(i)
         endif
c
c------- meson fields
         do imes = 1,4
            phi(i,imes) = zero
         enddo   ! imes
         vpstot(i,1) = vps(i,1)
         vpstot(i,2) = vps(i,2)
      enddo   ! iz
      enddo   ! iy
      enddo   ! ix

c
  100 format(a,2f10.4)
c
      if (lpr) then
      write(l6,*) '****** END START ***********************************'
      endif
c
      return
c-end START
      end

c======================================================================c
      subroutine startdel(lpr)
c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character*10 inputf,outputf,outputw,fdenful,fpot,outdel   ! ttttt1
      character*15 fdenpro,fdenneu                              ! ttttt2
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
      common /blokap/ nb,mu(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /initia/ vin,rin,ain,inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo


      if (lpr)
     &write(l6,*) ' ****** BEGIN STARTDEL ****************************'
c
c
c==== initialize or read pairing potential
      do ib = 1,nb2x
         do i = 1,nhhx
	        de(i,ib) = zero
	     enddo   ! i
      enddo   ! ib
c
c---- initial constant pairing field
      do it = 1,itx
         if(del(it).lt.1.d-5) del(it)=1.d-5
         do ib = 1,nb
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
            m  = ib + (it-1)*nbx
            do n = 1,nf
               de(n+(n-1)*nh,m) = del(it)
            enddo   ! n
         enddo   ! ib
      enddo   ! it
      if (lpr) write(l6,*) ' Pairing field has been calculated'
c
      if (lpr) then
         write(l6,*) ' ****** END STARTDEL ****************************'
      endif
      return
c-end-STARTDEL
      end

