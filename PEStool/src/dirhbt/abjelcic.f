! barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
!
!  From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)
c======================================================================c

      PROGRAM DIRHBT

c======================================================================c
      implicit real*8(a-h,o-z)

      call tictoc('tic',time);

c---- sets data
      call default()

c---- reads in data
      call reader(.true.);
#ifndef ORIGINAL
      call print_header();
#endif
c---- force-parameters
      call forces(.true.);

c---- Gauss-Hermite mesh points
      call gaush(.false.);

c---- oscillator basis for single particle states
      call base(.false.);

c---- preparations
      call prep(.true.);

c---- initialization of the potentials
      call inout(1,.false.);
      call dinout(1,.false.);
      call start(.false.);

c---- wave functions
      call gaupol(.false.);

c---- single-particle matrix elements
      call singf(.false.);

c---- pairing matrix elements
      call singd(.false.);

c---- meson propagators
      call greemes(.false.);

c---- Coulomb initialization
      call greecou(.false.);
#ifndef ORIGINAL
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(initialization): ',time,'s';
      write(6,*);
#endif
c---- iteration
      call tictoc('tic',time);
      call iter(.true.);
#ifndef ORIGINAL
      write(6,*);
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(iterations    ): ',time,'s';
#endif
c---- transformation to the canonical basis
      call tictoc('tic',time);
#ifdef ORIGINAL
      call canon(.true.);
#else
      call canon_abjelcic(.false.);
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(canon transf. ): ',time,'s';
#endif
c---- center-of-mass correction
      call tictoc('tic',time);
      call centmas();
#ifndef ORIGINAL
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(centmass corr.): ',time,'s';
#endif
c---- results
      !no need to waste memory when calculating large-scale PES...
      !call resu(.true.);
      !call plot(.false.);
      call expect(2,.true.);

c---- punching of potentials to tape
      !no need to waste memory when calculating large-scale PES...
      !call inout(2,.false.);
      !call dinout(2,.false.);
      
      call tictoc('tic',time);
      call inertia_abjelcic( .true. );
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(inertia       ): ',time,'s';
      
          
      end


#ifdef ORIGINAL
c======================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the triaxial Dirac-Hartree-Bogoliubov code
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      logical lpr,lprx
      character*2 nucnam
      character*14 text3
      character*27 text1,text2
c
      common /erwar / ea,rms,betg,gamg
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /initia/ vin,rin,ain,inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /pair  / del(2),spk(2),spk0(2)

      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
      if (lpr) then
      write(l6,*) '****** BEGIN ITER **********************************'
      endif
c
      ii=0
      call gamma()
      call broyden(.false.)
      do ite = 1,maxi
         ii = ite
c
         if (lpr) then
            write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
            if (l6.ne.6)
     &      write(6,102)  ii,'.It. si=',si,'  E/A=',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
 102        format(i3,a,f10.6,2(a,f7.3),a,f6.4,a,f5.1,a,f5.2)
         endif
c
c------- loop over neutrons and protons
         do it = 1,itx
c---------- diagonalization of the Dirac-equation
            call dirhb(it,.false.)
            call denssh(it,.false.)
         enddo   ! it
c
c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- density dependence new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb field
         call coulom(.false.)
c
c------- calculation of expectation values
         lprx = .false.
         call expect(1,lprx)
c
c------- potentials in oscillator space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         call gamma()
         do it = 1,itx
	        call delta(it,.false.)
c            si = max(si,abs(spk(it)-spk0(it)))
            spk0(it) = spk(it)
            if(abs(spk(it)).lt.0.0001) del(it)=0.d0
         enddo
         call broyden(.false.)
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 20
            if (ic.eq.2) goto 30
         endif
      enddo   ! ite
   20 write(6,100) nucnam,nmas,text1,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text3,si
      goto 40
c
   30 if (lpr) write(6,101) nucnam,nmas,text2,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ii,text3,si
c
  100 format(1x,68(1h*),/,1x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(a2,i4,a27,i4,a14,f17.10)
   40 if (lpr) then
      write(l6,*) '****** END ITER ************************************'
      if (l6.ne.6)
     & write(6,*) '****** END ITER ************************************'
      endif

c
      return
c-end-ITER
      end
#else
c======================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the triaxial Dirac-Hartree-Bogoliubov code
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      logical lpr,lprx
      character*2 nucnam
      character*14 text3
      character*27 text1,text2
c
      common /erwar / ea,rms,betg,gamg
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /initia/ vin,rin,ain,inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /pair  / del(2),spk(2),spk0(2)

      logical lambdaguess;

      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
      if (lpr) then
      write(l6,*) '****** BEGIN ITER **********************************'
      endif
c


      ii=0
      call gamma()
      call broyden(.false.)

      epsi = 5.d-3;

      si   = 1.d+20;
      do ite = 1,maxi
         ii = ite
c
         if (lpr) then
            write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
            if (l6.ne.6)
     &      write(6,102)  ii,'.It. si=',si,'  E/A=',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
         endif
c


         lambdaguess = .false.;
         if( ite .le. 4 ) lambdaguess = .true.;

         ZNtol = min( 1.d-3 , 1.d-3*si );

         do it = 1,itx
            call dirhbdenssh_abjelcic(it, ZNtol, lambdaguess, .false.);
         enddo



c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- density dependence new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb field
         call coulom(.false.)
c
c------- calculation of expectation values
         lprx = .false.
         call expect(1,lprx)
c
c------- potentials in oscillator space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         call gamma()
         do it = 1,itx
	        call delta(it,.false.)
c            si = max(si,abs(spk(it)-spk0(it)))
            spk0(it) = spk(it)
            if(abs(spk(it)).lt.0.0001) del(it)=0.d0
         enddo
         call broyden(.false.)
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 20
            if (ic.eq.2) goto 30
         endif
      enddo   ! ite
   20 write(6,100) nucnam,nmas,text1,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text3,si
      goto 40
c
   30 write(6,*)'First stage completed';
      ite1 = ite;






      ii=0
      call gamma()
      call broyden(.false.)

      epsi = 1.d-5;

      si   = 1.d+20;
      do ite = 1,maxi
         ii = ite
c
         if (lpr) then
            write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
            if (l6.ne.6)
     &      write(6,102)  ii,'.It. si=',si,'  E/A=',ea,
     &             ' R=',rms,' bet=',betg,' gam=',gamg,' mix=',xmix
         endif
c

         ZNtol = min( 1.d-3 , 1.d-3*si );
         do it = 1,itx
            call dirhbdensshfull_abjelcic(it, ZNtol, .false.);
         enddo



c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- density dependence new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb field
         call coulom(.false.)
c
c------- calculation of expectation values
         lprx = .false.
         call expect(1,lprx)
c
c------- potentials in oscillator space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         call gamma()
         do it = 1,itx
	        call delta(it,.false.)
c            si = max(si,abs(spk(it)-spk0(it)))
            spk0(it) = spk(it)
            if(abs(spk(it)).lt.0.0001) del(it)=0.d0
         enddo
         call broyden(.false.)
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 21
            if (ic.eq.2) goto 31
         endif
      enddo   ! ite
   21 write(6,100) nucnam,nmas,text1,ii+ite1,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii+ite1,text3,si
      goto 40
c
   31 if (lpr) write(6,101) nucnam,nmas,text2,ii+ite1,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ii+ite1,text3,si







  100 format(1x,68(1h*),/,1x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(a2,i4,a27,i4,a14,f17.10)
  102 format(i3,a,f10.6,2(a,f7.3),a,f6.4,a,f5.1,a,f5.2)
   40 if (lpr) then
      write(l6,*) '****** END ITER ************************************'
      if (l6.ne.6)
     & write(6,*) '****** END ITER ************************************'
      endif

c
      return
c-end-ITER
      end
#endif

c======================================================================c

      subroutine centmas

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      call cmcd
#ifdef ORIGINAL
      call cmcn
#else
      call cmcn_abjelcic
#endif

      return
      end

c======================================================================c

      subroutine assert( statement , error_message )

c======================================================================c

      IMPLICIT NONE;
      LOGICAL statement;
      CHARACTER( LEN = * ) error_message;

      if( statement .eqv. .false. ) then
          write(6,'(3a)') 'Error: ', error_message, '!';
          stop;
      endif

      return;
      end;

c======================================================================c

      subroutine tictoc( tic_or_toc , time )

c======================================================================c

      IMPLICIT NONE;
      CHARACTER( LEN = * ) tic_or_toc;
      DOUBLE PRECISION time;
      INTEGER it1, it2, iclock_rate, iclock_max;

      SAVE it1;

      if( tic_or_toc .eq. 'tic' ) then
          call system_clock( it1, iclock_rate, iclock_max );
      endif
      if( tic_or_toc .eq. 'toc' ) then
          call system_clock( it2, iclock_rate, iclock_max );
          time = DBLE(it2-it1) / DBLE(iclock_rate);
      endif

      return;
      end;

c======================================================================c

      subroutine print_header()

c======================================================================c

      IMPLICIT NONE;

      write(6,*);
      write(6,*)'*****************************************************';
      write(6,*)'* Speedup by A.Bjelcic, T.Niksic and Z.Drmac        *';
      write(6,*)'*                                                   *';
      write(6,*)'* This is the original DIRHBT code (2014) directly  *';
      write(6,*)'* from CPC program library with few modifications:  *';
      write(6,*)'* 1.) The routine iter    is deleted from dirhbt.f  *';
      write(6,*)'* 2.) The routine centmas is deleted from dirhbt.f  *';
      write(6,*)'* 3.) Output of the routine canon is suppressed     *';
      write(6,*)'* 4.) Main is deleted from dirhbt.f                 *';
      write(6,*)'* 5.) Ref. BLAS routines are deleted from dirhbt.f  *';
      write(6,*)'* 6.) Some minor bugs in dirhbt.f are fixed         *';
      write(6,*)'* 7.) File abjelcic.f is added                      *';
      write(6,*)'* 8.) Makefile is slightly modified                 *';
      write(6,*)'* 9.) (Open)BLAS & LAPACK are required              *';
      write(6,*)'*                                                   *';
      write(6,*)'* If you notice something weird post an issue on    *';
      write(6,*)'* GitHub repo: github.com/abjelcic/DIRHBspeedup.    *';
      write(6,*)'*                                                   *';
      write(6,*)'*****************************************************';
      write(6,*);

      return;
      end;









c======================================================================c

      subroutine dirhbdenssh_abjelcic( it , ZNtol , lambdaguess , lpr )

c======================================================================c
      implicit real*8(a-h,o-z)
      include 'dirhb.par'

      integer it;
      double precision ZNtol;
      logical lambdaguess;
      logical lpr;

      character tk*8
      character tb*5
      common /hamham/ hh(nhhx,nb2x)
      common /deldel/ de(nhhx,nb2x)
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /fermi / ala(2),tz(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
      common /pair  / del(2),spk(2),spk0(2)

      integer MAXL;
      parameter( MAXL = 150 );

      integer NBLOCK;
      parameter( NBLOCK = ( (N0FX+2) * (N0FX+3) * (N0FX+4) ) / 6 );
      double precision h         ( NBLOCK , NBLOCK , NBX );
      double precision diagh     ( NBLOCK ,          NBX );
      double precision Delta     ( NBLOCK , NBLOCK , NBX );
      double precision diagDelta ( NBLOCK ,          NBX );
      double precision tempmat   ( NBLOCK , NBLOCK       );
      double precision h2plDelta2( NBLOCK , NBLOCK       );
      double precision asymhDelta( NBLOCK , NBLOCK       );
      double complex   hDhD      ( NBLOCK , NBLOCK , NBX );
      double complex   hDhDlambda( NBLOCK , NBLOCK       );
      double precision Q1        ( NBLOCK , NBLOCK       );
      double precision Q2        ( NBLOCK , NBLOCK       );
      double precision hQ1       ( NBLOCK , NBLOCK       );
      double precision hQ2       ( NBLOCK , NBLOCK       );
      double precision DeltaQ1   ( NBLOCK , NBLOCK       );
      double precision DeltaQ2   ( NBLOCK , NBLOCK       );
      double complex   UVj       ( NBLOCK , NBLOCK , NBX );
      double precision U         ( NBLOCK , NBLOCK       );
      double precision V         ( NBLOCK , NBLOCK       );
      double precision qfU       ( NBLOCK , NBLOCK       );
      double precision qfV       ( NBLOCK , NBLOCK       );
      double precision rho       ( NBLOCK , NBLOCK       );
      double precision kappa     ( NBLOCK , NBLOCK       );

      double precision DEGENTOL;
      parameter( DEGENTOL = 1.e-4 );
      integer DEGENBLOCK;
      parameter( DEGENBLOCK = 30 );
      double precision eij (   DEGENBLOCK ,   DEGENBLOCK );
      double precision dij (   DEGENBLOCK ,   DEGENBLOCK );
      double precision edij( 2*DEGENBLOCK , 2*DEGENBLOCK );
      double complex   CSij(   DEGENBLOCK ,   DEGENBLOCK );

      external         DLAMCH;
      double precision DLAMCH;
      integer          INFO;
      integer          NFOUND;
      double precision ABSTOL;

      integer          ISUPPZhd( 2*max(1,NBLOCK) );
      integer          LWORKhd_ASK;
      integer          LWORKhd;
      integer          LRWORKhd_ASK;
      integer          LRWORKhd;
      integer          LIWORKhd_ASK;
      integer          LIWORKhd;
      parameter(       LWORKhd  = NBLOCK*NBLOCK );
      parameter(       LRWORKhd = NBLOCK*NBLOCK );
      parameter(       LIWORKhd = NBLOCK*NBLOCK );
      double complex   WORKhd ( LWORKhd  );
      double precision RWORKhd( LRWORKhd );
      integer          IWORKhd( LIWORKhd );
      double precision W( NBLOCK , NBX );
      double complex   Q1Q2j( NBLOCK , NBLOCK );

      integer          ISUPPZh( 2*max(1,NBLOCK) );
      integer          LWORKh_ASK;
      integer          LWORKh;
      integer          LIWORKh_ASK;
      integer          LIWORKh;
      parameter(       LWORKh  = NBLOCK*NBLOCK  );
      parameter(       LIWORKh = NBLOCK*NBLOCK  );
      double precision WORKh ( LWORKh  );
      integer          IWORKh( LIWORKh );
      double precision Wh( NBLOCK          , NBX );
      double precision Zh( NBLOCK , NBLOCK , NBX );

      integer          ISUPPZed( 2*max(1,2*DEGENBLOCK) );
      integer          LWORKed_ASK;
      integer          LWORKed;
      integer          LIWORKed_ASK;
      integer          LIWORKed;
      parameter(       LWORKed  = 4*DEGENBLOCK*DEGENBLOCK  );
      parameter(       LIWORKed = 4*DEGENBLOCK*DEGENBLOCK  );
      double precision WORKed ( LWORKed  );
      integer          IWORKed( LIWORKed );
      double precision Wed( 2*DEGENBLOCK                );
      double precision Zed( 2*DEGENBLOCK , 2*DEGENBLOCK );

      double precision QuasErg( NBLOCK*NBX );
      integer ndegenmax;




      ! Preparation phase
      ierg = 0;
      do ib = 1 , nb

          nf  = id(ib,1);
          ng  = id(ib,2);
          nh  = nf + ng;
          m   = ib + (it-1)*NBX;

          call assert( nh.le.NBLOCK , 'nh > NBLOCK' );

          do j = 1 , nh
              do i = j , nh
                   h(i,j,ib) = hh( i + (j-1)*nh , m );
              enddo
          enddo
          do j = 1 , nf
              do i = j , nf
                   Delta(i,j,ib) = de( i + (j-1)*nh , m );
              enddo
          enddo

          ! Spectral decomposition of h
          ABSTOL      = DLAMCH( 'Safe minimum' );
          LWORKh_ASK  = -1;
          LIWORKh_ASK = -1;
          call dsyevr( 'V','A','L',  nh,  h(1,1,ib),NBLOCK,
     &                 +0.d+0, +1.d+10,  1,nh,  ABSTOL, NFOUND,
     &                  Wh(1,ib),  Zh(1,1,ib),NBLOCK,
     &                  ISUPPZh,
     &                  WORKh,  LWORKh_ASK,
     &                  IWORKh, LIWORKh_ASK,
     &                  INFO                                       );
          call assert(              INFO.eq.0      ,'dsyevr fail'  );
          call assert(INT( WORKh(1)+0.5).le.LWORKh ,'LWORKh small' );
          call assert(INT(IWORKh(1)+0.5).le.LIWORKh,'LIWORKh small');

          call dsyevr( 'V','A','L',  nh,  h(1,1,ib),NBLOCK,
     &                 +0.d+0, +1.d+10,  1,nh, ABSTOL, NFOUND,
     &                  Wh(1,ib),  Zh(1,1,ib),NBLOCK,
     &                  ISUPPZh,
     &                  WORKh,  LWORKh,
     &                  IWORKh, LIWORKh,
     &                  INFO                                   );
          call assert(   INFO.eq.0  , 'dsyevr fail'            );
          call assert( NFOUND.eq.nh , 'NFOUND != nh in dsyevr' );


          call dsymm( 'L','L',  nf,nf,
     &                +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &                               Zh(1,ng+1,ib),NBLOCK,
     &                +0.0d+0,        tempmat(1,1),NBLOCK  );
          call dgemm( 'T','N',  nf,nf,nf,
     &                +1.0d+0,       Zh(1,ng+1,ib),NBLOCK,
     &                                tempmat(1,1),NBLOCK,
     &                +0.0d+0,       Delta(1,1,ib),NBLOCK  );

          do j = 1 , nf
                  diagh(j,ib) = Wh(ng+j,ib);
              diagDelta(j,ib) = Delta(j,j,ib);
          enddo


          ! h^2+Delta^2
          call dsyrk( 'L','N',  nf,nf,
     &               +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &               +0.0d+0,  h2plDelta2(1,1   ),NBLOCK  );
          do i = 1 , nf
              h2plDelta2(i,i) = h2plDelta2(i,i) + diagh(i,ib)**2.d+0;
          enddo
          ! Delta*h-h*Delta
          do j = 1 , nf
              do i = j , nf
                  asymhDelta(i,j) = Delta(i,j,ib)
     &                             * ( diagh(j,ib) - diagh(i,ib) );
              enddo
          enddo
          ! hDhD = ( h + 1jDelta )*( h - 1jDelta )
          do j = 1 , nf
              do i = j , nf
                  hDhD(i,j,ib) = COMPLEX( h2plDelta2(i,j)  ,
     &                                    asymhDelta(i,j)  );
              enddo
          enddo


          ! QuasErg
          do i = 1 , nf
              ierg = ierg + 1;
              QuasErg(ierg) = Wh( ng+i , ib );
          enddo

      enddo
      if( lambdaguess .eqv. .true. ) then
          call dlasrt( 'I' , ierg , QuasErg , INFO );
          call assert( INFO .eq. 0 , 'Sort failed' );

          NZ = INT( tz(it)+0.5d0 );
          call assert( MOD(NZ,2) .eq. 0 , 'Odd number of particles' );
          call assert( NZ/2+1 .le. ierg , 'n0f too small'           );

          ala(it) = 0.5d0*( QuasErg(NZ/2) + QuasErg(NZ/2+1) );
      endif

      ndegenmax = -1;
      dl        = 100.d0;
      xl        = ala(it) - dl;
      xh        = ala(it) + dl;
      al        = ala(it);

      ! Lambda iterations
      do lit = 1 , MAXL

          snold = sn;
          sn    = 0.d0;

          do ib = 1 , nb

              nf = id(ib,1);
              ng = id(ib,2);
              nh = nf + ng;
              m  = ib + (it-1)*NBX;

              !========================================================!
              != HFB eigensolve method starts here ====================!
              !========================================================!

              ! hDhDlambda = (h-lambdaI+1jDelta)*(h-lambdaI-1jDelta)
              do j = 1 , nf
                  do i = j , nf
                      hDhDlambda(i,j) = hDhD(i,j,ib);
                  enddo
                  hDhDlambda(j,j) = + hDhDlambda(j,j)
     &                              + al*al
     &                              - 2.0d+0*al*diagh(j,ib);
              enddo

              ! Spectral decomposition of hDhDlambda
              ABSTOL       = DLAMCH( 'Safe minimum' );
              LWORKhd_ASK  = -1;
              LRWORKhd_ASK = -1;
              LIWORKhd_ASK = -1;
              call zheevr( 'V','A','L',  nf,  hDhDlambda(1,1),NBLOCK,
     &                     +0.d+0, +1.d+10,  1,nf,  ABSTOL, NFOUND,
     &                      W(1,ib),  Q1Q2j,NBLOCK,
     &                      ISUPPZhd,
     &                      WORKhd,  LWORKhd_ASK,
     &                      RWORKhd, LRWORKhd_ASK,
     &                      IWORKhd, LIWORKhd_ASK,
     &                      INFO                                       );
              call assert(            INFO.eq. 0      ,'zheevr fail'   );
              call assert(REAL( WORKhd(1)).le.LWORKhd ,'LWORKhd small' );
              call assert(     RWORKhd(1) .le.LRWORKhd,'LRWORKhd small');
              call assert(     IWORKhd(1) .le.LIWORKhd,'LIWORKhd small');

              call zheevr( 'V','A','L',  nf,  hDhDlambda(1,1),NBLOCK,
     &                     +0.d+0, +1.d+10,  1,nf, ABSTOL,  NFOUND,
     &                      W(1,ib),  Q1Q2j,NBLOCK,
     &                      ISUPPZhd,
     &                      WORKhd,  LWORKhd,
     &                      RWORKhd, LRWORKhd,
     &                      IWORKhd, LIWORKhd,
     &                      INFO                                    );
              call assert( INFO  .eq. 0  , 'zheevr fail'            );
              call assert( NFOUND.eq. nf , 'NFOUND != nf in zheevr' );

              do j = 1 , nf
                  do i = 1 , nf
                      Q1(i,j) = DREAL( Q1Q2j(i,j) );
                      Q2(i,j) = DIMAG( Q1Q2j(i,j) );
                  enddo
              enddo

              ! hQ1, hQ2, DeltaQ1 and DeltaQ2
              do j = 1 , nf
                  do i = 1 , nf
                      hQ1(i,j) = ( diagh(i,ib) - al ) * Q1(i,j);
                      hQ2(i,j) = ( diagh(i,ib) - al ) * Q2(i,j);
                  enddo
              enddo
              call dsymm( 'L','L', nf,nf,
     &                     +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                                  Q1(1,1   ),NBLOCK,
     &                     +0.d+0, DeltaQ1(1,1   ),NBLOCK  );
              call dsymm( 'L','L', nf,nf,
     &                     +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                                  Q2(1,1   ),NBLOCK,
     &                     +0.d+0, DeltaQ2(1,1   ),NBLOCK  );


              i = 1;
              do while( i .le. nf )

                  E2i = W(i,ib);

                  do j = i , nf
                      if( DABS(W(j,ib)-E2i).gt.DEGENTOL ) EXIT;
                  enddo
                  j = j - 1;
                  ! E(i) = E(i+1) = ... = E(j-1) = E(j)
                  ndegen = j-i+1;
                  call assert(ndegen.le.DEGENBLOCK,'DEGENBLOCK small');
                  ndegenmax = max( ndegenmax , ndegen );

                  ! eij (can be extracted more efficiently - dsyr2k)
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                DeltaQ2(1,i),NBLOCK,
     &                        +0.d+0,     eij(1,1),DEGENBLOCK  );
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          eij(k1,k2) = eij(k1,k2) + eij(k2,k1);
                      enddo
                  enddo
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                    hQ1(1,i),NBLOCK,
     &                        +1.d+0,     eij(1,1),DEGENBLOCK );
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        -1.d+0,      Q2(1,i),NBLOCK,
     &                                    hQ2(1,i),NBLOCK,
     &                        +1.d+0,     eij(1,1),DEGENBLOCK );

                  ! dij (can be extracted more efficiently - dsyr2k)
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                    hQ2(1,i),NBLOCK,
     &                        +0.d+0,     dij(1,1),DEGENBLOCK );
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          dij(k1,k2) = - ( dij(k1,k2) + dij(k2,k1) );
                      enddo
                  enddo
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                DeltaQ1(1,i),NBLOCK,
     &                        +1.d+0,     dij(1,1),DEGENBLOCK );
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        -1.d+0,      Q2(1,i),NBLOCK,
     &                                DeltaQ2(1,i),NBLOCK,
     &                        +1.d+0,     dij(1,1),DEGENBLOCK );


                  ! Spectral decomposition of [eij,dij;dij,-eij]
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          edij(k1       ,k2       ) = + eij(k1,k2);
                          edij(k1+ndegen,k2+ndegen) = - eij(k1,k2);
                          edij(k1+ndegen,k2       ) = + dij(k1,k2);
                      enddo
                  enddo
                  do k2 = 1 , ndegen
                      do k1 = 1 , k2-1
                          edij(k1+ndegen,k2) = + dij(k2,k1);
                      enddo
                  enddo

                  ABSTOL       = DLAMCH( 'Safe minimum' );
                  LWORKed_ASK  = -1;
                  LIWORKed_ASK = -1;
                  call dsyevr( 'V','A','L', 2*ndegen,
     &                          edij(1,1),2*DEGENBLOCK,
     &                          0.d0, 1.d10, 1,2*ndegen, ABSTOL, NFOUND,
     &                          Wed,  Zed(1,1),2*DEGENBLOCK,
     &                          ISUPPZed,
     &                          WORKed,  LWORKed_ASK,
     &                          IWORKed, LIWORKed_ASK,
     &                          INFO                                  );
                  call assert(      INFO.eq.0       , 'dsyevr fail'   );
                  call assert( WORKed(1).le.LWORKed , 'LWORKed small' );
                  call assert(IWORKed(1).le.LIWORKed, 'LIWORKed small');

                  call dsyevr( 'V','A','L', 2*ndegen,
     &                          edij(1,1),2*DEGENBLOCK,
     &                          0.d0, 1.d10, 1,2*ndegen, ABSTOL, NFOUND,
     &                          Wed,  Zed(1,1),2*DEGENBLOCK,
     &                          ISUPPZed,
     &                          WORKed,  LWORKed,
     &                          IWORKed, LIWORKed,
     &                          INFO                                  );
                  call assert(  INFO.eq.0       ,'dsyevr fail'        );
                  call assert(NFOUND.eq.2*ndegen,'NFOUND != 2*ndegen ');

                  do k2 = 1 , ndegen
                      call assert(Wed(k2+ndegen).gt.0.d0,'Negative E+');
                      do k1 = 1 , ndegen
                          Ck1k1 = Zed( k1        , k2+ndegen );
                          Sk1k2 = Zed( k1+ndegen , k2+ndegen );
                          CSij(k1,k2) = COMPLEX( Ck1k1 , Sk1k2 );
                      enddo
                  enddo

                  call zgemm( 'N','N', nf,ndegen,ndegen,
     &                          COMPLEX(+1.d+0,+0.d+0) ,
     &                                Q1Q2j(1,i),NBLOCK,
     &                             CSij(1,1),DEGENBLOCK,
     &                          COMPLEX(+0.d+0,+0.d+0) ,
     &                               UVj(1,i,ib),NBLOCK  );


                  i = j + 1;
              enddo

              !========================================================!
              != HFB eigensolve method ends here ======================!
              !========================================================!

              ! Calculating number of particles Z/N for given lambda
              snib = +0.0d+0;
              do k = 1 , nf
                  do n = 1 , nf
                      snib = snib + DIMAG(UVj(n,k,ib))**2.d0;
                  enddo
              enddo
              snib = +2.0d+0 * snib; ! Time reversal symmetry

              sn = sn + snib;


          enddo !ib




          dn = sn - tz(it);
          if( abs(dn) .lt. ZNtol ) EXIT;

          if( lit .gt. 1 ) dd = (sn-snold)/(al-alold);
          alold = al;

          if( dn .lt. 0.d0 ) then
              xl = al;
          else
              xh = al;
          endif

          if( lit.eq.1 ) then
              if( dabs(dn) .le. 0.1d0 ) then
                  al = al - dn;
              else
                  al = al - 0.1d0*sign(1.d0,dn);
              endif
          else
              if( dabs(dd) .lt. 1.d-20 ) dd = 1.d-20;
              al = al - dn/dd;

              if( al.lt.xl .or. al.gt.xh ) al = 0.5d0*(xl+xh);
          endif


      enddo !lit

      if( lit .ge. MAXL ) then
          stop 'Too many lambda iterations, something went wrong!';
      endif

      ala(it) = al;
#ifdef VERBOSE
      write(6,'(a,i3,a,i3)')'Number of lambda iterations: ', lit,
     &                      ', maximum degen. size: ', ndegenmax;
#endif



      ! Output phase
      spk(it) = 0.d0;
      iaka    = 0;
      klp     = 0;
      kla     = 0;
      do ib = 1 , nb

          nf = id(ib,1);
          ng = id(ib,2);
          nh = nf + ng;
          m  = ib + (it-1)*NBX;


          do j = 1 , nf
              do i = 1 , nf
                  U(i,j) = DREAL( UVj(i,j,ib) );
                  V(i,j) = DIMAG( UVj(i,j,ib) );
              enddo
          enddo
          call dgemm('N','N',nf+ng,nf,nf, +1.0d+0, Zh(1,ng+1,ib),NBLOCK,
     &                                                    U(1,1),NBLOCK,
     &                                    +0.0d+0,      qfU(1,1),NBLOCK);
          call dgemm('N','N',nf+ng,nf,nf, +1.0d+0, Zh(1,ng+1,ib),NBLOCK,
     &                                                    V(1,1),NBLOCK,
     &                                    +0.0d+0,      qfV(1,1),NBLOCK);

          call dsyrk( 'L','N', nh,nf, +1.0d+0, qfV(1,1),NBLOCK,
     &                                +0.0d+0, rho(1,1),NBLOCK  );

          call dsyr2k( 'L','N', nf,nf, +0.5d+0,   qfV(1,1),NBLOCK,
     &                                            qfU(1,1),NBLOCK,
     &                                 +0.0d+0, kappa(1,1),NBLOCK  );

          ! Filling rosh array from denssh
          do n2 = 1 , nh
              do n1 = n2 , nh
                  rosh( n1 + (n2-1)*nh , m ) = 2.d0 * rho(n1,n2);
                  rosh( n2 + (n1-1)*nh , m ) = 2.d0 * rho(n1,n2);
              enddo
          enddo

          ! Filling aka array from denssh
          do n2 = 1 , nf
              do n1 = n2 , nf
                  iaka = iaka + 1;

                  aka(iaka,it) = 2.d0 * kappa(n1,n2);
                  if( n1 .ne. n2 ) then
                      aka(iaka,it) = 2.d0 * aka(iaka,it);
                  endif

                  if( n1 .eq. n2 ) then
                      spk(it) = spk(it) + 0.5d0*aka(iaka,it);
                  endif

              enddo
          enddo

          ! Filling fguv array for positive Fermi energies
          ka(ib,it) = klp;
          do k = 1 , nf
              klp = klp + 1;

              ibk(klp,it) = ib;
              write(tk(klp,it),'(i2,a5)')k,tb(ib);

              call assert( W(k,ib).gt.0.d0 , 'E^2 <= 0' );

              equ(klp,it) = DSQRT( W(k,ib) );
              do n = 1 , nh
                  fguv(n+ 0,klp,it) = qfU(n,k);
                  fguv(n+nh,klp,it) = qfV(n,k);
              enddo

          enddo
          kd(ib,it) = klp - ka(ib,it);

          ! Filling fguv array for negative Dirac energies
          ka(ib,it+2) = kla;
          do k = 1 , ng
              kla = kla + 1;

              ibk(kla,it+2) = ib;
              write(tk(kla,it+2),'(i2,a5)')ng+1-k,tb(ib);

              equ(kla,it+2) = Wh(ng-k+1,ib)- al;
              do n = 1 , nh
                  fguv(n+ 0,kla,it+2) = Zh(n,ng-k+1,ib);
                  fguv(n+nh,kla,it+2) = +0.d+0;
              enddo

          enddo
          kd(ib,it+2) = kla - ka(ib,it+2);

      enddo
      nk(it  ) = klp;
      nk(it+2) = kla;


      end

c======================================================================c

      subroutine dirhbdensshfull_abjelcic( it , ZNtol , lpr )

c======================================================================c
      implicit real*8(a-h,o-z)
      include 'dirhb.par'

      integer it;
      double precision ZNtol;
      logical lpr;

      character tk*8
      character tb*5
      common /hamham/ hh(nhhx,nb2x)
      common /deldel/ de(nhhx,nb2x)
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /waveuv/ fguv(nhfbx,nkx,4),equ(nkx,4)
      common /blokap/ nb,mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /fermi / ala(2),tz(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvfx,2)
      common /pair  / del(2),spk(2),spk0(2)

      integer MAXL;
      parameter( MAXL = 150 );

      integer NBLOCK;
      parameter( NBLOCK = ( (N0FX+2) * (N0FX+3) * (N0FX+4) ) / 6 );
      double precision h         ( NBLOCK , NBLOCK , NBX );
      double precision Delta     ( NBLOCK , NBLOCK , NBX );
      double precision tempmat   ( NBLOCK , NBLOCK       );
      double precision h2plDelta2( NBLOCK , NBLOCK       );
      double precision asymhDelta( NBLOCK , NBLOCK       );
      double complex   hDhD      ( NBLOCK , NBLOCK , NBX );
      double complex   hDhDlambda( NBLOCK , NBLOCK       );
      double precision Q1        ( NBLOCK , NBLOCK       );
      double precision Q2        ( NBLOCK , NBLOCK       );
      double precision hQ1       ( NBLOCK , NBLOCK       );
      double precision hQ2       ( NBLOCK , NBLOCK       );
      double precision DeltaQ1   ( NBLOCK , NBLOCK       );
      double precision DeltaQ2   ( NBLOCK , NBLOCK       );
      double complex   UVj       ( NBLOCK , NBLOCK , NBX );
      double precision U         ( NBLOCK , NBLOCK       );
      double precision V         ( NBLOCK , NBLOCK       );
      double precision rho       ( NBLOCK , NBLOCK       );
      double precision kappa     ( NBLOCK , NBLOCK       );

      double precision DEGENTOL;
      parameter( DEGENTOL = 1.e-4 );
      integer DEGENBLOCK;
      parameter( DEGENBLOCK = 30 );
      double precision eij (   DEGENBLOCK ,   DEGENBLOCK );
      double precision dij (   DEGENBLOCK ,   DEGENBLOCK );
      double precision edij( 2*DEGENBLOCK , 2*DEGENBLOCK );
      double complex   CSij(   DEGENBLOCK ,   DEGENBLOCK );

      external         DLAMCH;
      double precision DLAMCH;
      integer          INFO;
      integer          NFOUND;
      double precision ABSTOL;

      integer          ISUPPZhd( 2*max(1,NBLOCK) );
      integer          LWORKhd_ASK;
      integer          LWORKhd;
      integer          LRWORKhd_ASK;
      integer          LRWORKhd;
      integer          LIWORKhd_ASK;
      integer          LIWORKhd;
      parameter(       LWORKhd  = NBLOCK*NBLOCK );
      parameter(       LRWORKhd = NBLOCK*NBLOCK );
      parameter(       LIWORKhd = NBLOCK*NBLOCK );
      double complex   WORKhd ( LWORKhd  );
      double precision RWORKhd( LRWORKhd );
      integer          IWORKhd( LIWORKhd );
      double precision W( NBLOCK , NBX );
      double complex   Q1Q2j( NBLOCK , NBLOCK );

      integer          ISUPPZed( 2*max(1,2*DEGENBLOCK) );
      integer          LWORKed_ASK;
      integer          LWORKed;
      integer          LIWORKed_ASK;
      integer          LIWORKed;
      parameter(       LWORKed  = 4*DEGENBLOCK*DEGENBLOCK  );
      parameter(       LIWORKed = 4*DEGENBLOCK*DEGENBLOCK  );
      double precision WORKed ( LWORKed  );
      integer          IWORKed( LIWORKed );
      double precision Wed( 2*DEGENBLOCK                );
      double precision Zed( 2*DEGENBLOCK , 2*DEGENBLOCK );

      integer ndegenmax;




      ! Preparation phase
      do ib = 1 , nb

          nf  = id(ib,1);
          ng  = id(ib,2);
          nh  = nf + ng;
          m   = ib + (it-1)*NBX;

          call assert( nh.le.NBLOCK , 'nh > NBLOCK' );

          do j = 1 , nh
              do i = j , nh
                       h(i,j,ib) = hh( i + (j-1)*nh , m );
                   Delta(i,j,ib) = de( i + (j-1)*nh , m );
              enddo
          enddo
          do j = 1 , nh
              do i = 1 , j-1
                       h(i,j,ib) =     h(j,i,ib);
                   Delta(i,j,ib) = Delta(j,i,ib);
              enddo
          enddo

          ! h^2+Delta^2
          call dsyrk( 'L','N',  nh,nh,
     &               +1.0d+0,           h(1,1,ib),NBLOCK,
     &               +0.0d+0,  h2plDelta2(1,1   ),NBLOCK  );
          call dsyrk( 'L','N',  nf,nf,
     &               +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &               +1.0d+0,  h2plDelta2(1,1   ),NBLOCK  );
          ! Delta*h-h*Delta
          call dsymm( 'L','L', nf,nh,
     &               +1.0d+0 , Delta(1,1,ib),NBLOCK,
     &                             h(1,1,ib),NBLOCK,
     &               +0.0d+0 ,  tempmat(1,1),NBLOCK   );
          asymhDelta = 0.d0;
          do j = 1 , nh
              do i = 1 , nf
                  asymhDelta(i,j) = tempmat(i,j);
              enddo
          enddo
          do j = 1 , nf
              do i = 1 , nh
                  asymhDelta(i,j) = asymhDelta(i,j) - tempmat(j,i);
              enddo
          enddo
          ! hDhD = ( h + 1jDelta )*( h - 1jDelta )
          do j = 1 , nh
              do i = j , nh
                  hDhD(i,j,ib) = COMPLEX( h2plDelta2(i,j)  ,
     &                                    asymhDelta(i,j)  );
              enddo
          enddo

      enddo


      ndegenmax = -1;
      dl        = 100.d0;
      xl        = ala(it) - dl;
      xh        = ala(it) + dl;
      al        = ala(it);

      ! Lambda iterations
      do lit = 1 , MAXL

          snold = sn;
          sn    = 0.d0;

          do ib = 1 , nb

              nf = id(ib,1);
              ng = id(ib,2);
              nh = nf + ng;
              m  = ib + (it-1)*NBX;

              !========================================================!
              != HFB eigensolve method starts here ====================!
              !========================================================!

              ! hDhDlambda = (h-lambdaI+1jDelta)*(h-lambdaI-1jDelta)
              do j = 1 , nh
                  do i = j , nh
                      hDhDlambda(i,j) = +      hDhD(i,j,ib)
     &                                  - 2.d0*al*h(i,j,ib);
                  enddo
                  hDhDlambda(j,j) = hDhDlambda(j,j) + al*al;
              enddo

              ! Spectral decomposition of hDhDlambda
              ABSTOL       = DLAMCH( 'Safe minimum' );
              LWORKhd_ASK  = -1;
              LRWORKhd_ASK = -1;
              LIWORKhd_ASK = -1;
              call zheevr( 'V','A','L',  nh,  hDhDlambda(1,1),NBLOCK,
     &                     +0.d+0, +1.d+10,  1,nh,  ABSTOL, NFOUND,
     &                      W(1,ib),  Q1Q2j,NBLOCK,
     &                      ISUPPZhd,
     &                      WORKhd,  LWORKhd_ASK,
     &                      RWORKhd, LRWORKhd_ASK,
     &                      IWORKhd, LIWORKhd_ASK,
     &                      INFO                                       );
              call assert(            INFO.eq. 0      ,'zheevr fail'   );
              call assert(REAL( WORKhd(1)).le.LWORKhd ,'LWORKhd small' );
              call assert(     RWORKhd(1) .le.LRWORKhd,'LRWORKhd small');
              call assert(     IWORKhd(1) .le.LIWORKhd,'LIWORKhd small');

              call zheevr( 'V','A','L',  nh,  hDhDlambda(1,1),NBLOCK,
     &                     +0.d+0, +1.d+10,  1,nh, ABSTOL,  NFOUND,
     &                      W(1,ib),  Q1Q2j,NBLOCK,
     &                      ISUPPZhd,
     &                      WORKhd,  LWORKhd,
     &                      RWORKhd, LRWORKhd,
     &                      IWORKhd, LIWORKhd,
     &                      INFO                                    );
              call assert( INFO  .eq. 0  , 'zheevr fail'            );
              call assert( NFOUND.eq. nh , 'NFOUND != nh in zheevr' );

              do j = 1 , nh
                  do i = 1 , nh
                      Q1(i,j) = DREAL( Q1Q2j(i,j) );
                      Q2(i,j) = DIMAG( Q1Q2j(i,j) );
                  enddo
              enddo

              ! hQ1, hQ2, DeltaQ1 and DeltaQ2
              do j = 1 , nh
                  do i = j , nh
                      tempmat(i,j) = h(i,j,ib);
                  enddo
                  tempmat(j,j) = tempmat(j,j) - al;
              enddo
              call dsymm( 'L','L', nh,nh,
     &                     +1.d+0, tempmat(1,1),NBLOCK,
     &                                  Q1(1,1),NBLOCK,
     &                     +0.d+0,     hQ1(1,1),NBLOCK  );
              call dsymm( 'L','L', nh,nh,
     &                     +1.d+0, tempmat(1,1),NBLOCK,
     &                                  Q2(1,1),NBLOCK,
     &                     +0.d+0,     hQ2(1,1),NBLOCK  );
              call dsymm( 'L','L', nf,nh,
     &                     +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                                  Q1(1,1   ),NBLOCK,
     &                     +0.d+0, DeltaQ1(1,1   ),NBLOCK  );
              call dsymm( 'L','L', nf,nh,
     &                     +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                                  Q2(1,1   ),NBLOCK,
     &                     +0.d+0, DeltaQ2(1,1   ),NBLOCK  );


              i = 1;
              do while( i .le. nh )

                  E2i = W(i,ib);

                  do j = i , nh
                      if( DABS(W(j,ib)-E2i).gt.DEGENTOL ) EXIT;
                  enddo
                  j = j - 1;
                  ! E(i) = E(i+1) = ... = E(j-1) = E(j)
                  ndegen = j-i+1;
                  call assert(ndegen.le.DEGENBLOCK,'DEGENBLOCK small');
                  ndegenmax = max( ndegenmax , ndegen );

                  ! eij
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                DeltaQ2(1,i),NBLOCK,
     &                        +0.d+0,     eij(1,1),DEGENBLOCK  );
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          eij(k1,k2) = eij(k1,k2) + eij(k2,k1);
                      enddo
                  enddo
                  call dsyr2k( 'L','T', ndegen,nh,
     &                         +0.5d+0,  Q1(1,i),NBLOCK,
     &                                  hQ1(1,i),NBLOCK,
     &                         +1.0d+0, eij(1,1),DEGENBLOCK );
                  call dsyr2k( 'L','T', ndegen,nh,
     &                         -0.5d+0,  Q2(1,i),NBLOCK,
     &                                  hQ2(1,i),NBLOCK,
     &                         +1.0d+0, eij(1,1),DEGENBLOCK );

                  ! dij
                  call dgemm( 'T','N', ndegen,ndegen,nh,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                    hQ2(1,i),NBLOCK,
     &                        +0.d+0,     dij(1,1),DEGENBLOCK );
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          dij(k1,k2) = - ( dij(k1,k2) + dij(k2,k1) );
                      enddo
                  enddo
                  call dsyr2k( 'L','T', ndegen,nf,
     &                         +0.5d+0,      Q1(1,i),NBLOCK,
     &                                  DeltaQ1(1,i),NBLOCK,
     &                         +1.0d+0,     dij(1,1),DEGENBLOCK );
                  call dsyr2k( 'L','T', ndegen,nf,
     &                         -0.5d+0,      Q2(1,i),NBLOCK,
     &                                  DeltaQ2(1,i),NBLOCK,
     &                         +1.0d+0,     dij(1,1),DEGENBLOCK );



                  ! Spectral decomposition of [eij,dij;dij,-eij]
                  do k2 = 1 , ndegen
                      do k1 = k2 , ndegen
                          edij(k1       ,k2       ) = + eij(k1,k2);
                          edij(k1+ndegen,k2+ndegen) = - eij(k1,k2);
                          edij(k1+ndegen,k2       ) = + dij(k1,k2);
                      enddo
                  enddo
                  do k2 = 1 , ndegen
                      do k1 = 1 , k2-1
                          edij(k1+ndegen,k2) = + dij(k2,k1);
                      enddo
                  enddo

                  ABSTOL       = DLAMCH( 'Safe minimum' );
                  LWORKed_ASK  = -1;
                  LIWORKed_ASK = -1;
                  call dsyevr( 'V','A','L', 2*ndegen,
     &                          edij(1,1),2*DEGENBLOCK,
     &                          0.d0, 1.d10, 1,2*ndegen, ABSTOL, NFOUND,
     &                          Wed,  Zed(1,1),2*DEGENBLOCK,
     &                          ISUPPZed,
     &                          WORKed,  LWORKed_ASK,
     &                          IWORKed, LIWORKed_ASK,
     &                          INFO                                  );
                  call assert(      INFO.eq.0       , 'dsyevr fail'   );
                  call assert( WORKed(1).le.LWORKed , 'LWORKed small' );
                  call assert(IWORKed(1).le.LIWORKed, 'LIWORKed small');

                  call dsyevr( 'V','A','L', 2*ndegen,
     &                          edij(1,1),2*DEGENBLOCK,
     &                          0.d0, 1.d10, 1,2*ndegen, ABSTOL, NFOUND,
     &                          Wed,  Zed(1,1),2*DEGENBLOCK,
     &                          ISUPPZed,
     &                          WORKed,  LWORKed,
     &                          IWORKed, LIWORKed,
     &                          INFO                                  );
                  call assert(  INFO.eq.0       ,'dsyevr fail'        );
                  call assert(NFOUND.eq.2*ndegen,'NFOUND != 2*ndegen ');

                  do k2 = 1 , ndegen
                      call assert(Wed(k2+ndegen).gt.0.d0,'Negative E+');
                      do k1 = 1 , ndegen
                          Ck1k1 = Zed( k1        , k2+ndegen );
                          Sk1k2 = Zed( k1+ndegen , k2+ndegen );
                          CSij(k1,k2) = COMPLEX( Ck1k1 , Sk1k2 );
                      enddo
                  enddo

                  call zgemm( 'N','N', nh,ndegen,ndegen,
     &                          COMPLEX(+1.d+0,+0.d+0) ,
     &                                Q1Q2j(1,i),NBLOCK,
     &                             CSij(1,1),DEGENBLOCK,
     &                          COMPLEX(+0.d+0,+0.d+0) ,
     &                               UVj(1,i,ib),NBLOCK  );


                  i = j + 1;
              enddo

              !========================================================!
              != HFB eigensolve method ends here ======================!
              !========================================================!

              ! Calculating number of particles Z/N for given lambda
              snib = +0.0d+0;
              do k = 1 , nf
                  do n = 1 , nh
                      snib = snib + DIMAG(UVj(n,k,ib))**2.d0;
                  enddo
              enddo
              do k = 1 , ng
                  do n = 1 , nh
                      snib = snib + DREAL(UVj(n,k+nf,ib))**2.d0;
                  enddo
              enddo

              snib = +2.0d+0 * snib; ! Time reversal symmetry

              sn = sn + snib;


          enddo !ib




          dn = sn - tz(it);
          if( abs(dn) .lt. ZNtol ) EXIT;

          if( lit .gt. 1 ) dd = (sn-snold)/(al-alold);
          alold = al;

          if( dn .lt. 0.d0 ) then
              xl = al;
          else
              xh = al;
          endif

          if( lit.eq.1 ) then
              if( dabs(dn) .le. 0.1d0 ) then
                  al = al - dn;
              else
                  al = al - 0.1d0*sign(1.d0,dn);
              endif
          else
              if( dabs(dd) .lt. 1.d-20 ) dd = 1.d-20;
              al = al - dn/dd;

              if( al.lt.xl .or. al.gt.xh ) al = 0.5d0*(xl+xh);
          endif


      enddo !lit

      if( lit .ge. MAXL ) then
          stop 'Too many lambda iterations, something went wrong!';
      endif

      ala(it) = al;
#ifdef VERBOSE
      write(6,'(a,i3,a,i3)')'Number of lambda iterations: ', lit,
     &                      ', maximum degen. size: ', ndegenmax;
#endif




      ! Output phase
      spk(it) = 0.d0;
      iaka    = 0;
      klp     = 0;
      kla     = 0;
      do ib = 1 , nb

          nf = id(ib,1);
          ng = id(ib,2);
          nh = nf + ng;
          m  = ib + (it-1)*NBX;

          do j = 1 , nf
              do i = 1 , nh
                  U(i,j) = DREAL( UVj(i,j,ib) );
                  V(i,j) = DIMAG( UVj(i,j,ib) );
              enddo
          enddo
          do j = 1 , ng
              do i = 1 , nh
                  U(i,nf+j) = - DIMAG( UVj(i,nf+j,ib) );
                  V(i,nf+j) = + DREAL( UVj(i,nf+j,ib) );
              enddo
          enddo

          call dsyrk( 'L','N', nh,nh, +1.0d+0,   V(1,1),NBLOCK,
     &                                +0.0d+0, rho(1,1),NBLOCK  );

          call dsyr2k( 'L','N', nf,nh, +0.5d+0,     V(1,1),NBLOCK,
     &                                              U(1,1),NBLOCK,
     &                                 +0.0d+0, kappa(1,1),NBLOCK  );

          ! Filling rosh array from denssh
          do n2 = 1 , nh
              do n1 = n2 , nh
                  rosh( n1 + (n2-1)*nh , m ) = 2.d0 * rho(n1,n2);
                  rosh( n2 + (n1-1)*nh , m ) = 2.d0 * rho(n1,n2);
              enddo
          enddo

          ! Filling aka array from denssh
          do n2 = 1 , nf
              do n1 = n2 , nf
                  iaka = iaka + 1;

                  aka(iaka,it) = 2.d0 * kappa(n1,n2);
                  if( n1 .ne. n2 ) then
                      aka(iaka,it) = 2.d0 * aka(iaka,it);
                  endif

                  if( n1 .eq. n2 ) then
                      spk(it) = spk(it) + 0.5d0*aka(iaka,it);
                  endif

              enddo
          enddo

          ! Filling fguv array for positive Fermi energies
          ka(ib,it) = klp;
          do k = 1 , nf
              klp = klp + 1;

              ibk(klp,it) = ib;
              write(tk(klp,it),'(i2,a5)')k,tb(ib);

              call assert( W(k,ib).gt.0.d0 , 'E^2 <= 0' );

              equ(klp,it) = DSQRT( W(k,ib) );
              do n = 1 , nh
                  fguv(n+ 0,klp,it) = DREAL( UVj(n,k,ib) );
                  fguv(n+nh,klp,it) = DIMAG( UVj(n,k,ib) );
              enddo

          enddo
          kd(ib,it) = klp - ka(ib,it);

          ! Filling fguv array for negative Dirac energies
          ka(ib,it+2) = kla;
          do k = 1 , ng
              kla = kla + 1;

              ibk(kla,it+2) = ib;
              write(tk(kla,it+2),'(i2,a5)')ng+1-k,tb(ib);

              call assert( W(nf+k,ib).gt.0.d0 , 'E^2 <= 0' );

              equ(kla,it+2) = - DSQRT( W(nf+k,ib) );
              do n = 1 , nh
                  fguv(n+ 0,kla,it+2) = - DIMAG( UVj(n,nf+k,ib) );
                  fguv(n+nh,kla,it+2) = + DREAL( UVj(n,nf+k,ib) );
              enddo

          enddo
          kd(ib,it+2) = kla - ka(ib,it+2);

      enddo
      nk(it  ) = klp;
      nk(it+2) = kla;


      end

c======================================================================c

      subroutine cmcn_abjelcic

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpx,lpxx,lpy,lpyy,lpz,lpzz;
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



      double precision V1mat ( NHX , NKX );
      double precision V2mat ( NHX , NKX );
      double precision Px2mat( NHX , NHX );
      double precision Py2mat( NHX , NHX );
      double precision Pz2mat( NHX , NHX );
      double precision Tmat(max(NHX,NKX),max(NHX,NKX));
      double precision Ansx( NKX , NKX );
      double precision Ansy( NKX , NKX );
      double precision Ansz( NKX , NKX );


      if( icm .lt. 2 ) then
         ecm(1) = -0.75d0 * hom;
         ecm(2) = -0.75d0 * hom;
         ecm(3) = -0.75d0 * hom;
         return;
      endif



      faccm = -0.5d0 / (amu*amas);
      ax1   = +1.0d0 / ( b0*bx * DSQRT(2.d0) );
      ay1   = +1.0d0 / ( b0*by * DSQRT(2.d0) );
      az1   = +1.0d0 / ( b0*bz * DSQRT(2.d0) );

      do it = 1 , 2

          ecmn(it) = 0.d0;

          do ib1 = 1 , nb
              k11 = kacan(ib1,it)+1;
              k12 = kacan(ib1,it)+kdcan(ib1,it);
              kd1 = kdcan(ib1,it);
              nf1 = id(ib1,1);
              do k1 = k11 , k12
                  jk = k1 - k11 + 1;
                  do n = 1 , id(ib1,1)+id(ib1,2)
                      V1mat(n,jk) = fgcan(n,k1,it);
                  enddo
              enddo


              do ib2 = 1 , nb
                  k21 = kacan(ib2,it)+1;
                  k22 = kacan(ib2,it)+kdcan(ib2,it);
                  kd2 = kdcan(ib2,it);
                  nf2 = id(ib2,1);
                  do k2 = k21 , k22
                      jk = k2 - k21 + 1;
                      do n = 1 , id(ib2,1)+id(ib2,2)
                          V2mat(n,jk) = fgcan(n,k2,it);
                      enddo
                  enddo


                  Ansx = 0.0d0;
                  Ansy = 0.0d0;
                  Ansz = 0.0d0;

                  do ifg = 1 , 2
                      nd1 = id(ib1,ifg);
                      nd2 = id(ib2,ifg);
                      i01 = ia(ib1,ifg);
                      i02 = ia(ib2,ifg);
                      do n2 = 1 , nd2
                          nx2  = nxyz(i02+n2,1);
                          ny2  = nxyz(i02+n2,2);
                          nz2  = nxyz(i02+n2,3);
                          Nsh2 = nx2 + ny2 + nz2;
                          do n1 = 1 , nd1
                              nx1  = nxyz(i01+n1,1);
                              ny1  = nxyz(i01+n1,2);
                              nz1  = nxyz(i01+n1,3);
                              Nsh1 = nx1 + ny1 + nz1;

                              Px2mat(n1,n2) = 0.0d0;
                              Py2mat(n1,n2) = 0.0d0;
                              Pz2mat(n1,n2) = 0.0d0;

                              if( MOD(Nsh1+Nsh2,2) .eq. 0 ) CYCLE;

                              lpx  =          nx1 .eq. nx2;
                              lpy  =          ny1 .eq. ny2;
                              lpz  =          nz1 .eq. nz2;
                              lpxx = abs(nx1-nx2) .eq.   1;
                              lpyy = abs(ny1-ny2) .eq.   1;
                              lpzz = abs(nz1-nz2) .eq.   1;

                              if( lpxx .and. lpy  .and. lpz  ) then
                              Px2mat(n1,n2) = ax1
     &                                        * (-1)**max(nx1,nx2)
     &                                        * SQRT(DBLE(max(nx1,nx2)))
     &                                        * (-1)**(ny2+1)
     &                                        * (-1)**(ifg+1);
                              endif
                              if( lpx  .and. lpyy .and. lpz  ) then
                              Py2mat(n1,n2) = ay1
     &                                        * SQRT(DBLE(max(ny1,ny2)))
     &                                        * (-1)**(ifg+1);
                              endif
                              if( lpx  .and. lpy  .and. lpzz ) then
                              Pz2mat(n1,n2) = az1
     &                                        * SQRT(DBLE(max(nz1,nz2)))
     &                                        * SIGN(1,nz2-nz1);
                              endif


                          enddo
                      enddo


                  ! Ansx
                  call dgemm( 'N' , 'N' ,
     &                        nd1, kd2, nd2,
     &                        1.d0,
     &                        Px2mat(1,1)                , NHX,
     &                        V2mat( (ifg-1)*nf2+1 , 1 ) , NHX,
     &                        0.d0,
     &                        Tmat(1,1)                  , max(NHX,NKX)
     &                       );
                  call dgemm( 'T' , 'N' ,
     &                        kd1, kd2, nd1,
     &                        1.d0,
     &                        V1mat( (ifg-1)*nf1+1 , 1 ) , NHX,
     &                        Tmat(1,1)                  , max(NHX,NKX),
     &                        1.d0,
     &                        Ansx(1,1)                  , NKX
     &                       );


                  ! Ansy
                  call dgemm( 'N' , 'N' ,
     &                        nd1, kd2, nd2,
     &                        1.d0,
     &                        Py2mat(1,1)                , NHX,
     &                        V2mat( (ifg-1)*nf2+1 , 1 ) , NHX,
     &                        0.d0,
     &                        Tmat(1,1)                  , max(NHX,NKX)
     &                       );
                  call dgemm( 'T' , 'N' ,
     &                        kd1, kd2, nd1,
     &                        1.d0,
     &                        V1mat( (ifg-1)*nf1+1 , 1 ) , NHX,
     &                        Tmat(1,1)                  , max(NHX,NKX),
     &                        1.d0,
     &                        Ansy(1,1)                  , NKX
     &                       );


                  ! Ansz
                  call dgemm( 'N' , 'N' ,
     &                        nd1, kd2, nd2,
     &                        1.d0,
     &                        Pz2mat(1,1)                , NHX,
     &                        V2mat( (ifg-1)*nf2+1 , 1 ) , NHX,
     &                        0.d0,
     &                        Tmat(1,1)                  , max(NHX,NKX)
     &                       );
                  call dgemm( 'T' , 'N' ,
     &                        kd1, kd2, nd1,
     &                        1.d0,
     &                        V1mat( (ifg-1)*nf1+1 , 1 ) , NHX,
     &                        Tmat(1,1)                  , max(NHX,NKX),
     &                        1.d0,
     &                        Ansz(1,1)                  , NKX
     &                       );

                  enddo


                  do k2 = k21 , k22
                      call assert(vvcan(k2,it).ge.0.d0,'vk^2 < 0');
                      call assert(vvcan(k2,it).le.1.d0,'vk^2 > 1');
                      vk2 = DSQRT(        vvcan(k2,it) );
                      uk2 = DSQRT( 1.d0 - vvcan(k2,it) );
                      do k1 = k11 , k12
                          call assert(vvcan(k1,it).ge.0.d0,'vk^2 < 0');
                          call assert(vvcan(k1,it).le.1.d0,'vk^2 > 1');
                          vk1 = DSQRT(        vvcan(k1,it) );
                          uk1 = DSQRT( 1.d0 - vvcan(k1,it) );

                          vfac = + vk1*vk1*vk2*vk2
     &                           + vk1*uk1*vk2*uk2;

                          ik = k1 - k11 + 1;
                          jk = k2 - k21 + 1;
                          Ans  = + Ansx(ik,jk)**2.d0
     &                           + Ansy(ik,jk)**2.d0
     &                           + Ansz(ik,jk)**2.d0;

                          ecmn(it) = ecmn(it) + 2.d0 * vfac * Ans;
                      enddo
                  enddo


              enddo
          enddo


          ecmn(it) = -hbc*faccm * ecmn(it);
          ecm(it)  = ecmd(it) + ecmn(it);
      enddo




      ecmn(3) = ecmn(1) + ecmn(2);
      ecm (3) = ecm (1) + ecm (2);

      return;
      end;

c======================================================================c

      subroutine canon_abjelcic(lpr)

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

      data ash/100.0/



      double precision Umat ( NHX , NKX );
      double precision Vmat ( NHX , NKX );
      double precision rho  ( NHX , NHX );
      double precision kappa( NHX , NHX );
      double precision Tmat1( NHX , NHX );
      double precision Tmat2( NHX , NHX );
      double precision Tmat3( NHX , NHX );
      integer ISUPPZ( 2*NHX );
      parameter( LWORKrho = NHX*NHX );
      double precision WORK(LWORKrho);
      parameter( LIWORKrho = NHX*NHX );
      integer IWORK( LIWORKrho );




      if(lpr) then
      write(l6,*) ' ****** BEGIN CANON *******************************';
      endif






      do it = 1 , 2





         if (lpr) write(l6,100) tit(it)
  100    format(' single-particle energies and gaps ',1x,
     &          'in the canonical basis: ',a,/,1x,66(1h-))
         if (lpr) write(l6,1000) ' ',' N pi ','  [ nx ny nz]',
     &                           'smax','eecan','vvcan','decan'
1000     format(a2,a6,a13,2a11,2x,2a11)






         klp = 0
         kla = 0
         do ib = 1,nb
	        mul  = mu(ib);
            nf   = id(ib,1);
            ng   = id(ib,2);
            nh   = nf + ng;
            nhfb = nh + nh;
            i0f  = ia(ib,1);
            i0g  = ia(ib,2);
	        mf   = ib + (it-1)*nbx;
            kf   = kd(ib,it  );
            kg   = kd(ib,it+2);
	        k0f  = ka(ib,it  );
	        k0g  = ka(ib,it+2);



            do k = k0f+1 , k0f+kf
                do n = 1 , nh
                    Umat( n , k-k0f ) = fguv(n   ,k,it);
                    Vmat( n , k-k0f ) = fguv(n+nh,k,it);
                enddo
            enddo
            call dsyrk( 'L','N',  nh,kf,
     &                  1.d0,  Vmat(1,1),NHX,
     &                  0.d0,   rho(1,1),NHX  );
            call dsyr2k( 'L','N', nf,kf,
     &                   0.5d0,  Vmat(1,1),NHX,
     &                           Umat(1,1),NHX,
     &                   0.0d0, kappa(1,1),NHX  );
            do k = k0g+1 , k0g+kg
                do n = 1 , nh
                    Vmat( n , k-k0g ) = fguv(n+nh,k,it+2);
                enddo
            enddo
            call dsyrk( 'L','N',  nh,kg,
     &                   ash,  Vmat(1,1),NHX,
     &                  1.d0,   rho(1,1),NHX  );

            do n2 = 1 , nh
                do n1 = 1 , nh

                    if( n1 .ge. n2 ) then
                        aa( n1 + (n2-1)*nh ) = rho(n1,n2);
                    else
                        aa( n1 + (n2-1)*nh ) = rho(n2,n1);
                    endif

                    if( n1.le.nf .and. n2.le.nf ) then
                        if( n1 .ge. n2 ) then
                            aaka( n1 + (n2-1)*nh ) = kappa(n1,n2);
                        else
                            aaka( n1 + (n2-1)*nh ) = kappa(n2,n1);
                        endif
                    else
                        aaka( n1 + (n2-1)*nh ) = 0.d0;
                    endif

                enddo
            enddo








        ABSTOL = DLAMCH( 'Safe minimum' );
        call dsyevr( 'V','A','L', nh,
     &               aa,nh,
     &               0.d0, 1.d9,
     &               1, nh,
     &               ABSTOL, NFOUND,
     &               v2,
     &               dd,nh,
     &               ISUPPZ,
     &               WORK,  -1,
     &               IWORK, -1,
     &               INFO                                             );
        call assert( INFO  .eq.0  , 'INFO   != 0 '                    );
        call assert(INT(WORK(1)+0.5d0).le.LWORKrho ,'LWORK  too small');
        call assert(          IWORK(1).le.LIWORKrho,'LIWORK too small');
        call dsyevr( 'V','A','L', nh,
     &               aa,nh,
     &               0.d0, 1.d9,
     &               1, nh,
     &               ABSTOL, NFOUND,
     &               v2,
     &               dd,nh,
     &               ISUPPZ,
     &               WORK,  LWORKrho,
     &               IWORK, LIWORKrho,
     &               INFO                                             );
        call assert( INFO  .eq.0  , 'INFO   != 0 '                    );
        call assert( NFOUND.eq.nh , 'NFOUND != nh'                    );







          eps = 1.d-10;
          call degen_abjelcic(nh,nh,v2,dd,hh(1,mf),eb,eps,aa,z);



          do k = 1 , nh
	          cmax = 0.d0;
	          do i = 1 , nh
	              if( DABS(dd( i + (k-1)*nh )) .gt. DABS(cmax) ) then
                      cmax = dd( i + (k-1)*nh );
                  endif
              enddo
	          if( cmax .lt. 0.d0 ) then
	              do i = 1,nh
		              dd( i + (k-1)*nh ) = - dd( i + (k-1)*nh );
	              enddo
	          endif
          enddo




          call dsymm( 'L','L', nh,nh,
     &                1.d0,  hh(1,mf) , nh,
     &                       dd       , nh,
     &                0.d0, Tmat1(1,1), NHX  );
          call dsymm( 'L','L', nh,nh,
     &                1.d0,  de(1,mf) , nh,
     &                       dd       , nh,
     &                0.d0, Tmat2(1,1), NHX  );
          call dsymm( 'L','L', nh,nh,
     &                1.d0,  aaka     , nh,
     &                       dd       , nh,
     &                0.d0, Tmat3(1,1), NHX  );
          do k = 1 , nh
               h(k) = ddot(nh, Tmat1(1,k),1, dd(1+(k-1)*nh),1);
               d(k) = ddot(nh, Tmat2(1,k),1, dd(1+(k-1)*nh),1);
              ak(k) = ddot(nh, Tmat3(1,k),1, dd(1+(k-1)*nh),1);
          enddo





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


      if(lpr) then
      write(l6,*) ' ****** END CANON *********************************';
      endif


      return
      end

c======================================================================c

      subroutine degen_abjelcic( na,n, ea, dd, bb, eb, eps, zz,z )

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

      !Note that eb is actually not valid in the original code...
      !because, when ndegen=1, then even the original code doesn't
      !calculate eb(.)

      implicit real*8 (a-h,o-z)

      dimension bb(na,n),dd(na,n),ea(n),eb(n)
      dimension zz(na,n),z(n)

      double precision    bbdd(na,n);
      double precision ddtbbdd(na,n);
      double precision   ddnew(na,n);

      do i = 1 , n-1
          call assert( ea(i) .le. ea(i+1) , 'eigenvalues not sorted!' );
      enddo
      do i = 1 , n
          do j = 1 , i-1
              call assert(DABS(bb(i,j)-bb(j,i)).le.1.d-10,'not symm');
          enddo
      enddo


      iev = 1;
      do while( iev .le. n )

          ev = ea(iev);
          do jev = iev , n
              if( DABS(ea(jev)-ev).gt.eps ) EXIT;
          enddo
          jev = jev - 1;
          ! ea(iev) = ea(iev+1) = ... = ea(jev-1) = ea(jev)
          ndegen = jev-iev+1;

          if( ndegen .eq. 1 ) then
            iev = iev + 1;
            CYCLE;
          endif

          call dsymm( 'L','L', n,ndegen,
     &                 1.d0,   bb(1,1)  ,na,
     &                         dd(1,iev),na,
     &                 0.d0, bbdd(1,1)  ,na     );
          call dgemm( 'T','N', ndegen,n,ndegen,
     &                 1.d0,      dd(1,iev),na,
     &                          bbdd(1,1)  ,na,
     &                 0.d0, ddtbbdd(1,1)  ,na  );

          call sdiag( na,ndegen , ddtbbdd , eb(iev) , ddtbbdd , z , 1 );

          call dgemm( 'N','N', n,ndegen,ndegen,
     &                 1.d0,      dd(1,iev) ,na,
     &                         ddtbbdd(1,1) ,na,
     &                 0.d0,   ddnew(1,iev) ,na  );
          do j = iev , jev
              do i = 1 , n
                  dd(i,j) = ddnew(i,j);
              enddo
          enddo



          iev = jev + 1;
      enddo


      return;
      end;












c======================================================================c

      subroutine inertia_abjelcic( lpr )

c======================================================================c
c                                                                      c
c     Calculates inertia (B1,B2,B3) and mass parameters (Bbb,Bbg,Bgg)  c
c                                                                      c
c----------------------------------------------------------------------c
      implicit double precision(a-h,o-z)

      include 'dirhb.par'
      logical lpr;

      character tb*5                  ! blokap
      character tt*8                  ! bloqua
      character tp*1,tis*1,tit*8,tl*1 ! textex

      common /blokap/ nb,mu(nbx),tb(nbx)
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /basdef/ beta0,gamma0,bx,by,bz
      common /baspar/ hom,hb0,b0
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nxyz(ntx,3),ns(ntx),np(ntx),tt(ntx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo
      common /erwar / ea,rms,betg,gamg
      common /fermi / ala(2),tz(2)
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4),akacan(nkx,4)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)

      ! new common blocks
      ! etrrqq block is weird, it does not appear in dirhbt code at all
      common /etrrqq/ etot,rr2,q0p,q2p
      common /jexcan/ ajxcan(nkx,2),ajycan(nkx,2),ajzcan(nkx,2)



      ! local memory
      double precision, dimension(:,:,:,:), allocatable :: qxx;
      double precision, dimension(:,:,:,:), allocatable :: qyy;
      double precision, dimension(:,:,:,:), allocatable :: qzz;
      double precision, dimension(:,:,:,:), allocatable :: qxy;
      double precision, dimension(:,:,:,:), allocatable :: qxz;
      double precision, dimension(:,:,:,:), allocatable :: qyz;

      double precision, dimension(:,:,:,:), allocatable :: ajx;
      double precision, dimension(:,:,:,:), allocatable :: ajy;
      double precision, dimension(:,:,:,:), allocatable :: ajz;

      double precision, dimension(:),       allocatable :: vt;
      double precision, dimension(:),       allocatable :: ut;
      double precision, dimension(:),       allocatable :: ek;
      double precision, dimension(:),       allocatable :: eqp;
      integer,          dimension(:),       allocatable :: kk;

      double precision, dimension(:,:),     allocatable :: qq0;
      double precision, dimension(:,:),     allocatable :: qq2;
      double precision, dimension(:,:),     allocatable :: qqx;
      double precision, dimension(:,:),     allocatable :: qqy;
      double precision, dimension(:,:),     allocatable :: qqz;
      double precision, dimension(:,:),     allocatable :: ajjx;
      double precision, dimension(:,:),     allocatable :: ajjy;
      double precision, dimension(:,:),     allocatable :: ajjz;

      double precision, dimension(:,:),     allocatable :: s1;
      double precision, dimension(:,:),     allocatable :: s2;
      double precision, dimension(:,:),     allocatable :: s3;
      double precision, dimension(:,:),     allocatable :: ss1;
      double precision, dimension(:,:),     allocatable :: ss2;
      double precision, dimension(:,:),     allocatable :: ss3;
      double precision, dimension(:,:),     allocatable :: sp1;
      double precision, dimension(:,:),     allocatable :: sn1;
      double precision, dimension(:,:),     allocatable :: sn2;

      double precision, dimension(:),       allocatable :: BB1;
      double precision, dimension(:),       allocatable :: BB2;
      double precision, dimension(:),       allocatable :: BB3;

      logical lx0,ly0,lz0, lx1,ly1,lz1, lx2,ly2,lz2;
      double precision M1(2,2);
      double precision M3(2,2);
      double precision invM1_M3_invM1(2,2);










      if(lpr) then
      write(l6,*) ' ****** BEGIN INERTIA *****************************';
      endif










      ! constant factors
      pi = 4.d0 * DATAN(1.d0);

      fac0 = DSQRT(  5.d0/(16.d0*pi) ); !Q20
      fac2 = DSQRT( 15.d0/(16.d0*pi) ); !sqrt(2)*Q22
      fact = hom**2;
      facm = 9.d0*2.0736d0 * (amas**(10.d0/3.d0)) / (16.d0*pi*pi*b0**4);

      sxx = ( 1.d0/DSQRT(2.d0) * bx )**2;
      syy = ( 1.d0/DSQRT(2.d0) * by )**2;
      szz = ( 1.d0/DSQRT(2.d0) * bz )**2;
      sx  = ( 1.d0/DSQRT(2.d0) * bx );
      sy  = ( 1.d0/DSQRT(2.d0) * by );
      sz  = ( 1.d0/DSQRT(2.d0) * bz );










      ! allocation and initialization of local memory
      allocate( BB1(3) ); BB1 = 0.d0;
      allocate( BB2(3) ); BB2 = 0.d0;
      allocate( BB3(3) ); BB3 = 0.d0;


      nalloc = 0;
      do ib = 1 , nb
          do ifg = 1 , 2
              nalloc = max( nalloc , id(ib,ifg) );
          enddo
      enddo
      allocate( qxx( nalloc , nalloc , 2 , nb ) ); qxx = 0.d0;
      allocate( qyy( nalloc , nalloc , 2 , nb ) ); qyy = 0.d0;
      allocate( qzz( nalloc , nalloc , 2 , nb ) ); qzz = 0.d0;

      allocate( qxy( nalloc , nalloc , 2 , nb ) ); qxy = 0.d0;
      allocate( qxz( nalloc , nalloc , 2 , nb ) ); qxz = 0.d0;
      allocate( qyz( nalloc , nalloc , 2 , nb ) ); qyz = 0.d0;

      allocate( ajx( nalloc , nalloc , 2 , nb ) ); ajx = 0.d0;
      allocate( ajy( nalloc , nalloc , 2 , nb ) ); ajy = 0.d0;
      allocate( ajz( nalloc , nalloc , 2 , nb ) ); ajz = 0.d0;


      kalloc = 0;
      do it = 1 , 2
          kalloc = max( kalloc , nkcan(it) );
      enddo
      allocate(   vt( kalloc          ) ); vt   = 0.d0;
      allocate(   ut( kalloc          ) ); ut   = 0.d0;
      allocate(   ek( kalloc          ) ); ek   = 0.d0;
      allocate(  eqp( kalloc          ) ); eqp  = 0.d0;
      allocate(   kk( kalloc          ) ); kk   = 0;

      allocate(  qq0( kalloc , kalloc ) ); qq0  = 0.d0;
      allocate(  qq2( kalloc , kalloc ) ); qq2  = 0.d0;
      allocate(  qqx( kalloc , kalloc ) ); qqx  = 0.d0;
      allocate(  qqy( kalloc , kalloc ) ); qqy  = 0.d0;
      allocate(  qqz( kalloc , kalloc ) ); qqz  = 0.d0;
      allocate( ajjx( kalloc , kalloc ) ); ajjx = 0.d0;
      allocate( ajjy( kalloc , kalloc ) ); ajjy = 0.d0;
      allocate( ajjz( kalloc , kalloc ) ); ajjz = 0.d0;


      nkalloc = 0;
      do it = 1 , 2
          do ib = 1 , nb
              nkalloc = max( nkalloc , kdcan(ib,it) );
          enddo
      enddo
      allocate(  s1( nkalloc , nkalloc ) ); s1  = 0.d0;
      allocate(  s2( nkalloc , nkalloc ) ); s2  = 0.d0;
      allocate(  s3( nkalloc , nkalloc ) ); s3  = 0.d0;
      allocate( ss1( nkalloc , nkalloc ) ); ss1 = 0.d0;
      allocate( ss2( nkalloc , nkalloc ) ); ss2 = 0.d0;
      allocate( ss3( nkalloc , nkalloc ) ); ss3 = 0.d0;
      allocate( sp1( nkalloc , nkalloc ) ); sp1 = 0.d0;
      allocate( sn1( nkalloc , nkalloc ) ); sn1 = 0.d0;
      allocate( sn2( nkalloc , nkalloc ) ); sn2 = 0.d0;


      aM001 = 0.d0;
      aM021 = 0.d0;
      aM221 = 0.d0;
      aM002 = 0.d0;
      aM022 = 0.d0;
      aM222 = 0.d0;
      aM003 = 0.d0;
      aM023 = 0.d0;
      aM223 = 0.d0;
      aMp12 = 0.d0;
      aMn12 = 0.d0;
      aMn22 = 0.d0;
      aMp13 = 0.d0;
      aMn13 = 0.d0;
      aMn23 = 0.d0;










      ! calculating qxx, qyy, qzz, qxy, qxz, qyz
      qxx = 0.d0;
      qyy = 0.d0;
      qzz = 0.d0;
      qxy = 0.d0;
      qxz = 0.d0;
      qyz = 0.d0;
      do ib = 1 , nb
          do ifg = 1 , 2

              do n2 = 1 , id(ib,ifg)
                  nx2  = nxyz(ia(ib,ifg)+n2,1);
                  ny2  = nxyz(ia(ib,ifg)+n2,2);
                  nz2  = nxyz(ia(ib,ifg)+n2,3);
                  Nsh2 = nx2 + ny2 + nz2;

                  do n1 = n2 , id(ib,ifg)
                      nx1  = nxyz(ia(ib,ifg)+n1,1);
                      ny1  = nxyz(ia(ib,ifg)+n1,2);
                      nz1  = nxyz(ia(ib,ifg)+n1,3);
                      Nsh1 = nx1 + ny1 + nz1;

                      if( n1 .eq. n2 ) then
                          qxx(n1,n1,ifg,ib) = (2.d0*DBLE(nx1)+1.d0)*sxx;
                          qyy(n1,n1,ifg,ib) = (2.d0*DBLE(ny1)+1.d0)*syy;
                          qzz(n1,n1,ifg,ib) = (2.d0*DBLE(nz1)+1.d0)*szz;
                          CYCLE;
                      endif

                      if( abs( Nsh1 - Nsh2 ) .gt. 2 ) CYCLE;

                      lx0 = abs(nx1-nx2) .eq. 0;
                      ly0 = abs(ny1-ny2) .eq. 0;
                      lz0 = abs(nz1-nz2) .eq. 0;
                      lx1 = abs(nx1-nx2) .eq. 1;
                      ly1 = abs(ny1-ny2) .eq. 1;
                      lz1 = abs(nz1-nz2) .eq. 1;
                      lx2 = abs(nx1-nx2) .eq. 2;
                      ly2 = abs(ny1-ny2) .eq. 2;
                      lz2 = abs(nz1-nz2) .eq. 2;

                      if( lx2 .and. ly0 .and. lz0 ) then
                          txx = + sxx;
                          txx = txx * DSQRT( DBLE( MIN(nx1,nx2)+1 ) );
                          txx = txx * DSQRT( DBLE( MAX(nx1,nx2)   ) );

                          qxx(n1,n2,ifg,ib) = txx;
                      endif

                      if( lx0 .and. ly2 .and. lz0 ) then
                          tyy = - syy;
                          tyy = tyy * DSQRT( DBLE( MIN(ny1,ny2)+1 ) );
                          tyy = tyy * DSQRT( DBLE( MAX(ny1,ny2)   ) );

                          qyy(n1,n2,ifg,ib) = tyy;
                      endif

                      if( lx0 .and. ly0 .and. lz2 ) then
                          tzz = + szz;
                          tzz = tzz * DSQRT( DBLE( MIN(nz1,nz2)+1 ) );
                          tzz = tzz * DSQRT( DBLE( MAX(nz1,nz2)   ) );

                          qzz(n1,n2,ifg,ib) = tzz;
                      endif

                      if( lx1 .and. ly1 .and. lz0 ) then
                          txy = (-1)**(nx2+ny2+1) * sx * sy;
                          txy = txy * DSQRT( DBLE( MAX(nx1,nx2) ) );
                          txy = txy * DSQRT( DBLE( MAX(ny1,ny2) ) );
                          if( ny1-ny2 .eq. -1 ) txy = -1.d0 * txy;

                          qxy(n1,n2,ifg,ib) = txy;
                      endif

                      if( lx1 .and. ly0 .and. lz1 ) then
                          txz = (-1)**(ifg-1+nx2+ny2+1) * sx * sz;
                          txz = txz * DSQRT( DBLE( MAX(nx1,nx2) ) );
                          txz = txz * DSQRT( DBLE( MAX(nz1,nz2) ) );

                          qxz(n1,n2,ifg,ib) = txz;
                      endif

                      if( lx0 .and. ly1 .and. lz1 ) then
                          tyz = (-1)**(ifg-1) * sy * sz;
                          tyz = tyz * DSQRT( DBLE( MAX(ny1,ny2) ) );
                          tyz = tyz * DSQRT( DBLE( MAX(nz1,nz2) ) );
                          if( ny1-ny2 .eq. -1 ) tyz = -1.d0 * tyz;

                          qyz(n1,n2,ifg,ib) = tyz;
                      endif

                  enddo
              enddo

              ! upper triangle
              do n2 = 1 , id(ib,ifg)
                  do n1 = 1 , n2-1
                      qxx(n1,n2,ifg,ib) = + qxx(n2,n1,ifg,ib);
                      qyy(n1,n2,ifg,ib) = + qyy(n2,n1,ifg,ib);
                      qzz(n1,n2,ifg,ib) = + qzz(n2,n1,ifg,ib);
                      qxy(n1,n2,ifg,ib) = - qxy(n2,n1,ifg,ib);
                      qxz(n1,n2,ifg,ib) = - qxz(n2,n1,ifg,ib);
                      qyz(n1,n2,ifg,ib) = - qyz(n2,n1,ifg,ib);
                  enddo
              enddo

          enddo
      enddo

      ! calculating ajx, ajy, ajz
      ajx = 0.d0;
      ajy = 0.d0;
      ajz = 0.d0;
      do ib = 1 , nb
          do ifg = 1 , 2

              do n2 = 1 , id(ib,ifg)
                  nx2  = nxyz(ia(ib,ifg)+n2,1);
                  ny2  = nxyz(ia(ib,ifg)+n2,2);
                  nz2  = nxyz(ia(ib,ifg)+n2,3);

                  do n1 = n2 , id(ib,ifg)
                      nx1  = nxyz(ia(ib,ifg)+n1,1);
                      ny1  = nxyz(ia(ib,ifg)+n1,2);
                      nz1  = nxyz(ia(ib,ifg)+n1,3);

                      if( n1 .eq. n2 ) then
                          ajx(n1,n1,ifg,ib) = 0.5d0 * (-1)**(nx1+1);
                          ajy(n1,n1,ifg,ib) = 0.5d0 * (-1)**(ny1+1);
                          ajz(n1,n1,ifg,ib) = 0.5d0 * (-1)**(nx1+ny1+1);
                          CYCLE;
                      endif

                      lx0 = abs(nx1-nx2) .eq. 0;
                      ly0 = abs(ny1-ny2) .eq. 0;
                      lz0 = abs(nz1-nz2) .eq. 0;
                      lx1 = abs(nx1-nx2) .eq. 1;
                      ly1 = abs(ny1-ny2) .eq. 1;
                      lz1 = abs(nz1-nz2) .eq. 1;

                      if( lx0 .and. ly1 .and. lz1 ) then
                          s = + 1.d0;
                          s = s * DSQRT( DBLE( MAX(ny1,ny2) ) );
                          s = s * DSQRT( DBLE( MAX(nz1,nz2) ) );

                          sum = 0.5d0 * ( by/bz + bz/by );
                          dif = 0.5d0 * ( by/bz - bz/by );
                          if( nz2.eq.nz1+1 .and. ny2.eq.ny1+1 ) then
                              s = + dif * s;
                          endif
                          if( nz2.eq.nz1-1 .and. ny2.eq.ny1-1 ) then
                              s = + dif * s;
                          endif
                          if( nz2.eq.nz1+1 .and. ny2.eq.ny1-1 ) then
                              s = - sum * s;
                          endif
                          if( nz2.eq.nz1-1 .and. ny2.eq.ny1+1 ) then
                              s = - sum * s;
                          endif


                          ajx(n1,n2,ifg,ib) = s;
                      endif

                      if( lx1 .and. ly0 .and. lz1 ) then
                          s = + 1.d0;
                          s = s * DSQRT( DBLE( MAX(nx1,nx2) ) );
                          s = s * DSQRT( DBLE( MAX(nz1,nz2) ) );

                          sum = 0.5d0 * ( bz/bx + bx/bz );
                          dif = 0.5d0 * ( bz/bx - bx/bz );
                          if( nx2.eq.nx1+1 .and. nz2.eq.nz1+1 ) then
                              s = + dif * s;
                          endif
                          if( nx2.eq.nx1-1 .and. nz2.eq.nz1-1 ) then
                              s = - dif * s;
                          endif
                          if( nx2.eq.nx1+1 .and. nz2.eq.nz1-1 ) then
                              s = + sum * s;
                          endif
                          if( nx2.eq.nx1-1 .and. nz2.eq.nz1+1 ) then
                              s = - sum * s;
                          endif


                          ajy(n1,n2,ifg,ib) = (-1)**(nx2+ny2+1) * s;
                      endif

                      if( lx1 .and. ly1 .and. lz0 ) then
                          s = + 1.d0;
                          s = s * DSQRT( DBLE( MAX(nx1,nx2) ) );
                          s = s * DSQRT( DBLE( MAX(ny1,ny2) ) );

                          sum = 0.5d0 * ( bx/by + by/bx );
                          dif = 0.5d0 * ( bx/by - by/bx );
                          if( nx2.eq.nx1+1 .and. ny2.eq.ny1+1 ) then
                              s = + dif * s;
                          endif
                          if( nx2.eq.nx1-1 .and. ny2.eq.ny1-1 ) then
                              s = + dif * s;
                          endif
                          if( nx2.eq.nx1+1 .and. ny2.eq.ny1-1 ) then
                              s = + sum * s;
                          endif
                          if( nx2.eq.nx1-1 .and. ny2.eq.ny1+1 ) then
                              s = + sum * s;
                          endif


                          ajz(n1,n2,ifg,ib) = (-1)**(nx2+ny2+1) * s;
                      endif

                  enddo
              enddo

              ! upper triangle
              do n2 = 1 , id(ib,ifg)
                  do n1 = 1 , n2-1
                      ajx(n1,n2,ifg,ib) = + ajx(n2,n1,ifg,ib);
                      ajy(n1,n2,ifg,ib) = + ajy(n2,n1,ifg,ib);
                      ajz(n1,n2,ifg,ib) = + ajz(n2,n1,ifg,ib);
                  enddo
              enddo

          enddo
      enddo










      ! calculating
      ! BB1   , BB2   , BB3   ,
      ! aM003 , aM223 , aM023 ,
      ! aM002 , aM222 , aM022 ,
      ! aM001 , aM221 , aM021 ,
      ! aMp13 , aMn13 , aMn23 ,
      ! aMp12 , aMn12 , aMn22
      do it = 1 , 2


          ! calculating non-used arrays because the original code does
          km = 0;
          do k = 1 , nkcan(it)
              lam = ala(it);

              vt(k)  = DSQRT( ABS(     vvcan(k,it)) );
              ut(k)  = DSQRT( ABS(1.d0-vvcan(k,it)) );
              ek(k)  = eecan(k,it);
              eqp(k) = DSQRT( (eecan(k,it)-lam)**2 + decan(k,it)**2 );

              km     = km + 1;
              kk(km) = k;

              call assert( k.le.NKX , 'NKX too small' );
          enddo
          do ib = 1 , nb
              do k = kacan(ib,it)+1 , kacan(ib,it)+kdcan(ib,it)
                  call assert( ib .eq. ibk(k,it) , 'ibk(k,it) wrong' );
              enddo
          enddo



          ! calculating qq0, qq2, qqx, qqy, qqz, ajjx, ajjy, ajjz
          qq0  = 0.d0;
          qq2  = 0.d0;
          qqx  = 0.d0;
          qqy  = 0.d0;
          qqz  = 0.d0;
          ajjx = 0.d0;
          ajjy = 0.d0;
          ajjz = 0.d0;
          do ib = 1 , nb
              k1 = kacan(ib,it)+1;
              k2 = kacan(ib,it)+kdcan(ib,it);
              kd = kdcan(ib,it);
              nf = id(ib,1);

              s1  = 0.d0;
              s2  = 0.d0;
              s3  = 0.d0;
              ss1 = 0.d0;
              ss2 = 0.d0;
              ss3 = 0.d0;
              sp1 = 0.d0;
              sn1 = 0.d0;
              sn2 = 0.d0;
              do ifg = 1 , 2
                  nd = id(ib,ifg);

                  ! Performs C = C + alpha * U^T*A*U transformation of (anti)symmetric A
                  !
                  ! U = fgcan( (ifg-1)*nf + 1...(ifg-1)*nf + nd , k1...k2 , it )
                  !
                  ! (symmetric A)
                  ! A = qxx( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C =  s1( 1...kd , 1...kd )
                  ! A = qyy( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C =  s2( 1...kd , 1...kd )
                  ! A = qzz( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C =  s3( 1...kd , 1...kd )
                  ! A = ajx( 1...nd , 1...nd , ifg , ib ), alpha = (-1)^(ifg+1), C = ss1( 1...kd , 1...kd )
                  ! A = ajy( 1...nd , 1...nd , ifg , ib ), alpha = (-1)^(ifg+1), C = ss2( 1...kd , 1...kd )
                  ! A = ajz( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C = ss3( 1...kd , 1...kd )
                  !
                  ! (anti-symmetric A)
                  ! A = qyz( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C = sp1( 1...kd , 1...kd )
                  ! A = qxz( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C = sn1( 1...kd , 1...kd )
                  ! A = qxy( 1...nd , 1...nd , ifg , ib ), alpha = 1.0         , C = sn2( 1...kd , 1...kd )


                  ! s1
                  call  symmUtAU( nd , kd , 1.d0                     ,
     &                            qxx(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            s1(1,1)                   , nkalloc );
                  ! s2
                  call  symmUtAU( nd , kd , 1.d0                     ,
     &                            qyy(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            s2(1,1)                   , nkalloc );
                  ! s3
                  call  symmUtAU( nd , kd , 1.d0                     ,
     &                            qzz(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            s3(1,1)                   , nkalloc );
                  ! ss1
                  call  symmUtAU( nd , kd , DBLE((-1)**(ifg+1))      ,
     &                            ajx(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            ss1(1,1)                  , nkalloc );
                  ! ss2
                  call  symmUtAU( nd , kd , DBLE((-1)**(ifg+1))      ,
     &                            ajy(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            ss2(1,1)                  , nkalloc );
                  ! ss3
                  call  symmUtAU( nd , kd , 1.d0                     ,
     &                            ajz(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            ss3(1,1)                  , nkalloc );


                  ! sp1
                  call asymmUtAU( nd , kd , 1.d0                     ,
     &                            qyz(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            sp1(1,1)                  , nkalloc );
                  ! sn1
                  call asymmUtAU( nd , kd , 1.d0                     ,
     &                            qxz(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            sn1(1,1)                  , nkalloc );
                  ! sn2
                  call asymmUtAU( nd , kd , 1.d0                     ,
     &                            qxy(1,1,ifg,ib)           , nalloc ,
     &                            fgcan((ifg-1)*nf+1,k1,it) , NHX    ,
     &                            sn2(1,1)                  , nkalloc );

              enddo



              do j = 1 , kd
                  do i = 1 , kd

                  ki = k1 + (i-1);
                  kj = k1 + (j-1);

                   qq0(ki,kj) = fac0 * ( 2.d0*s3(i,j)-s1(i,j)-s2(i,j) );
                   qq2(ki,kj) = fac2 * ( s1(i,j)-s2(i,j) );
                   qqx(ki,kj) = 2.d0 * sp1(i,j);
                   qqy(ki,kj) = 2.d0 * sn1(i,j);
                   qqz(ki,kj) = 2.d0 * sn2(i,j);
                  ajjx(ki,kj) = ss1(i,j);
                  ajjy(ki,kj) = ss2(i,j);
                  ajjz(ki,kj) = ss3(i,j);

                  ! In the original code, one encounters the following:
                  !
                  ! qqx(k1,k2) = 2.*sp1
                  ! qqx(k2,k1) = 2.*sp1
                  ! qqy(k1,k2) = 2.*sn1
                  ! qqy(k2,k1) = 2.*sn1
                  ! qqz(k1,k2) = 2.*sn2
                  ! qqz(k2,k1) = 2.*sn2
                  !
                  ! where sp1,sn1,sn2 are obtained with:
                  !
                  ! sp1 = sp1 + fgt*qyz(n1,n2,ifg,ib1)
                  ! sn1 = sn1 + fgt*qxz(n1,n2,ifg,ib1)
                  ! sn2 = sn2 + fgt*qxy(n1,n2,ifg,ib1)
                  !
                  ! Since qxy,qxz,qyz are anti-symmetric matrices,
                  ! sp1,sn1,sn2 matrices are also anti-symmetric, and
                  ! thus qqx,qqy,qqz are also anti-symmetric.
                  ! In the original code, one loops over
                  ! only the lower triangle of (k1,k2) pairs, and in
                  ! the end symmetrizes matrices qqx,qqy,qqz.
                  ! However, this presents no problem, since qqx,qqy,qqz
                  ! matrices are used only in the squared form:
                  !
                  ! aMp13 = aMp13 + tuve3*qqx(k1,k2)**2
                  ! aMn13 = aMn13 + tuve3*qqy(k1,k2)**2
                  ! aMn23 = aMn23 + tuve3*qqz(k1,k2)**2

                  ! aMp12 = aMp12 + tuve2*qqx(k1,k2)**2
                  ! aMn12 = aMn12 + tuve2*qqy(k1,k2)**2
                  ! aMn22 = aMn22 + tuve2*qqz(k1,k2)**2
                  !
                  ! and thus, I won't symmetrize matrices qqx,qqy,qqz.


                  enddo
              enddo
              call assert( k1+kd-1 .le. kalloc , 'kalloc too small' );

              do i = 1 , kd
                  ki = k1 + (i-1);
                  ajxcan( ki , it ) = ajjx( ki , ki );
                  ajycan( ki , it ) = ajjy( ki , ki );
                  ajzcan( ki , it ) = ajjz( ki , ki );
              enddo
              call assert( k1+kd-1 .le. NKX , 'NKX too small' );


          enddo




          do ib = 1 , nb
              do k2 = kacan(ib,it)+1 , kacan(ib,it)+kdcan(ib,it)
                  call assert(vvcan(k2,it).ge.0.d0,'vk^2 < 0');
                  call assert(vvcan(k2,it).le.1.d0,'vk^2 > 1');
                  vk2  = DSQRT(        vvcan(k2,it) );
                  uk2  = DSQRT( 1.d0 - vvcan(k2,it) );
                  ek2  = eecan(k2,it);
                  dk2  = decan(k2,it);
                  eqp2 = DSQRT( (ek2-ala(it))**2 + dk2**2 );

                  do k1 = kacan(ib,it)+1 , kacan(ib,it)+kdcan(ib,it)
                      call assert(vvcan(k1,it).ge.0.d0,'vk^2 < 0');
                      call assert(vvcan(k1,it).le.1.d0,'vk^2 > 1');
                      vk1  = DSQRT(        vvcan(k1,it) );
                      uk1  = DSQRT( 1.d0 - vvcan(k1,it) );
                      ek1  = eecan(k1,it);
                      dk1  = decan(k1,it);
                      eqp1 = DSQRT( (ek1-ala(it))**2 + dk1**2 );


                      denom1 = ( eqp1 + eqp2 )**1 + 1.d-8;
                      denom2 = ( eqp1 + eqp2 )**2 + 1.d-8;
                      denom3 = ( eqp1 + eqp2 )**3 + 1.d-8;


                      tuv   = 2.d0 * ( uk1*vk2 + vk1*uk2 )**2;
                      tuvc  = 2.d0 * ( uk1*vk2 - vk1*uk2 )**2 / denom1;
                      tuve1 = tuv / denom1;
                      tuve2 = tuv / denom2;
                      tuve3 = tuv / denom3;



                      BB1(it) = BB1(it) + tuvc * ajjx(k1,k2)**2;
                      BB2(it) = BB2(it) + tuvc * ajjy(k1,k2)**2;
                      BB3(it) = BB3(it) + tuvc * ajjz(k1,k2)**2;
                      ! write(*,*) it , k1 , k2 , tuvc , ajjz(k1,k2)

                      aM003 = aM003 + tuve3 * qq0(k1,k2)**2;
                      aM223 = aM223 + tuve3 * qq2(k1,k2)**2;
                      aM023 = aM023 + tuve3 * qq0(k1,k2)*qq2(k2,k1);

                      aM002 = aM002 + tuve2 * qq0(k1,k2)**2;
                      aM222 = aM222 + tuve2 * qq2(k1,k2)**2;
                      aM022 = aM022 + tuve2 * qq0(k1,k2)*qq2(k2,k1);

                      aM001 = aM001 + tuve1 * qq0(k1,k2)**2;
                      aM221 = aM221 + tuve1 * qq2(k1,k2)**2;
                      aM021 = aM021 + tuve1 * qq0(k1,k2)*qq2(k2,k1);

                      aMp13 = aMp13 + tuve3 * qqx(k1,k2)**2;
                      aMn13 = aMn13 + tuve3 * qqy(k1,k2)**2;
                      aMn23 = aMn23 + tuve3 * qqz(k1,k2)**2;

                      aMp12 = aMp12 + tuve2 * qqx(k1,k2)**2;
                      aMn12 = aMn12 + tuve2 * qqy(k1,k2)**2;
                      aMn22 = aMn22 + tuve2 * qqz(k1,k2)**2;


                  enddo
              enddo
          enddo


      ! write(*,*) it , BB1(it);

      enddo










      ! calculating B00,B02,B22
      M1(1,1) = aM001;
      M1(1,2) = aM021;
      M1(2,1) = aM021;
      M1(2,2) = aM221;

      M3(1,1) = aM003;
      M3(1,2) = aM023;
      M3(2,1) = aM023;
      M3(2,2) = aM223;
      call invB_A_invB( M3 , M1 , invM1_M3_invM1 );

      B00 = facm * invM1_M3_invM1(1,1);
      B02 = facm * invM1_M3_invM1(1,2);
      B20 = facm * invM1_M3_invM1(2,1);
      B22 = facm * invM1_M3_invM1(2,2);

      !call assert( ABS(B20-B02)/ABS(B20) .le. 1.d-10 , 'B20 != B02' );
      ! write(77,*) 'aM001=',aM001
      ! write(77,*) 'aM021=',aM021
      ! write(77,*) 'aM221=',aM221
      ! write(77,*) 'aM002=',aM002
      ! write(77,*) 'aM022=',aM022
      ! write(77,*) 'aM222=',aM222
      ! write(77,*) 'aM003=',aM003
      ! write(77,*) 'aM023=',aM023
      ! write(77,*) 'aM223=',aM223
      ! write(77,*) 'aMp13=',aMp13
      ! write(77,*) 'aMn13=',aMn13
      ! write(77,*) 'aMn23=',aMn23



      ! translate B00,B02,B22 to Bbb,Bbg,Bgg
      betax  = betg;
      gammax = gamg/180.d0 * pi;
      gammat = gammax;
      if( betax .lt. 1.d-3 ) gammat = 0.d0;

      cgcg = DCOS(gammat) * DCOS(gammat);
      sgsg = DSIN(gammat) * DSIN(gammat);
      sgcg = DSIN(gammat) * DCOS(gammat);

      Bbb = B00*(+cgcg) + B22*(+sgsg) + B02*(+2.d0*sgcg);
      Bbg = B00*(-sgcg) + B22*(+sgcg) + B02*(+cgcg-sgsg);
      Bgg = B00*(+sgsg) + B22*(+cgcg) + B02*(-2.d0*sgcg);



      ! translate BB1,BB2,BB3 to B1,B2,B3
      BB1(3) = BB1(1) + BB1(2);
      BB2(3) = BB2(1) + BB2(2);
      BB3(3) = BB3(1) + BB3(2);

      B1 = 0.d0;
      B2 = 0.d0;
      B3 = 0.d0;
      if( betax .lt. 1.d-3 ) then

         B1 = Bbb;
         B2 = Bbb;
         B3 = Bbb;

      else

         B1 = BB1(3)/( 2.d0*betax * DSIN( gammax - 2.d0/3.d0*pi ) )**2;
         B2 = BB2(3)/( 2.d0*betax * DSIN( gammax - 4.d0/3.d0*pi ) )**2;
         B3 = BB3(3)/( 2.d0*betax * DSIN( gammax - 6.d0/3.d0*pi ) )**2;

         if( ABS( gammax - pi/3.d0 ) .le. 1.d-4 ) B2 = Bgg;
         if( ABS( gammax - 0.d0    ) .le. 1.d-4 ) B3 = Bgg;

      endif



      ! calculating zero point energy
      vvib = ( aM003*aM222 + aM223*aM002 - 2.d0*aM023*aM022 )
     &       / ( 4.d0 * ( aM003*aM223 - aM023*aM023 ) );

      vrot = ( aMp12/aMp13 + aMn12/aMn13 + aMn22/aMn23 ) / 4.d0;










      ! printing part


      ! This printing is weird, {betac,gammac} are not
      ! even visible in this routine, and the common
      ! block "etrrqq" containing {etot,rr2,q0p,q2p}
      ! does not even exist in the original DIRHBT code.
      ! I'll just comment it out for now and ignore it...
!     if( lpr ) then
!         write(l6,400) betac , gammac*180.d0/pi ,
!    &                  etot , vvib , vrot ,
!    &                  B1   , B2   , B3   ,
!    &                  Bbb  , Bbg  , Bgg  ,
!    &                  rr2  , q0p  , q2p  ;
!     endif
! 400 format( f8.4 , f8.2 , 15f12.4 )



      if( lpr ) then

          write(l6,'(a)')
     &    '.....................................';

          write(l6,'(a,3f15.6)')
     &    ' Moment of inertia Jx................', BB1(1),BB1(2),BB1(3);

          write(l6,'(a,1f15.6)')
     &    '...................Bx................', B1;

          write(l6,'(a,3f15.6)')
     &    ' Moment of inertia Jy................', BB2(1),BB2(2),BB2(3);

          write(l6,'(a,1f15.6)')
     &    '...................By................', B2;

          write(l6,'(a,3f15.6)')
     &    ' Moment of inertia Jz................', BB3(1),BB3(2),BB3(3);

          write(l6,'(a,1f15.6)')
     &    '...................Bz................', B3;

          write(l6,'(a,3f15.6)')
     &    ' Mass parameters (Bbb,Bbg,Bgg).......', Bbb,Bbg,Bgg;

          write(l6,'(a,3f15.6)')
     &    ' Mass parameters (B00,B02,B22).......', B00,B02,B22;

          write(l6,'(a,1f15.6)')
     &    ' Rotational correction...............', vrot;

          write(l6,'(a,1f15.6)')
     &    ' Vibrational correction..............', vvib;

          write(l6,'(a)')
     &    '.....................................';



          ! printout for particles
          do it = 1 , 2

              write(l6,'(a,a)')
     &        ' single-particle exp in the canonical basis: ' , tit(it);
              write(l6,'(a)')
     &        '-------------------------------------------------------';

              do ib = 1 , nb

                  k1 = kacan(ib,it)+1;
                  k2 = kacan(ib,it)+kdcan(ib,it);

                  do k = k1 , k2
                      ajxc = ajxcan(k,it);
                      ajyc = ajycan(k,it);
                      ajzc = ajzcan(k,it);

                      write(l6,'(2i4,3f10.5)')
     &                it , k , ajxc , ajyc , ajzc;

                  enddo

              enddo

          enddo



      endif



      ! This printing is also weird, because the
      ! common block "etrrgg" containing {etot,rr2,q0p,q2p}
      ! does not even exist in the original DIRHBT code.
      ! I'll just comment it out for now and ignore it...
!     open( 49 , file = 'cpar.out'  , status = 'unknown' );
!     open( 47 , file = 'cparo.out' , status = 'unknown' );
!
!     write( 47 , 100 );
!     write( 49 , 200 );
!
!100  format(2x,'beta',4x,'gamm',10x,'Eb',10x,'Vvib',10x,'Vrot',10x,
!    &       'Bx',10x,'By',10x,'Bz',10x,'B00',10x,'B02',10x,'B22',
!    &       8x,'Rc^2',10x,'Q20(eb)',10x,'Q22(eb)');
!200  format(2x,'beta',4x,'gamm',10x,'Eb',10x,'B1',10x,'B2',10x,'B3',
!    &       10x,'Bbb',10x,'Bbg',10x,'Bgg',10x,'Rc^2',10x,'Q20(eb)',
!    &       10x,'Q22(eb)');
!
!     write(47,400) betg   , gamg   ,
!    &              etot   , vvib   , vrot   ,
!    &              BB1(3) , BB2(3) , BB3(3) ,
!    &              B00    , B02    , B22    ,
!    &              rr2    , q0p    , q2p    ;
!
!     write(49,400) betg   , gamg   ,
!    &              etot-vvib-vrot  ,
!    &              B1      , B2    , B3     ,
!    &              Bbb    , Bbg    , Bgg    ,
!    &              rr2    , q0p    , q2p    ;
!
!     close(47);
!     close(49);










      ! deallocation of local memory
      deallocate( qxx  );
      deallocate( qyy  );
      deallocate( qzz  );
      deallocate( qxy  );
      deallocate( qxz  );
      deallocate( qyz  );

      deallocate( ajx  );
      deallocate( ajy  );
      deallocate( ajz  );

      deallocate( vt   );
      deallocate( ut   );
      deallocate( ek   );
      deallocate( eqp  );
      deallocate( kk   );

      deallocate( qq0  );
      deallocate( qq2  );
      deallocate( qqx  );
      deallocate( qqy  );
      deallocate( qqz  );
      deallocate( ajjx );
      deallocate( ajjy );
      deallocate( ajjz );

      deallocate( BB1  );
      deallocate( BB2  );
      deallocate( BB3  );

      deallocate( s1  );
      deallocate( s2  );
      deallocate( s3  );
      deallocate( ss1 );
      deallocate( ss2 );
      deallocate( ss3 );
      deallocate( sp1 );
      deallocate( sn1 );
      deallocate( sn2 );










      if(lpr) then
      write(l6,*) ' ****** END INERTIA *******************************';
      endif

      return;
      end;

c======================================================================c

      subroutine  symmUtAU( n,k , alpha , A,LDA , U,LDU ,  C,LDC )

c======================================================================c
      ! Performs an update: C <-- C + alpha*( U^T * A * U )
      ! [in]:  A     is n x n symmetric matrix
      ! [in]:  U     is n x k matrix
      ! [in]:  C     is k x k symmetric matrix on input
      ! [in]:  alpha is real  scalar
      !
      ! [out]: C     is k x k symmetric matrix on output

      implicit none;

      integer n;
      integer k;
      integer LDA;
      integer LDU;
      integer LDC;

      double precision alpha;

      double precision A(LDA,n);
      double precision U(LDU,k);
      double precision C(LDC,k);

      double precision Acpy(n,n);
      double precision tmp1(n,k);
      double precision tmp2(k,k);
      integer i,j;


      call assert( n .le. LDA , 'LDA too small' );
      call assert( n .le. LDU , 'LDU too small' );
      call assert( k .le. LDC , 'LDC too small' );


      do j = 1 , n
          do i = 1 , j-1
              Acpy(i,j) = A(i,j);
          enddo
      enddo
      do i = 1 , n
          Acpy(i,i) = 0.5d0 * A(i,i);
      enddo
      do j = 1 , k
          do i = 1 , n
              tmp1(i,j) = U(i,j);
          enddo
      enddo

      ! tmp1 = alpha * UpperTriangle(A) * U
      call dtrmm( 'L','U','N','N', n,k, alpha , Acpy,n, tmp1,n );

      ! tmp2 = U^T * tmp1 = alpha * U^T * UpperTriangle(A) * U
      call dgemm( 'T','N', k,k,n, 1.d0, U,LDU, tmp1,n, 0.d0, tmp2,k );

      ! C <-- C + tmp2 + tmp2^T = C + alpha*( U^T * A * U )
      do j = 1 , k
          do i = 1 , k
              C(i,j) = C(i,j) + tmp2(i,j) + tmp2(j,i);
          enddo
      enddo


      return;
      end;

c======================================================================c

      subroutine asymmUtAU( n,k , alpha , A,LDA , U,LDU ,  C,LDC )

c======================================================================c
      ! Performs an update: C <-- C + alpha*( U^T * A * U )
      ! [in]:  A     is n x n anti-symmetric matrix
      ! [in]:  U     is n x k matrix
      ! [in]:  C     is k x k anti-symmetric matrix on input
      ! [in]:  alpha is real  scalar
      !
      ! [out]: C     is k x k anti-symmetric matrix on output

      implicit none;

      integer n;
      integer k;
      integer LDA;
      integer LDU;
      integer LDC;

      double precision alpha;

      double precision A(LDA,n);
      double precision U(LDU,k);
      double precision C(LDC,k);

      double precision Acpy(n,n);
      double precision tmp1(n,k);
      double precision tmp2(k,k);
      integer i,j;


      call assert( n .le. LDA , 'LDA too small' );
      call assert( n .le. LDU , 'LDU too small' );
      call assert( k .le. LDC , 'LDC too small' );


      do j = 1 , n
          do i = 1 , j-1
              Acpy(i,j) = A(i,j);
          enddo
      enddo
      do i = 1 , n
          Acpy(i,i) = 0.5d0 * A(i,i);
      enddo
      do j = 1 , k
          do i = 1 , n
              tmp1(i,j) = U(i,j);
          enddo
      enddo

      ! tmp1 = alpha * UpperTriangle(A) * U
      call dtrmm( 'L','U','N','N', n,k, alpha , Acpy,n, tmp1,n );

      ! tmp2 = U^T * tmp1 = alpha * U^T * UpperTriangle(A) * U
      call dgemm( 'T','N', k,k,n, 1.d0, U,LDU, tmp1,n, 0.d0, tmp2,k );

      ! C <-- C + tmp2 - tmp2^T = C + alpha*( U^T * A * U )
      do j = 1 , k
          do i = 1 , k
              C(i,j) = C(i,j) + tmp2(i,j) - tmp2(j,i);
          enddo
      enddo


      return;
      end;

c======================================================================c

      subroutine inverse2x2( A , invA )

c======================================================================c
      ! Calculates the inverse of a 2 x 2 regular matrix A
      implicit none;

      double precision A(2,2);
      double precision invA(2,2);

      double precision det;

      det = A(1,1)*A(2,2) - A(1,2)*A(2,1);

      invA(1,1) = + A(2,2) / det;
      invA(1,2) = - A(1,2) / det;
      invA(2,1) = - A(2,1) / det;
      invA(2,2) = + A(1,1) / det;

      return;
      end;

c======================================================================c

      subroutine invB_A_invB( A , B , Ans )

c======================================================================c
      ! For 2 x 2 matrices A and B, returns Ans = inv(B) * A * inv(B)
      implicit none;

      double precision   A(2,2);
      double precision   B(2,2);
      double precision Ans(2,2);

      double precision invB(2,2);

      integer i,j,k1,k2;



      call inverse2x2( B , invB );

      do i = 1 , 2
          do j = 1 , 2

              Ans(i,j) = 0.d0;
              do k1 = 1 , 2
              do k2 = 1 , 2
                  Ans(i,j) = Ans(i,j) + invB(i,k1)*A(k1,k2)*invB(k2,j);
              enddo
              enddo

          enddo
      enddo


      return;
      end;


