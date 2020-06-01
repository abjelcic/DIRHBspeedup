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
      call canon(.true.);
#ifndef ORIGINAL
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
      call resu(.true.);

      call plot(.false.);

c---- punching of potentials to tape
      call inout(2,.false.);
      call dinout(2,.false.);

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

c=====================================================================c
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
      write(6,*)'* Speedup by A.Bjelcic                              *';
      write(6,*)'*                                                   *';
      write(6,*)'* This is the original DIRHBT code (2014) directly  *';
      write(6,*)'* from CPC program library with few modifications:  *';
      write(6,*)'* 1.) The routine iter    is deleted from dirhbt.f  *';
      write(6,*)'* 2.) The routine centmas is deleted from dirhbt.f  *';
      write(6,*)'* 3.) Main is deleted from dirhbt.f                 *';
      write(6,*)'* 4.) Ref. BLAS routines are deleted from dirhbt.f  *';
      write(6,*)'* 5.) Some minor bugs in dirhbt.f are fixed         *';
      write(6,*)'* 6.) File abjelcic.f is added                      *';
      write(6,*)'* 7.) Makefile is slightly modified                 *';
      write(6,*)'* 8.) (Open)BLAS & LAPACK are required              *';
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
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        +1.d+0,      Q1(1,i),NBLOCK,
     &                                    hQ1(1,i),NBLOCK,
     &                        +1.d+0,     eij(1,1),DEGENBLOCK );
                  call dgemm( 'T','N', ndegen,ndegen,nf,
     &                        -1.d+0,      Q2(1,i),NBLOCK,
     &                                    hQ2(1,i),NBLOCK,
     &                        +1.d+0,     eij(1,1),DEGENBLOCK );

                  ! dij
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


c=====================================================================c

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




































