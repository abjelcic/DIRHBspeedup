! barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
!
!  From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)

c======================================================================c
c
      PROGRAM DIRHBZ
c
c======================================================================c
c     Relativistic Hartree-Bogoliubov theory in a axially symmetric basis
c     Main part of the code
c----------------------------------------------------------------------c
c
c-------------------------------------
      implicit real*8 (a-h,o-z)
      common /mathco/ zero,one,two,half,third,pi

      call tictoc('tic',time);
c-------------------------------------
c
c
c---- sets data
      call default(.false.)
c
c---- reads in data
      call reader(.true.)
#ifndef ORIGINAL
      call print_header();
#endif
c
c---- force-parameters
      call forces(.true.)
c
c---- Gauss-Hermite mesh points
      call gaush(two,.false.)
      call gausl(.false.)
c
c---- oscillator basis for single particle states
      call base(.false.)
c
c---- preparations
      call prep(.false.)
c
c---- wavefunctions at Gauss-Meshpoints
      call gaupol(.false.)
c
c---- initialization of the potentials
      call inout(1,.false.)
c
c---- initialization of the pairing field
      call dinout(1,.false.)
      call start(.false.)
c
c---- single-particle matix elements
      call singf(.false.)
c
c---- preparation of pairing matrix elements
#ifndef ORIGINAL
      write(6,*)'Calculation of two dimensional Talmi Moshinsky';
      write(6,*)'brackets can take some time, hang on...';
      write(6,*);
#endif
      call singd(.false.)
c
c---- coulomb and meson propagators
      call greecou(.false.)
      call greemes(.false.)
#ifndef ORIGINAL
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(initialization): ',time,'s';
      write(6,*);
#endif
      call tictoc('tic',time);
c---- iteration
      call iter(.true.)
#ifndef ORIGINAL
      write(6,*);
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(iterations    ): ',time,'s';
#endif
c
c---- transformation to the canonical basis
      call tictoc('tic',time);
      call canon(.true.)
#ifndef ORIGINAL
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(canon transf. ): ',time,'s';
#endif
c
c---- center of mass correction
      call tictoc('tic',time);
      call centmas(.false.)
#ifndef ORIGINAL
      call tictoc('toc',time);
      write(6,'(a,f10.3,a)')'  Elapsed time(centmass corr.): ',time,'s';
#endif
c
c---- results
      call resu(.true.)
      call plot(.false.)
c
c---- punching of potentials to tape  dis.wel
      call inout(2,.false.)
      call dinout(2,.false.)

c-end-DIZ
      end


#ifdef ORIGINAL
c=====================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the spherical Dirac program
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
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nneu,npro,nmas,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      data spk_min/0.1/
c
c
      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
c      CALL zeit(1,'ITER_____',6,.false.)
c
      write(l6,*) '****** BEGIN ITER **********************************'

      ii=0
      do it=1,2
         call gamma(it)
      enddo
      call broyden(.false.)
c
      do ite = 1,maxi
         ii = ite
         write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
         write( 6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
c
c
c
c------- loop over neutrons and protons
         do it = 1,2
c
c---------- calculation of the mean field Gamma
c           call gamma(it)
c
c---------- diagonalization of the Dirac-Bogolibov equation
            call dirhb(it,.false.)
c
c---------- calculation of densities in oscillator basis
            call denssh(it,.false.)
c
         enddo   ! it
c
c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb potential
         call coulom(.false.)
c
c------- calculation of expectation values
         call expect(.false.)
c
c------- potentials in r-space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         do it = 1,2
            call gamma(it)
            call delta(it,.false.)
         enddo

         call broyden(.false.)
c
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 20
            if (ic.eq.2) goto 30
         endif
c
      enddo   ! ite
   20 write(6,100) nucnam,nmas,text1,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text3,si
      goto 40
c
   30 write(6,101) nucnam,nmas,text2,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ii,text3,si
c
   40 write(l6,*) '****** END ITER ************************************'
      write( 6,*) '****** END ITER ************************************'
c     read*
c
c      CALL zeit(2,'ITER_____',6,.false.)
c
  100 format(1x,68(1h*),/,2x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(1x,a2,i4,a27,i4,a14,f17.10)
  102 format(i3,a,f10.6,3(a,f7.3),a,f5.2)
c
      return
c-end-ITER
      end
#else
c=====================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the spherical Dirac program
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
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nneu,npro,nmas,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      data spk_min/0.1/
c
c
      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
c      CALL zeit(1,'ITER_____',6,.false.)
c
      write(l6,*) '****** BEGIN ITER **********************************'


      ii=0
      do it=1,2
         call gamma(it)
      enddo
      call broyden(.false.)
c

      epsi = +5.0d-3;
      si   = +1.d+20;
      do ite = 1,maxi
         ii = ite
         write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
         write( 6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
c
c
c
c------- loop over neutrons and protons
         do it = 1,2
c
c---------- calculation of the mean field Gamma
c           call gamma(it)
c
c---------- diagonalization of the Dirac-Bogolibov equation
            call dirhb_abjelcic(it,.false.)
c
c---------- calculation of densities in oscillator basis
            call denssh(it,.false.)
c
         enddo   ! it


c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb potential
         call coulom(.false.)
c
c------- calculation of expectation values
         call expect(.false.)
c
c------- potentials in r-space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         do it = 1,2
            call gamma(it)
            call delta(it,.false.)
         enddo

         call broyden(.false.)
c
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 20
            if (ic.eq.2) goto 30
         endif
c
      enddo   ! ite


   20 write(6,100) nucnam,nmas,text1,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text3,si
      goto 40
c
   30 write(6,*)'First stage completed';
      ite1 = ite;






      ii=0
      do it=1,2
         call gamma(it)
      enddo
      call broyden(.false.)
c

      epsi = 1.d-6;
      si   = 1.d+20;
      do ite = 1,maxi
         ii = ite
         write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
         write( 6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
c
c
c
c------- loop over neutrons and protons
         do it = 1,2
c
c---------- calculation of the mean field Gamma
c           call gamma(it)
c
c---------- diagonalization of the Dirac-Bogolibov equation
            call dirhbfull_abjelcic(it,.false.)
c
c---------- calculation of densities in oscillator basis
            call denssh(it,.false.)
c
         enddo   ! it
c
c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb potential
         call coulom(.false.)
c
c------- calculation of expectation values
         call expect(.false.)
c
c------- potentials in r-space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         do it = 1,2
            call gamma(it)
            call delta(it,.false.)
         enddo

         call broyden(.false.)
c
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 21
            if (ic.eq.2) goto 31
         endif
c
      enddo   ! ite


   21 write(6,100) nucnam,nmas,text1,ite1+ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ite1+ii,text3,si
      goto 40
c
   31 write(6,101) nucnam,nmas,text2,ite1+ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ite1+ii,text3,si
c
















   40 write(l6,*) '****** END ITER ************************************'
      write( 6,*) '****** END ITER ************************************'
c     read*
c
c      CALL zeit(2,'ITER_____',6,.false.)
c
  100 format(1x,68(1h*),/,2x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(1x,a2,i4,a27,i4,a14,f17.10)
  102 format(i3,a,f10.6,3(a,f7.3),a,f5.2)
c
      return
c-end-ITER
      end
#endif

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
      write(6,*)'* This is the original DIRHBZ code (2014) directly  *';
      write(6,*)'* from CPC program library with few modifications:  *';
      write(6,*)'* 1.) The subroutine iter is deleted from dirhbz.f  *';
      write(6,*)'* 2.) The function talmos2 is deleted from dirhbz.f *';
      write(6,*)'* 3.) Main is deleted from dirhbz.f                 *';
      write(6,*)'* 4.) Ref. BLAS routines are deleted from dirhbz.f  *';
      write(6,*)'* 5.) Some minor bugs in dirhbz.f are fixed         *';
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

#ifdef ORIGINAL
c=======================================================================c
c
      real*8 function talmos2(im1,in1,in2,im3,in3,in4)
c
c=======================================================================c
c
c	2d-moshinsky bracket:
c
c	<n1 m1, n2 -m1 | n3 m3, n4 -m3>
c
c	radial quantum number start from zero: n=0,1,2,....
c	orbital angular momentum m1=-m2>0,m3=-m4>0
c
c-----------------------------------------------------------------------c
C
c     iv(n)  =  (-1)**n
c     fak(n) =  n!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
C
C-----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      common /gfviv / iv(-IGFV:IGFV)
      common /gfvfak/ fak(0:IGFV)
      common /gfvfi / fi(0:IGFV)
      common /gfvwf / wf(0:IGFV)
      common /gfvwfi/ wfi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
c
      im2 = -im1
      im4 = -im3
c
      talmos2 = zero
c
      nn12 = 2*in1+iabs(im1)+2*in2+iabs(im2)
      nn34 = 2*in3+iabs(im3)+2*in4+iabs(im4)
      if (im1.lt.0)     stop 'in TALMOS2: m1 < 0'
      if (im3.lt.0)     stop 'in TALMOS2: m3 < 0'
      if (nn12.ne.nn34) return
c
      s34 = zero
c
      if (im1.lt.im2) then
         n1 = in2
         n2 = in1
         m1 = im2
         m2 = im1
      else
         n1 = in1
         n2 = in2
	 m1 = im1
         m2 = im2
      endif
      n3 = in3
      n4 = in4
      m3 = im3
      m4 = im4
      nn1 = 2*n1 + abs(m1)
      nn2 = 2*n2 + abs(m2)
      nn3 = 2*n3 + abs(m3)
      nn4 = 2*n4 + abs(m4)
      if (m3.gt.m4) then
         if (n1+n2.ne.n3+n4+m2-m4) return
      else
         if (n1+n2.ne.n3+n4+m2-m3) return
      endif
      prout = iv(n3+n4-n1-n2)/sqrt(two**(nn3+nn4))*
     &	wf(n1)*wf(n1+abs(m1))*
     &	wf(n2)*wf(n2+abs(m2))*
     &	wfi(n3)*wfi(n3+abs(m3))*
     &	wfi(n4)*wfi(n4+abs(m4))
c
      sn3 = zero
      sn4 = zero
      sm3 = zero
      sm4 = zero
c
      do i3 = 0,n3
      do j3 = 0,n3
      do k3 = 0,n3
c
         l3 = n3 - i3 - j3 - k3
c
      do i4 = 0,n4
      do j4 = 0,n4
      do 20 k4 = 0,n4
c
         l4 = n4 - i4 - j4 - k4
c
         if (l3.lt.0.or.l4.lt.0) goto 20
c
         do it3 = 0,abs(m3)
         do 10 it4 = 0,abs(m4)
c
         if (m3.gt.m4) then
            if (i3+i4+j3+j4+it3.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4-it3+it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it3) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         else
            if (i3+i4+j3+j4+it4.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4+it3-it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it4) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         endif ! m3,m4

   10 continue   ! it4 modified by tamara
      enddo   ! it3
c
   20 continue   ! k2
      enddo   ! j2
      enddo   ! i2
c
      enddo   ! k1
      enddo   ! j1
      enddo   ! i1
c
      talmos2 = s34*prout
c
      return
c-end-TALMOS2
      end
#else
c=======================================================================c
c
      real*8 function talmos2(im1,in1,in2,im3,in3,in4)
c
c=======================================================================c
c
c	2d-moshinsky bracket:
c
c	<n1 m1, n2 -m1 | n3 m3, n4 -m3>
c
c	radial quantum number start from zero: n=0,1,2,....
c	orbital angular momentum m1=-m2>0,m3=-m4>0
c
c-----------------------------------------------------------------------c
C
c     iv(n)  =  (-1)**n
c     fak(n) =  n!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
C
C-----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      common /gfviv / iv(-IGFV:IGFV)
      common /gfvfak/ fak(0:IGFV)
      common /gfvfi / fi(0:IGFV)
      common /gfvwf / wf(0:IGFV)
      common /gfvwfi/ wfi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
c
      im2 = -im1
      im4 = -im3
c
      talmos2 = zero
c
      nn12 = 2*in1+iabs(im1)+2*in2+iabs(im2)
      nn34 = 2*in3+iabs(im3)+2*in4+iabs(im4)
      if (im1.lt.0)     stop 'in TALMOS2: m1 < 0'
      if (im3.lt.0)     stop 'in TALMOS2: m3 < 0'
      if (nn12.ne.nn34) return
c
      s34 = zero
c
      if (im1.lt.im2) then
         n1 = in2
         n2 = in1
         m1 = im2
         m2 = im1
      else
         n1 = in1
         n2 = in2
	 m1 = im1
         m2 = im2
      endif
      n3 = in3
      n4 = in4
      m3 = im3
      m4 = im4
      nn1 = 2*n1 + abs(m1)
      nn2 = 2*n2 + abs(m2)
      nn3 = 2*n3 + abs(m3)
      nn4 = 2*n4 + abs(m4)
      if (m3.gt.m4) then
         if (n1+n2.ne.n3+n4+m2-m4) return
      else
         if (n1+n2.ne.n3+n4+m2-m3) return
      endif
      prout = iv(n3+n4-n1-n2)/sqrt(two**(nn3+nn4))*
     &	wf(n1)*wf(n1+abs(m1))*
     &	wf(n2)*wf(n2+abs(m2))*
     &	wfi(n3)*wfi(n3+abs(m3))*
     &	wfi(n4)*wfi(n4+abs(m4))
c
      sn3 = zero
      sn4 = zero
      sm3 = zero
      sm4 = zero
c
      do i3 = 0,n3
      do j3 = 0,n3-i3 !abjelcic
      do k3 = 0,n3-i3-j3 !abjelcic
c
         l3 = n3 - i3 - j3 - k3
c
      do i4 = 0,n4
      do j4 = 0,n4-i4 !abjelcic
      do 20 k4 = 0,n4-i4-j4 !abjelcic
c
         l4 = n4 - i4 - j4 - k4
c
         if (l3.lt.0.or.l4.lt.0) goto 20
c
         do it3 = 0,abs(m3)
         do 10 it4 = 0,abs(m4)
c
         if (m3.gt.m4) then
            if (i3+i4+j3+j4+it3.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4-it3+it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it3) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         else
            if (i3+i4+j3+j4+it4.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4+it3-it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it4) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         endif ! m3,m4

   10 continue   ! it4 modified by tamara
      enddo   ! it3
c
   20 continue   ! k2
      enddo   ! j2
      enddo   ! i2
c
      enddo   ! k1
      enddo   ! j1
      enddo   ! i1
c
      talmos2 = s34*prout
c
      return
c-end-TALMOS2
      end
#endif











c======================================================================c

      subroutine dirhb_abjelcic(it,lpr)

c======================================================================c
c
c     solves the RHB-Equation
c     IT    = 1 for neutrons
c     IT    = 2 for protons
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr,lprl
c
      character*1 bbb
      character*8 tbb(NHBX)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tb*6                                            ! blokap
      character tt*8                                            ! quaosc
      character nucnam*2                                        ! nucnuc
c
      dimension hb(NHBQX),e(NHBX),ez(NHBX)
c
      common /blodir/ ka(NBX,4),kd(NBX,4)
      common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      common /bloosc/ ia(NBX,2),id(NBX,2)
      common /deldel/ de(NHHX,NB2X)
      common /fermi / ala(2),tz(2)
      common /gamgam/ hh(NHHX,NB2X)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /quaosc/ nt,nnn(NTX,5),tt(NTX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)
c
      data maxl/200/,epsl/1.d-8/,bbb/'-'/,lprl/.false./
c


c---------------------------------------------------------------------abjelcic-inserted-start
      integer NBLOCK;
      parameter( NBLOCK = ((N0FX+2)*(N0FX+3))/2 );
      double precision h         ( NBLOCK , NBLOCK , NBX );
      double precision diagh     ( NBLOCK ,          NBX );
      double precision Delta     ( NBLOCK , NBLOCK , NBX );
      double precision tempmat   ( NBLOCK , NBLOCK       );
      double precision diagDelta ( NBLOCK ,          NBX );
      double precision h2plDelta2( NBLOCK , NBLOCK       );
      double precision asymhDelta( NBLOCK , NBLOCK       );
      double complex         hDhD( NBLOCK , NBLOCK , NBX );
      double complex   hDhDlambda( NBLOCK , NBLOCK       );
      double precision Q1        ( NBLOCK , NBLOCK       );
      double precision Q2        ( NBLOCK , NBLOCK       );
      double precision ediag     ( NBLOCK                );
      double precision ddiag     ( NBLOCK                );
      double precision hmilambda ( NBLOCK , NBLOCK       );
      double precision hQ1       ( NBLOCK , NBLOCK       );
      double precision DeltaQ1   ( NBLOCK , NBLOCK       );
      double precision trhQ2     ( NBLOCK , NBLOCK       );
      double precision trDeltaQ2 ( NBLOCK , NBLOCK       );
      double precision cosdiag   ( NBLOCK                );
      double precision sindiag   ( NBLOCK                );
      double precision Epldiag   ( NBLOCK          , NBX );
      double precision Upl       ( NBLOCK , NBLOCK , NBX );
      double precision Vpl       ( NBLOCK , NBLOCK , NBX );
      double precision qUpl      ( NBLOCK , NBLOCK       );
      double precision qVpl      ( NBLOCK , NBLOCK       );
      external         ddot;
      double precision ddot;
      double precision e1,e2,e3;
      double precision d1,d2,d3;
      double precision ei,di,Epli;


      external         DLAMCH;
      double precision DLAMCH;
      double precision ABSTOL;
      integer          NFOUND;
      double precision W(NBLOCK);
      double complex   Z(NBLOCK,NBLOCK);
      integer          ISUPPZ( 2*max(1,NBLOCK) );

      integer          LWORK_ASK;
      integer          LWORK;
      parameter( LWORK = NBLOCK*NBLOCK );
      double complex   WORK(LWORK);

      integer          LRWORK_ASK;
      integer          LRWORK;
      parameter( LRWORK = NBLOCK*NBLOCK );
      double precision RWORK(LRWORK);

      integer          LIWORK_ASK;
      integer          LIWORK;
      parameter( LIWORK = NBLOCK*NBLOCK );
      integer          IWORK(LIWORK);

      integer          INFO;


      double precision Wh(NBLOCK);
      double precision Zh(NBLOCK,NBLOCK,NBX);

      integer          LWORKh_ASK;
      integer          LWORKh;
      parameter( LWORKh = NBLOCK*NBLOCK );
      double precision WORKh(LWORKh);

      integer          LIWORKh_ASK;
      integer          LIWORKh;
      parameter( LIWORKh = NBLOCK*NBLOCK );
      integer          IWORKh(LIWORKh);
c---------------------------------------------------------------------abjelcic-inserted-end





      if (lpr) then
      write(l6,*) ' ****** BEGIN DIRHB ********************************'
      endif

      dl    = 100.d0
      xh    = ala(it) + dl
      xl    = ala(it) - dl
      al    = ala(it)
c
c======================================================================c
      do lit = 1,maxl         ! loop over lambda-iteration
c======================================================================c
         snold = sn
         sn  = zero
         klp = 0
         kla = 0
c======================================================================c
      do ib = 1,nb            ! loop over differnt blocks
c======================================================================c
         mul  = mb(ib)
         nf   = id(ib,1)
         ng   = id(ib,2)
         nh   = nf + ng
         nhb  = nh + nh
         m    = ib + (it-1)*NBX







c---------------------------------------------------------------------abjelcic-inserted-start
         if( lit .eq. 1 ) then

         call assert( nh.le.NBLOCK , 'nh > NBLOCK in dirhb()' );


         ! Spectral decomposition of symmetric matrix h
         ! Keep in mind that h gets destroyed
         do j = 1 , nh
             do i = j , nh
                  h(i,j,ib) = hh( i + (j-1)*nh , m );
             enddo
         enddo
         ABSTOL     = DLAMCH( 'Safe minimum' );
         LWORKh_ASK  = -1;
         LIWORKh_ASK = -1;
         call dsyevr('V','A','L',  nh,  h(1,1,ib),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nh, ABSTOL,  NFOUND,
     &                Wh,  Zh(1,1,ib),NBLOCK,
     &                ISUPPZ,
     &                WORKh,  LWORKh_ASK,
     &                IWORKh, LIWORKh_ASK,
     &                INFO                                     );
         call assert(INFO .eq. 0  , 'INFO != 0 in dirhb()' );
         call assert(INT( WORKh(1)+0.5d0).le.LWORKh ,'LWORKh small' );
         call assert(INT(IWORKh(1)+0.5d0).le.LIWORKh,'LIWORKh small');

         call dsyevr('V','A','L',  nh,  h(1,1,ib),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nh, ABSTOL,  NFOUND,
     &                Wh,  Zh(1,1,ib),NBLOCK,
     &                ISUPPZ,
     &                WORKh,  LWORKh,
     &                IWORKh, LIWORKh,
     &                INFO                                     );
         call assert(   INFO.eq.0  ,   'INFO != 0  in dirhb()' );
         call assert( NFOUND.eq.nh , 'NFOUND != nh in dirhb()' );


         do j = 1 , nf
             do i = j , nf
                 Delta(i,j,ib) = de( i + (j-1)*nh , m );
             enddo
         enddo
         call dsymm('L','L',  nf,nf,
     &              +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &                             Zh(1,ng+1,ib),NBLOCK,
     &              +0.0d+0,        tempmat(1,1),NBLOCK  );
         call dgemm('T','N',  nf,nf,nf,
     &              +1.0d+0,       Zh(1,ng+1,ib),NBLOCK,
     &                              tempmat(1,1),NBLOCK,
     &              +0.0d+0,       Delta(1,1,ib),NBLOCK  );

         do j = 1 , nf
                 diagh(j,ib) = Wh(ng+j);
             diagDelta(j,ib) = Delta(j,j,ib);
         enddo


         ! h^2+Delta^2
         call dsyrk( 'L','N',  nf,nf,
     &              +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &              +0.0d+0,  h2plDelta2(1,1   ),NBLOCK );
         do i = 1 , nf
             h2plDelta2(i,i) = h2plDelta2(i,i) + diagh(i,ib)**2.d+0;
         enddo
         !Delta*h-h*Delta
         do j = 1 , nf
             do i = j , nf
                 asymhDelta(i,j) = Delta(i,j,ib)
     &                            * ( diagh(j,ib) - diagh(i,ib) );
             enddo
         enddo
         ! hDhD = ( h + 1jDelta )*( h - 1jDelta )
         do j = 1 , nf
             do i = j , nf
                 hDhD(i,j,ib) = COMPLEX( h2plDelta2(i,j)  ,
     &                                   asymhDelta(i,j) );
             enddo
         enddo




         ! Filling arrays required by the rest of the DIRHB code
         ka(ib,it+2) = kla;
         do k = 1 , ng
             kla = kla + 1;

             equ(kla,it+2) = Wh(ng-k+1)-al;
             do n = 1 , nh
                 fguv(n+ 0,kla,it+2) = Zh(n,ng-k+1,ib);
                 fguv(n+nh,kla,it+2) = +0.d+0;
             enddo

         enddo
         kd(ib,it+2) = kla - ka(ib,it+2)



         endif







         !==============================================================!
         != HFB eigensolve method starts here ==========================!
         !==============================================================!

         ! hDhDlambda = ( h-lambdaI + 1jDelta )*( h-lambdaI - 1jDelta )
         do j = 1 , nf
             do i = j , nf
                 hDhDlambda(i,j) = hDhD(i,j,ib);
             enddo
             hDhDlambda(j,j) = hDhDlambda(j,j) + al*al
     &                                         - 2.0d+0*al*diagh(j,ib);
         enddo

         ! Spectral decomposition of hermitian matrix hDhDlambda
         ! Keep in mind that hDhDlambda gets destroyed
         ABSTOL     = DLAMCH( 'Safe minimum' );
         LWORK_ASK  = -1;
         LRWORK_ASK = -1;
         LIWORK_ASK = -1;
         call zheevr('V','A','L',  nf,  hDhDlambda(1,1),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nf, ABSTOL,  NFOUND,
     &                W,  Z,NBLOCK,
     &                ISUPPZ,
     &                WORK,  LWORK_ASK,
     &                RWORK, LRWORK_ASK,
     &                IWORK, LIWORK_ASK,
     &                INFO                                     );
         call assert( INFO .eq. 0  , 'INFO != 0 in dirhb()' );
         call assert( INT( WORK(1)+0.5d0).le.LWORK ,'LWORK too small' );
         call assert( INT(RWORK(1)+0.5d0).le.LRWORK,'LRWORK too small');
         call assert( INT(IWORK(1)+0.5d0).le.LIWORK,'LIWORK too small');

         call zheevr('V','A','L',  nf,  hDhDlambda(1,1),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nf, ABSTOL,  NFOUND,
     &                W,  Z,NBLOCK,
     &                ISUPPZ,
     &                WORK,  LWORK,
     &                RWORK, LRWORK,
     &                IWORK, LIWORK,
     &                INFO                                     );
         call assert(   INFO.eq.0  ,   'INFO != 0  in dirhb()' );
         call assert( NFOUND.eq.nf , 'NFOUND != nf in dirhb()' );

         ! Q1 and Q2
         do j = 1 , nf
             do i = 1 , nf
                 Q1(i,j) = DREAL(Z(i,j));
                 Q2(i,j) = DIMAG(Z(i,j));

                 trDeltaQ2(i,j) = Q2(i,j);
             enddo
         enddo

         ! DeltaQ1 and trDeltaQ2
         call dsymm('L','L', nf,nf,
     &               +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                            Q1(1,1   ),NBLOCK,
     &               +0.d+0, DeltaQ1(1,1   ),NBLOCK  );
         call dtrmm('L','L','N','N', nf,nf,
     &               +1.0d+0,     Delta(1,1,ib),NBLOCK,
     &                        trDeltaQ2(1,1   ),NBLOCK );

         ! ediag and ddiag
         do i = 1 , nf
             e1 = +2.0d0 * ddot( nf, DeltaQ1(1,i),1, Q2(1,i),1 );
             d1 = +0.0d+0;
             do j = 1 , nf
                 d1 = d1 + (diagh(j,ib)-al)*Q2(j,i)*Q1(j,i);
             enddo
             d1 = -2.0d0 * d1;

             e2 = +0.0d+0;
             do j = 1 , nf
                 e2 = e2 + (diagh(j,ib)-al)*Q1(j,i)*Q1(j,i);
             enddo
             d2 = ddot( nf, DeltaQ1(1,i),1, Q1(1,i),1 );

             e3 = +0.0d+0;
             do j = 1 , nf
                 e3 = e3 + (diagh(j,ib)-al)*Q2(j,i)*Q2(j,i);
             enddo
             e3 = -e3;
             d3 = +0.0d+0;
             do j = 1 , nf
                 d3 = d3 +  diagDelta(j,ib)*Q2(j,i)*Q2(j,i);
             enddo
             d3 = d3 - 2.0d0 * ddot( nf, trDeltaQ2(1,i),1, Q2(1,i),1 );


             ediag(i) = e1 + e2 + e3;
             ddiag(i) = d1 + d2 + d3;
         enddo

         ! Ediag, cosdiag and sindiag
         do i = 1 , nf
             ei   = ediag(i);
             di   = ddiag(i);
             Epli = DSQRT( ei**2.d0 + di**2.d0 );

             !call assert(ABS(Epli-W(i)).lt.1.D-7,'eigenvalues error');

             Epldiag(i,ib) = Epli;
             cosdiag(i) = (Epli+ei)/DSQRT((Epli+ei)**2.d0+di**2.d0);
             sindiag(i) = (di     )/DSQRT((Epli+ei)**2.d0+di**2.d0);
         enddo

         do j = 1 , nf
             do i = 1 , nf
                 Upl(i,j,ib) = Q1(i,j)*cosdiag(j) - Q2(i,j)*sindiag(j);
                 Vpl(i,j,ib) = Q1(i,j)*sindiag(j) + Q2(i,j)*cosdiag(j);
             enddo
         enddo

         !==============================================================!
         != HFB eigensolve method ends here ============================!
         !==============================================================!


         ! Calculating number of particles Z/N for given lambda
         snib = +0.0d+0;
         do k = 1 , nf
             do n = 1 , nf
                 snib = snib + Vpl(n,k,ib)**2.d0;
             enddo
         enddo
         snib = +2.0d0 * snib; ! Time reversal symmetry

         sn = sn + snib;
c---------------------------------------------------------------------abjelcic-inserted-end










c======================================================================c
      enddo   ! ib
c======================================================================c


      if (lit.gt.1) dd = (sn - snold)/(al - alold)
c------- calculation of a new lambda-value
      alold = al
      dn    = sn - tz(it)
      if (dn.lt.zero) then
          xl = al
      else
          xh = al
      endif
      if (lit.eq.1) then
         if(dabs(dn).le.0.1d0) then
            al = al - dn
         else
             al = al - 0.1d0*sign(one,dn)
         endif
      else
c           secant method
         if (dd.eq.zero) dd = 1.d-20
         al    = al - dn/dd
         if (al.lt.xl.or.al.gt.xh) then
c              bisection
            al = half*(xl+xh)
            bbb = 'B'
         endif
      endif

      if( dabs(sn-tz(it)).lt.min(0.1d0,0.01d0*si) ) goto 30;
      if(abs(al-alold).lt.epsl) goto 30

      if (lprl.or.lit.gt.10) then
          write(l6,113) lit,'. L-Iteration: ',bbb,alold,dn,al
          write(6,113)  lit,'. L-Iteration: ',bbb,alold,dn,al
          bbb = ' '
      endif
c
c---- end of lambda-loop
      enddo




      write(l6,'(a,i4,a)')
     &     ' Lambda-Iteration interupted after',lit-1,' steps'
      stop
   30 if (lprl) then
         write(l6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
c        write(6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
      endif
      ala(it) = al
#ifdef VERBOSE
      write(6,'(a,i3)')'Number of lambda iterations: ', lit;
#endif








      ! Filling arrays required by the rest of the DIRHB code
      klp = 0;
      do ib = 1 , nb
         nf = id(ib,1);
         ng = id(ib,2);
         nh = nf + ng;

         call dgemm('N','N',nf+ng,nf,nf, +1.0d+0, Zh(1,ng+1,ib),NBLOCK,
     &                                              Upl(1,1,ib),NBLOCK,
     &                                   +0.0d+0,     qUpl(1,1),NBLOCK);
         call dgemm('N','N',nf+ng,nf,nf, +1.0d+0, Zh(1,ng+1,ib),NBLOCK,
     &                                              Vpl(1,1,ib),NBLOCK,
     &                                   +0.0d+0,     qVpl(1,1),NBLOCK);

         ka(ib,it) = klp;
         do k = 1 , nf
             klp = klp + 1;

             equ(klp,it) = + Epldiag(k,ib);
             do n = 1 , nh
                 fguv(n+ 0,klp,it) = qUpl(n,k);
                 fguv(n+nh,klp,it) = qVpl(n,k);
             enddo

         enddo
         kd(ib,it) = klp - ka(ib,it);

      enddo







      if (lpr) then
      write(l6,*) ' ****** END DIRHB **********************************'
      endif
c
  101 format(i4,a,i4,3f13.8)
  113 format(i4,a,a1,3f13.8)
c
      return
C-end-DIRHB
      end

c======================================================================c

      subroutine dirhbfull_abjelcic(it,lpr)

c======================================================================c
c
c     solves the RHB-Equation
c     IT    = 1 for neutrons
c     IT    = 2 for protons
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr,lprl
c
      character*1 bbb
      character*8 tbb(NHBX)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tb*6                                            ! blokap
      character tt*8                                            ! quaosc
      character nucnam*2                                        ! nucnuc
c
      dimension hb(NHBQX),e(NHBX),ez(NHBX)
c
      common /blodir/ ka(NBX,4),kd(NBX,4)
      common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      common /bloosc/ ia(NBX,2),id(NBX,2)
      common /deldel/ de(NHHX,NB2X)
      common /fermi / ala(2),tz(2)
      common /gamgam/ hh(NHHX,NB2X)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /quaosc/ nt,nnn(NTX,5),tt(NTX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)
c
      data maxl/200/,epsl/1.d-8/,bbb/'-'/,lprl/.false./
c


c---------------------------------------------------------------------abjelcic-inserted-start
      integer NBLOCK;
      parameter( NBLOCK = ((N0FX+2)*(N0FX+3))/2 );
      double precision h         ( NBLOCK , NBLOCK , NBX );
      double precision diagh     ( NBLOCK ,          NBX );
      double precision Delta     ( NBLOCK , NBLOCK , NBX );
      double precision diagDelta ( NBLOCK ,          NBX );
      double precision h2plDelta2( NBLOCK , NBLOCK       );
      double precision asymhDelta( NBLOCK , NBLOCK       );
      double complex         hDhD( NBLOCK , NBLOCK , NBX );
      double complex   hDhDlambda( NBLOCK , NBLOCK       );
      double precision Q1        ( NBLOCK , NBLOCK       );
      double precision Q2        ( NBLOCK , NBLOCK       );
      double precision ediag     ( NBLOCK                );
      double precision ddiag     ( NBLOCK                );
      double precision hmilambda ( NBLOCK , NBLOCK       );
      double precision hQ1       ( NBLOCK , NBLOCK       );
      double precision DeltaQ1   ( NBLOCK , NBLOCK       );
      double precision trhQ2     ( NBLOCK , NBLOCK       );
      double precision trDeltaQ2 ( NBLOCK , NBLOCK       );
      double precision cosdiag   ( NBLOCK                );
      double precision sindiag   ( NBLOCK                );
      double precision Epldiag   ( NBLOCK                );
      double precision Upl       ( NBLOCK , NBLOCK       );
      double precision Vpl       ( NBLOCK , NBLOCK       );
      external         ddot;
      double precision ddot;
      double precision e1,e2,e3;
      double precision d1,d2,d3;
      double precision ei,di,Epli;


      external         DLAMCH;
      double precision DLAMCH;
      double precision ABSTOL;
      integer          NFOUND;
      double precision W(NBLOCK);
      double complex   Z(NBLOCK,NBLOCK);
      integer          ISUPPZ( 2*max(1,NBLOCK) );

      integer          LWORK_ASK;
      integer          LWORK;
      parameter( LWORK = NBLOCK*NBLOCK );
      double complex   WORK(LWORK);

      integer          LRWORK_ASK;
      integer          LRWORK;
      parameter( LRWORK = NBLOCK*NBLOCK );
      double precision RWORK(LRWORK);

      integer          LIWORK_ASK;
      integer          LIWORK;
      parameter( LIWORK = NBLOCK*NBLOCK );
      integer          IWORK(LIWORK);

      integer          INFO;
c---------------------------------------------------------------------abjelcic-inserted-end





      if (lpr) then
      write(l6,*) ' ****** BEGIN DIRHB ********************************'
      endif

      dl    = 100.d0
      xh    = ala(it) + dl
      xl    = ala(it) - dl
      al    = ala(it)
c
c======================================================================c
      do lit = 1,maxl         ! loop over lambda-iteration
c======================================================================c
         snold = sn
         sn  = zero
         klp = 0
         kla = 0
c======================================================================c
      do ib = 1,nb            ! loop over differnt blocks
c======================================================================c
         mul  = mb(ib)
         nf   = id(ib,1)
         ng   = id(ib,2)
         nh   = nf + ng
         nhb  = nh + nh
         m    = ib + (it-1)*NBX







c---------------------------------------------------------------------abjelcic-inserted-start
         if( lit .eq. 1 ) then

             call assert( nh.le.NBLOCK , 'nh > NBLOCK in dirhb()' );

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
                     diagh(j,ib) =     h(j,j,ib);
                 diagDelta(j,ib) = Delta(j,j,ib);
             enddo

             ! h^2+Delta^2
             call dsyrk( 'L','N',  nh,nh,
     &                  +1.0d+0,           h(1,1,ib),NBLOCK,
     &                  +0.0d+0,  h2plDelta2(1,1   ),NBLOCK );
             call dsyrk( 'L','N',  nh,nh,
     &                  +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &                  +1.0d+0,  h2plDelta2(1,1   ),NBLOCK );

             !Delta*h-h*Delta
             call dsymm('L','L',  nh,nh,
     &                  +1.0d+0,       Delta(1,1,ib),NBLOCK,
     &                                     h(1,1,ib),NBLOCK,
     &                  +0.0d+0,  asymhDelta(1,1   ),NBLOCK  );
             do j = 1 , nh
                 do i = j , nh
                     asymhDelta(i,j) = + asymhDelta(i,j)
     &                                 - asymhDelta(j,i);
                 enddo
             enddo

             ! hDhD = ( h + 1jDelta )*( h - 1jDelta )
             do j = 1 , nh
                 do i = j , nh
                     hDhD(i,j,ib) = COMPLEX( h2plDelta2(i,j)  ,
     &                                       asymhDelta(i,j) );
                 enddo
             enddo

         endif

         !==============================================================!
         != HFB eigensolve method starts here ==========================!
         !==============================================================!

         ! hDhDlambda = ( h-lambdaI + 1jDelta )*( h-lambdaI - 1jDelta )
         do j = 1 , nh
             do i = j , nh
                 hDhDlambda(i,j) = hDhD(i,j,ib) - 2.0d0*al*h(i,j,ib)
             enddo
             hDhDlambda(j,j) = hDhDlambda(j,j) + al*al;
         enddo

         ! Spectral decomposition of hermitian matrix hDhDlambda
         ! Keep in mind that hDhDlambda gets destroyed
         ABSTOL     = DLAMCH( 'Safe minimum' );
         LWORK_ASK  = -1;
         LRWORK_ASK = -1;
         LIWORK_ASK = -1;
         call zheevr('V','A','L',  nh,  hDhDlambda(1,1),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nh, ABSTOL,  NFOUND,
     &                W,  Z,NBLOCK,
     &                ISUPPZ,
     &                WORK,  LWORK_ASK,
     &                RWORK, LRWORK_ASK,
     &                IWORK, LIWORK_ASK,
     &                INFO                                     );
         call assert( INFO .eq. 0  , 'INFO != 0 in dirhb()' );
         call assert( INT( WORK(1)+0.5d0).le.LWORK ,'LWORK too small' );
         call assert( INT(RWORK(1)+0.5d0).le.LRWORK,'LRWORK too small');
         call assert( INT(IWORK(1)+0.5d0).le.LIWORK,'LIWORK too small');

         call zheevr('V','A','L',  nh,  hDhDlambda(1,1),NBLOCK,
     &               +0.d+0, +1.d+10,  1,nh, ABSTOL,  NFOUND,
     &                W,  Z,NBLOCK,
     &                ISUPPZ,
     &                WORK,  LWORK,
     &                RWORK, LRWORK,
     &                IWORK, LIWORK,
     &                INFO                                     );
         call assert(   INFO.eq.0  ,   'INFO != 0  in dirhb()' );
         call assert( NFOUND.eq.nh , 'NFOUND != nh in dirhb()' );

         ! Q1 and Q2
         do j = 1 , nh
             do i = 1 , nh
                 Q1(i,j) = DREAL(Z(i,j));
                 Q2(i,j) = DIMAG(Z(i,j));

                 trhQ2(i,j)     = Q2(i,j);
                 trDeltaQ2(i,j) = Q2(i,j);

                 hmilambda(i,j) = h(i,j,ib);
             enddo
             hmilambda(j,j) = hmilambda(j,j) - al;
         enddo

         ! hQ1, DeltaQ1, trhQ2 and trDeltaQ2
         call dsymm('L','L', nh,nh,
     &               +1.d+0, hmilambda(1,1),NBLOCK,
     &                              Q1(1,1),NBLOCK,
     &               +0.d+0,       hQ1(1,1),NBLOCK  );
         call dsymm('L','L', nh,nh,
     &               +1.d+0,   Delta(1,1,ib),NBLOCK,
     &                            Q1(1,1   ),NBLOCK,
     &               +0.d+0, DeltaQ1(1,1   ),NBLOCK  );
         call dtrmm('L','L','N','N', nh,nh,
     &               +1.0d+0, hmilambda(1,1),NBLOCK,
     &                            trhQ2(1,1),NBLOCK );
         call dtrmm('L','L','N','N', nh,nh,
     &               +1.0d+0,     Delta(1,1,ib),NBLOCK,
     &                        trDeltaQ2(1,1   ),NBLOCK );

         ! ediag and ddiag
         do i = 1 , nh
             e1 = +2.0d0 * ddot( nh, DeltaQ1(1,i),1, Q2(1,i),1 );
             d1 = -2.0d0 * ddot( nh,     hQ1(1,i),1, Q2(1,i),1 );

             e2 = ddot( nh,     hQ1(1,i),1, Q1(1,i),1 );
             d2 = ddot( nh, DeltaQ1(1,i),1, Q1(1,i),1 );

             e3 = +0.0d+0;
             d3 = +0.0d+0;
             do j = 1 , nh
                 e3 = e3 + (diagh(j,ib)-al)*Q2(j,i)*Q2(j,i);
                 d3 = d3 +  diagDelta(j,ib)*Q2(j,i)*Q2(j,i);
             enddo
             e3 = e3 - 2.0d0 * ddot( nh,     trhQ2(1,i),1, Q2(1,i),1 );
             d3 = d3 - 2.0d0 * ddot( nh, trDeltaQ2(1,i),1, Q2(1,i),1 );


             ediag(i) = e1 + e2 + e3;
             ddiag(i) = d1 + d2 + d3;
         enddo

         ! Ediag, cosdiag and sindiag
         do i = 1 , nh
             ei   = ediag(i);
             di   = ddiag(i);
             Epli = DSQRT( ei**2.d0 + di**2.d0 );

             !call assert(ABS(Epli-W(i)).lt.1.D-7,'eigenvalues error');

               Epldiag(i) = Epli;
             cosdiag(i) = (Epli+ei)/DSQRT((Epli+ei)**2.d0+di**2.d0);
             sindiag(i) = (di     )/DSQRT((Epli+ei)**2.d0+di**2.d0);
         enddo

         do j = 1 , nh
             do i = 1 , nh
                 Upl(i,j) = Q1(i,j)*cosdiag(j) - Q2(i,j)*sindiag(j);
                 Vpl(i,j) = Q1(i,j)*sindiag(j) + Q2(i,j)*cosdiag(j);
             enddo
         enddo

         !==============================================================!
         != HFB eigensolve method ends here ============================!
         !==============================================================!


         ! Filling arrays required by the rest of the DIRHB code
         ka(ib,it) = klp;
         do k = 1 , nf
             klp = klp + 1;

             equ(klp,it) = + Epldiag(k);
             do n = 1 , nh
                 fguv(n+ 0,klp,it) = Upl(n,k);
                 fguv(n+nh,klp,it) = Vpl(n,k);
             enddo

         enddo
         kd(ib,it) = klp - ka(ib,it);

         ka(ib,it+2) = kla;
         do k = 1 , ng
             kla = kla + 1;

             equ(kla,it+2) = - Epldiag(nf+k);
             do n = 1 , nh
                 fguv(n+ 0,kla,it+2) = - Vpl(n,nf+k);
                 fguv(n+nh,kla,it+2) = + Upl(n,nf+k);
             enddo

         enddo
         kd(ib,it+2) = kla - ka(ib,it+2)


         ! Calculating number of particles Z/N for given lambda
         snib = +0.0d+0;
         do k = 1 , nf
             do n = 1 , nh
                 snib = snib + (+Vpl(n,k+ 0))**2.d0;
             enddo
         enddo
         do k = 1 , ng
             do n = 1 , nh
                 snib = snib + (+Upl(n,k+nf))**2.d0;
             enddo
         enddo
         snib = +2.0d0 * snib; ! Time reversal symmetry

         sn = sn + snib;
c---------------------------------------------------------------------abjelcic-inserted-end










c======================================================================c
      enddo   ! ib
c======================================================================c


      if (lit.gt.1) dd = (sn - snold)/(al - alold)
c------- calculation of a new lambda-value
      alold = al
      dn    = sn - tz(it)
      if (dn.lt.zero) then
          xl = al
      else
          xh = al
      endif
      if (lit.eq.1) then
         if(dabs(dn).le.0.1d0) then
            al = al - dn
         else
             al = al - 0.1d0*sign(one,dn)
         endif
      else
c           secant method
         if (dd.eq.zero) dd = 1.d-20
         al    = al - dn/dd
         if (al.lt.xl.or.al.gt.xh) then
c              bisection
            al = half*(xl+xh)
            bbb = 'B'
         endif
      endif

      if( dabs(sn-tz(it)).lt.min(0.1d0,0.01d0*si) ) goto 30;
      if(abs(al-alold).lt.epsl) goto 30

      if (lprl.or.lit.gt.10) then
          write(l6,113) lit,'. L-Iteration: ',bbb,alold,dn,al
          write(6,113)  lit,'. L-Iteration: ',bbb,alold,dn,al
          bbb = ' '
      endif
c
c---- end of lambda-loop
      enddo




      write(l6,'(a,i4,a)')
     &     ' Lambda-Iteration interupted after',lit-1,' steps'
      stop
   30 if (lprl) then
         write(l6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
c        write(6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
      endif
      ala(it) = al
#ifdef VERBOSE
      write(6,'(a,i3)')'Number of lambda iterations: ', lit;
#endif

      if (lpr) then
      write(l6,*) ' ****** END DIRHB **********************************'
      endif
c
  101 format(i4,a,i4,3f13.8)
  113 format(i4,a,a1,3f13.8)
c
      return
C-end-DIRHB
      end














