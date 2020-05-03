c======================================================================c

      subroutine base(lpr)

c======================================================================c
c
c     determines the basis in spherical oscillators for Dirac solution 
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      common /basnnn/ n0f,n0b
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)
      common /bosqua/ no
      common /gfviv / iv(0:igfv)
      common /sdimos/ nrm,nlm,nrbm
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN BASE *********************************'
c
      n0fx0 = n0fx
      if (n0f.gt.n0fx0) stop ' in BASE: n0f too large'
      if (is.eq.2.and.n0fx0.lt.20) stop 'in BASE: n0fx < 20'
c
      nrm = 1 + (n0f+1)/2
      nlm = n0f+1
      if (nrm.gt.ndx) stop ' in BASE: ndx too small'
      if (nrm.gt.nrx) stop ' in BASE: nrx too small '
      if (nlm.gt.nlx) stop ' in BASE: nlx too small '
c
c----------------------------------------------------------------------c
c     construction of the different kappa-blocks
c----------------------------------------------------------------------c
      kapmax = n0f + 1
      ib  = 0
      il  = 0
      ik  = 0
      ik1 = 0
      ntz = 0
      nfm = 0
      ngm = 0
c
c
	do jj = 1,n0f+1
		do ipp = 1,2
		ibb = 2*(jj-1) + ipp
		if (mod(jj,2).eq.ipp-1) then
			ik1 = +jj
			l = jj
			la = jj-1
		else
			ik1 = -jj
			l = jj-1
			la = jj
		endif
		ijb(ibb) = jj
		ilb(ibb,1) = l
		ilb(ibb,2) = la
		ipb(ibb) = ipp 
		ikb(ibb) = ik1
		enddo	! ip
	enddo	! j

      do j = 1,kapmax
         do ipar = 1,-1,-2
	        lf = j - (1-ipar*iv(j))/2
            kappa = j * (2*(lf-j)+1)
	        lg    = 2*j - lf -1 

            if (lf.le.n0f) then
               ip = 1 + mod(lf,2)
               nf = (n0f-lf)/2 + 1
               ng = (n0f+1-lg)/2 + 1
               ib = ib + 1
               if (ib.gt.nbx) stop ' in BASE: nbx too small'
               write(tb(ib),'(a1,i2,2h/2)') tl(lf),j+j-1
c
	           kb(ib)    = kappa
	           mb(ib)    = 2*iabs(kappa)
               id(ib,1)  = nf
               id(ib,2)  = ng
	           ia(ib,1)  = il 
	           ia(ib,2)  = il + id(ib,1) 
c
               ntz = ntz + 2*j
c
               do ifg = 1,2
                  ne = id(ib,ifg)
                  l  = lfgkap(kappa,ifg)
                  do n = 1,ne
                     il = il + 1
                     if (il .gt.ntx) stop ' in BASE: ntx too small'
                     nr(il) = n
                     nl(il) = l
                     nj(il) = j
                     kk(il) = kappa
		             np(il) = ifg
                     write(tt(il),'(i2,a1,i3,2h/2)') n,tl(l),j+j-1
                     nn     = 2*(n-1)+l
                  enddo   ! n 
               enddo   ! ifg 
	           ik = ik + max(nf,ng)
	           nfm = max(nfm,nf)
	           ngm = max(ngm,ng)
            endif
c
   10    enddo   ! ivkap
      enddo   ! j
      nb  = ib 
      nt  = il
      nk  = ik
      no  = n0b/2 + 1
      if (nk.gt.nkx)   stop 'in BASE: nkx too small'
      if (nfm.gt.nfx)  stop 'in BASE: nfx too small '
      if (ngm.gt.ngx)  stop 'in BASE: ngx too small '
c
c
c
c
c----------------------------------------------------------------------c
c     Printout
c----------------------------------------------------------------------c
      if (lpr) then
         do i = 1,nt
            nn     =  2*(nr(i)-1) + nl(i)
            write(l6,'(i3,i3,1x,a,i4)') i,nn,tt(i),np(i)
         enddo
         do ib = 1,nb
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng 
            j  = iabs(kb(ib))
            write(l6,'(/,a,4i4)') tb(ib),id(ib,1),id(ib,2),
     &                                   ia(ib,1),ia(ib,2)
            do i = ia(ib,1)+1,ia(ib,1)+nh                             

               nn     =  2*(nr(i)-1) + nl(i)
               write(l6,102) i,'   NN = ',nn,
     &         '   nr = ',nr(i),'   l =',nl(i),'   j =',j+j-1,tt(i)
               if (i.eq.ia(ib,1)+nf) write(l6,'(3x,61(1h-))')

            enddo   ! i       
         enddo   ! ib
c
         write(l6,'(/,a,2i4)') ' Number of blocks: nb  = ',nb,nbx
         write(l6,100) ' Number of levels  nt  = ',nt,ntx 
         write(l6,100) ' Number of levels  nk  = ',nk,nkx 
         write(l6,100) ' Maximal n:        nrm = ',nrm,nrx
         write(l6,100) ' Maximal l:        nlm = ',nlm,nlx
         write(l6,100) ' Maximal nf            = ',nfm,nfx
         write(l6,100) ' Maximal ng            = ',ngm,ngx
      endif
c    
c
c----------------------------------------------------------------------c
c     Construction of two-body pairs (i1,i2)                      
c----------------------------------------------------------------------c
c---- only f-pairs
      il = 0
      do ib = 1,nb
         ipos(ib) = il
         nf = id(ib,1)
         do n2 =  1,nf
         do n1 = n2,nf
            il = il + 1
            nib(il)   = ib
            nni(1,il) = n1
            nni(2,il) = n2
         enddo  ! n1
         enddo  ! n2
      enddo  ! ib
      mv = il
      if (lpr) write(l6,100) ' number of pairs f  mv = ',mv,mvx
      if (mv.gt.mvx) then
	 write(6,*) 'mv =',mv,' mvx = ',mvx
	 stop ' in BASE: mvx too small'
      endif
c
  100 format(a,2i6)
  102 format(i4,a,i2,a,i2,a,i2,a,i2,3h/2 ,a) 
c
      if (lpr)
     &write(l6,*) ' ****** END BASE ***********************************'
c 
c-end-BASE
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
      parameter (nn = 2*(MVX+MVXG)+2*MVX+2)
      parameter (mm = 7)      
c
      common /baspar/ hom,hb0,b0
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx) 
      common /physco/ hbc,alphi,r0
      common /deldel/ de(nhhx,nb2x) 
      common /gamgam/ hh(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
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
                  if(del(it).lt.1.d-5) then 
                     vou(ipos)=zero
                     vin(ipos)=zero
                  endif
               enddo
            enddo            
         enddo
         ipos=ipos+1
         vou(ipos)=ala(it)-vin(ipos)
      enddo   
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
c
c
      if (lpr) then
      write(l6,*) '****** END BROYDE **************************'
      endif
c
      return
c-end-BROYDEN
      end

      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTRF computes the factorization of a real symmetric matrix A using
*  the Bunch-Kaufman diagonal pivoting method.  The form of the
*  factorization is
*
*     A = U*D*U**T  or  A = L*D*L**T
*
*  where U (or L) is a product of permutation and unit upper (lower)
*  triangular matrices, and D is symmetric and block diagonal with
*  1-by-1 and 2-by-2 diagonal blocks.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, the block diagonal matrix D and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK >=1.  For best performance
*          LWORK >= N*NB, where NB is the block size returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
*                has been completed, but the block diagonal matrix D is
*                exactly singular, and division by zero will occur if it
*                is used to solve a system of equations.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', then A = U*D*U', where
*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    v    0   )   k-s
*     U(k) =  (   0    I    0   )   s
*             (   0    0    I   )   n-k
*                k-s   s   n-k
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*
*  If UPLO = 'L', then A = L*D*L', where
*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    0     0   )  k-1
*     L(k) =  (   0    I     0   )  s
*             (   0    v     I   )  n-k-s+1
*                k-1   s  n-k-s+1
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASYF, DSYTF2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )
     $   NB = N
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U' using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        KB, where KB is the number of columns factorized by DLASYF;
*        KB is either NB or NB-1, or K for the last block
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 40
*
         IF( K.GT.NB ) THEN
*
*           Factorize columns k-kb+1:k of A and use blocked code to
*           update columns 1:k-kb
*
            CALL DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK,
     $                   IINFO )
         ELSE
*
*           Use unblocked code to factorize columns 1:k of A
*
            CALL DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*        Decrease K and return to the start of the main loop
*
         K = K - KB
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L' using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        KB, where KB is the number of columns factorized by DLASYF;
*        KB is either NB or NB-1, or N-K+1 for the last block
*
         K = 1
   20    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( K.LE.N-NB ) THEN
*
*           Factorize columns k:k+kb-1 of A and use blocked code to
*           update columns k+kb:n
*
            CALL DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ),
     $                   WORK, LDWORK, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns k:n of A
*
            CALL DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + K - 1
*
*        Adjust IPIV
*
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSE
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
*
*        Increase K and return to the start of the main loop
*
         K = K + KB
         GO TO 20
*
      END IF
*
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DSYTRF
*
      END
      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTRI computes the inverse of a real symmetric indefinite matrix
*  A using the factorization A = U*D*U**T or A = L*D*L**T computed by
*  DSYTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the details of the factorization are stored
*          as an upper or lower triangular matrix.
*          = 'U':  Upper triangular, form is A = U*D*U**T;
*          = 'L':  Lower triangular, form is A = L*D*L**T.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the block diagonal matrix D and the multipliers
*          used to obtain the factor U or L as computed by DSYTRF.
*
*          On exit, if INFO = 0, the (symmetric) inverse of the original
*          matrix.  If UPLO = 'U', the upper triangular part of the
*          inverse is formed and the part of A below the diagonal is not
*          referenced; if UPLO = 'L' the lower triangular part of the
*          inverse is formed and the part of A above the diagonal is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D
*          as determined by DSYTRF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*               inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KP, KSTEP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSWAP, DSYMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U'.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ),
     $                     1 )
               A( K, K+1 ) = A( K, K+1 ) -
     $                       DDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL DCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) -
     $                         DDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            CALL DSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL DSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L'.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   50    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 60
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
               A( K, K-1 ) = A( K, K-1 ) -
     $                       DDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ),
     $                       1 )
               CALL DCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) -
     $                         DDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            IF( KP.LT.N )
     $         CALL DSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL DSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of DSYTRI
*
      END
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DX(I)
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     +            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
**
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  ILAENV returns an INTEGER
*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      SUBROUTINE DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASYF computes a partial factorization of a real symmetric matrix A
*  using the Bunch-Kaufman diagonal pivoting method. The partial
*  factorization has the form:
*
*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
*        ( 0  U22 ) (  0   D  ) ( U12' U22' )
*
*  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
*        ( L21  I ) (  0  A22 ) (  0    I   )
*
*  where the order of D is at most NB. The actual order is returned in
*  the argument KB, and is either NB or NB-1, or N if N <= NB.
*
*  DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
*  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
*  A22 (if UPLO = 'L').
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NB      (input) INTEGER
*          The maximum number of columns of the matrix A that should be
*          factored.  NB should be at least 2 to allow for 2-by-2 pivot
*          blocks.
*
*  KB      (output) INTEGER
*          The number of columns of A that were actually factored.
*          KB is either NB-1 or NB, or N if N <= NB.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, A contains details of the partial factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If UPLO = 'U', only the last KB elements of IPIV are set;
*          if UPLO = 'L', only the first KB elements are set.
*
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W.  LDW >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*               has been completed, but the block diagonal matrix D is
*               exactly singular.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     $                   KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1,
     $                   ROWMAX, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Factorize the trailing columns of A using the upper triangle
*        of A and working backwards, and compute the matrix W = U12*D
*        for use in updating A11
*
*        K is the main loop index, decreasing from N in steps of 1 or 2
*
*        KW is the column of W which corresponds to column K of A
*
         K = N
   10    CONTINUE
         KW = NB + K - N
*
*        Exit from loop
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
*
*        Copy column K of A to column KW of W and update it
*
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA,
     $                  W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( W( K, KW ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column KW-1 of W and update it
*
               CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                     W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N )
     $            CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),
     $                        LDA, W( IMAX, KW+1 ), LDW, ONE,
     $                        W( 1, KW-1 ), 1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = ABS( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, KW-1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column KW-1 of W to column KW
*
                  CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            KKW = NB + KK - N
*
*           Updated column KP is already stored in column KKW of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
*
*              Interchange rows KK and KP in last KK columns of A and W
*
               CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column KW of W now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Store U(k) in column k of A
*
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-by-2 pivot block D(k): columns KW and KW-1 of W now
*              hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
               IF( K.GT.2 ) THEN
*
*                 Store U(k) and U(k-1) in columns k and k-1 of A
*
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
   30    CONTINUE
*
*        Update the upper triangle of A11 (= A(1:k,1:k)) as
*
*        A11 := A11 - U12*D*U12' = A11 - U12*W'
*
*        computing blocks of NB columns at a time
*
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
*
*           Update the upper triangle of the diagonal block
*
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', JJ-J+1, N-K, -ONE,
     $                     A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,
     $                     A( J, JJ ), 1 )
   40       CONTINUE
*
*           Update the rectangular superdiagonal block
*
            CALL DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE,
     $                  A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE,
     $                  A( 1, J ), LDA )
   50    CONTINUE
*
*        Put U12 in standard form by partially undoing the interchanges
*        in columns k+1:n
*
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF( JP.NE.JJ .AND. J.LE.N )
     $      CALL DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LE.N )
     $      GO TO 60
*
*        Set KB to the number of columns factorized
*
         KB = N - K
*
      ELSE
*
*        Factorize the leading columns of A using the lower triangle
*        of A and working forwards, and compute the matrix W = L21*D
*        for use in updating A22
*
*        K is the main loop index, increasing from 1 in steps of 1 or 2
*
         K = 1
   70    CONTINUE
*
*        Exit from loop
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
*
*        Copy column K of A to column K of W and update it
*
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA,
     $               W( K, 1 ), LDW, ONE, W( K, K ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( W( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column K+1 of W and update it
*
               CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ),
     $                     1 )
               CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),
     $                     LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = ABS( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, K+1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column K+1 of W to column K
*
                  CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
*
*           Updated column KP is already stored in column KK of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
*
*              Interchange rows KK and KP in first KK columns of A and W
*
               CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k of W now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
*              Store L(k) in column k of A
*
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = ONE / A( K, K )
                  CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k+1 of W now hold
*
*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
*
*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
*              of L
*
               IF( K.LT.N-1 ) THEN
*
*                 Store L(k) and L(k+1) in columns k and k+1 of A
*
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 70
*
   90    CONTINUE
*
*        Update the lower triangle of A22 (= A(k:n,k:n)) as
*
*        A22 := A22 - L21*D*L21' = A22 - L21*W'
*
*        computing blocks of NB columns at a time
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           Update the lower triangle of the diagonal block
*
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE,
     $                     A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,
     $                     A( JJ, JJ ), 1 )
  100       CONTINUE
*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                     K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW,
     $                     ONE, A( J+JB, J ), LDA )
  110    CONTINUE
*
*        Put L21 in standard form by partially undoing the interchanges
*        in columns 1:k-1
*
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF( JP.NE.JJ .AND. J.GE.1 )
     $      CALL DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GE.1 )
     $      GO TO 120
*
*        Set KB to the number of columns factorized
*
         KB = K - 1
*
      END IF
      RETURN
*
*     End of DLASYF
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for xHSEQR and its subroutines. It is called whenever 
*       ILAENV is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ISPEC  (input) integer scalar
*              ISPEC specifies which tunable parameter IPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to xLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        IPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents TTQRE from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  During the multi-shift QR sweep,
*                            xLAQR5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  During the multi-shift QR sweep,
*                            xLAQR5 and/or xLAQRaccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            xLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (If xTRMM is slower than xGEMM, then
*                        IPARMQ(ISPEC=16)=1 may be more efficient than
*                        IPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       NAME    (input) character string
*               Name of the calling subroutine
*
*       OPTS    (input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (input) INTEGER
*       IHI     (input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORK   (input) integer scalar
*               The amount of workspace available.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xLAQR3 and xLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by IPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*                        Default: 75. (Must be at least 11.)
*
*       IPARMQ(ISPEC=13) Recommended deflation window size.
*                        This depends on ILO, IHI and NS, the
*                        number of simultaneous shifts returned
*                        by IPARMQ(ISPEC=15).  The default for
*                        (IHI-ILO+1).LE.500 is NS.  The default
*                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*
*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                        a multi-shift QR iteration.
*
*                        If IHI-ILO+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       NS =   2+
*                               30               60       NS =   4+
*                               60              150       NS =  10
*                              150              590       NS =  **
*                              590             3000       NS =  64
*                             3000             6000       NS = 128
*                             6000             infinity   NS = 256
*
*                    (+)  By default matrices of this order are
*                         passed to the implicit double shift routine
*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*                         values of NS are used only in case of a rare
*                         xLAHQR failure.
*
*                    (**) The asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       IPARMQ(ISPEC=16) Select structured matrix multiply.
*                        (See ISPEC=16 above for details.)
*                        Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         IPARMQ = 0
         IF( NS.GE.KACMIN )
     $      IPARMQ = 1
         IF( NS.GE.K22MIN )
     $      IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      IF (B(L,J).NE.ZERO) THEN
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      END IF
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B' + beta*C
*
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*B(J,L)
                          DO 150 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  150                     CONTINUE
                      END IF
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     uses unrolled loops for increments equal one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  END IF
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  DSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END
      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTF2 computes the factorization of a real symmetric matrix A using
*  the Bunch-Kaufman diagonal pivoting method:
*
*     A = U*D*U'  or  A = L*D*L'
*
*  where U (or L) is a product of permutation and unit upper (lower)
*  triangular matrices, U' is the transpose of U, and D is symmetric and
*  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, the block diagonal matrix D and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*               has been completed, but the block diagonal matrix D is
*               exactly singular, and division by zero will occur if it
*               is used to solve a system of equations.
*
*  Further Details
*  ===============
*
*  09-29-06 - patch from
*    Bobby Cheng, MathWorks
*
*    Replace l.204 and l.372
*         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*    by
*         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
*
*  01-01-96 - Based on modifications by
*    J. Lewis, Boeing Computer Services Company
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
*         Company
*
*  If UPLO = 'U', then A = U*D*U', where
*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    v    0   )   k-s
*     U(k) =  (   0    I    0   )   s
*             (   0    0    I   )   n-k
*                k-s   s   n-k
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*
*  If UPLO = 'L', then A = L*D*L', where
*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    0     0   )  k-1
*     L(k) =  (   0    I     0   )  s
*             (   0    v     I   )  n-k-s+1
*                k-1   s  n-k-s+1
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1,
     $                   ROWMAX, T, WK, WKM1, WKP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTF2', -INFO )
         RETURN
      END IF
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U' using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 70
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
*
*           Column K is zero or contains a NaN: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the leading
*              submatrix A(1:k,1:k)
*
               CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           Update the leading submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Perform a rank-1 update of A(1:k-1,1:k-1) as
*
*              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
*
               R1 = ONE / A( K, K )
               CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
*
*              Store U(k) in column k
*
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k-1 now hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
*              Perform a rank-2 update of A(1:k-2,1:k-2) as
*
*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
*
               IF( K.GT.2 ) THEN
*
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
*
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*WK -
     $                              A( I, K-1 )*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
*
               END IF
*
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L' using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2
*
         K = 1
   40    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 70
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
*
*           Column K is zero or contains a NaN: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the trailing
*              submatrix A(k:n,k:n)
*
               IF( KP.LT.N )
     $            CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           Update the trailing submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
               IF( K.LT.N ) THEN
*
*                 Perform a rank-1 update of A(k+1:n,k+1:n) as
*
*                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
*
                  D11 = ONE / A( K, K )
                  CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1,
     $                       A( K+1, K+1 ), LDA )
*
*                 Store L(k) in column K
*
                  CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-by-2 pivot block D(k)
*
               IF( K.LT.N-1 ) THEN
*
*                 Perform a rank-2 update of A(k+2:n,k+2:n) as
*
*                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
*
*                 where L(k) and L(k+1) are the k-th and (k+1)-th
*                 columns of L
*
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
*
                  DO 60 J = K + 2, N
*
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
*
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*WK -
     $                              A( I, K+1 )*WKP1
   50                CONTINUE
*
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
*
   60             CONTINUE
               END IF
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 40
*
      END IF
*
   70 CONTINUE
*
      RETURN
*
*     End of DSYTF2
*
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  Purpose
*  =======
*
*  DSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR  .
*
      END
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      LOGICAL FUNCTION DISNAN(DIN)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN      (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
      LOGICAL FUNCTION DLAISNAN(DIN1,DIN2)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN1,DIN2
*     ..
*
*  Purpose
*  =======
*
*  This routine is not for general use.  It exists solely to avoid
*  over-optimization in DISNAN.
*
*  DLAISNAN checks for NaNs by comparing its two arguments for
*  inequality.  NaN is the only floating-point value where NaN != NaN
*  returns .TRUE.  To check for NaNs, pass the same variable as both
*  arguments.
*
*  Strictly speaking, Fortran does not allow aliasing of function
*  arguments. So a compiler must assume that the two arguments are
*  not the same variable, and the test will not be optimized away.
*  Interprocedural or whole-program optimization may delete this
*  test.  The ISNAN functions will be replaced by the correct
*  Fortran 03 intrinsic once the intrinsic is widely available.
*
*  Arguments
*  =========
*
*  DIN1     (input) DOUBLE PRECISION
*  DIN2     (input) DOUBLE PRECISION
*          Two numbers to compare for inequality.
*
*  =====================================================================
*
*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) INTEGER 
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END

c=====================================================================c
  
      subroutine canon(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpr

      character*8 tbb(nhbx)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tk*8,tx*8                                       ! blolev
      character tt*8                                            ! bloqua
      character tb*5                                            ! blokap
c
c
      dimension aa(nhhx),dd(nhhx),v2(nhx),z(nhx),eb(nhx),h(nhx),d(nhx)
      
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)      
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /gfvsq / sq(0:igfv)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /gamgam/ hh(nhhx,nb2x)
      common /deldel/ de(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(nhbx,nkx,4),equ(nkx,4)
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /baspar/ hom,hb0,b0
      common /wavec/ wave(ngh2,nkx,2),dwave(ngh2,nkx,2)
      
c 
      data ash/100.d0/ 
      if(lpr)
     &write(l6,*) ' ****** BEGIN CANON ********************************'
c
c
c---- BCS structure in the canonical basis:
c
      do it = 1,2
         write(l6,100) tit(it)
  100    format(' single-particle energies and gaps ',1x,
     &          'in the canonical basis: ',a,/,1x,66(1h-))
         write(l6,101) 'j pi','(nr,l,j)','smax','Ecan',
     &                  'vvcan','decan'
  101    format(4x,a,5x,a,2x,a,7x,a,7x,a,7x,a)
c
c---- loop of the j-blocks

         klp = 0
         kla = 0
         do ib = 1,nb
            kap = kb(ib)
            j   = iabs(kap)
	        mul = 2*j
            lf  = lfkap(kap)
            lg  = lgkap(kap)
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            nhfb= nh + nh
            i0f = ia(ib,1)
            i0g = ia(ib,2)
	        mf  = ib + (it-1)*nbx
            kf  = kd(ib,it)
            kg  = kd(ib,it+2)
	        k0f = ka(ib,it)
	        k0g = ka(ib,it+2)

c------- transformation to the canonical basis
c------- calculation of the generalized density V*VT     
            do n2 = 1,nh
               do n1 = 1,nh
                  s = zero
                  do k = k0f+1,k0f+kf
                     s = s + fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
                  enddo				  
                  s1 = zero
                  do k = k0g+1,k0g+kg
                     s1 = s1 + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
                  enddo       
                  aa(n1+(n2-1)*nh) = s + ash*s1
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
               do n2 = 1,nh
                  h2 = zero
                  d2 = zero
                  do n1 = 1,nh
                     h2 = h2 + dd(n1+(k-1)*nh)*hh(n1+(n2-1)*nh,mf)
                     d2 = d2 + dd(n1+(k-1)*nh)*de(n1+(n2-1)*nh,mf)
                  enddo
                  hk = hk + h2*dd(n2+(k-1)*nh) 
                  dk = dk + d2*dd(n2+(k-1)*nh) 
               enddo
               h(k) = hk
               d(k) = dk
            enddo   ! k

c------- reordering according to the energy h(k) 
            call ordx(nh,h,d,v2,dd)
                        
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
               mulcan(kla,it+2)=mul
               do i=1,nh
                  fgcan(i,kla,it+2)=dd(i+(k-1)*nh)
               enddo
            enddo
            kdcan(ib,it+2) = kla - kacan(ib,it+2)
            
c------- wave functions in the coordinate space    
            k1=kacan(ib,it)+1
            k2=kacan(ib,it)+kdcan(ib,it)
            do k = k1,k2
               do ih=1,ngh
                  wx = xh(ih)*sqrt(wh(ih))*b0**(3./2.)*sqrt(4*pi)
                  wxd= wx*b0
                  s = 0.0
                  sd = 0.0
                  do j = 1, nf    !loop over oscilator states
                     s = s + fgcan(j,k,it)*rnl(j,lf,ih)
                     sd = sd + fgcan(j,k,it)*rnl1(j,lf,ih)
                  enddo                                     
                  wave(ih,k,it) = s/wx
                  dwave(ih,k,it) = sd/wxd
                  s = 0.0
                  sd = 0.0
                  do j = 1, ng
                     s = s + fgcan(nf+j,k,it)*rnl(j,lg,ih)
                     sd = sd + fgcan(nf+j,k,it)*rnl1(j,lg,ih)
                  enddo
                  wave(ngh+ih,k,it) = s/wx
                  dwave(ngh+ih,k,it) = sd/wxd 
               enddo  
               suml=zero
               sums=zero
               do ih=1,ngh
                  wx = wdcor(ih)
                  suml=suml+wave(ih,k,it)**2*wx
                  sums=sums+wave(ngh+ih,k,it)**2*wx 
               enddo
               snorm = one/sqrt(suml+sums)
               do ih=1,ngh2
                  wave(ih,k,it) = wave(ih,k,it)*snorm
                  dwave(ih,k,it) = dwave(ih,k,it)*snorm
               enddo
            enddo
c------- printout for particles
            if (ib.eq.1) e0 = h(ng+1)
            k1 = kacan(ib,it)+1
            k2 = kacan(ib,it)+kdcan(ib,it)
            do k = k1,k2
               e1 = eecan(k,it)
               v1 = vvcan(k,it)
               d1 = decan(k,it)
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
               ipar = ipb(ib)
c
               if (lpr) write(l6,103) 2*ijb(ib)-1,tp(ipar),
     &                    tb(ib),dx,e1,v1,d1
  103          format(' j =',i2,'/2',a1,'  (',a8,')',
     &                f6.2,4f12.4)
           enddo  ! k        
c---- end loop over j-blocks
   10    enddo !ib
         nkcan(it)   = klp
         nkcan(it+2) = kla
c------- printout for anti-particles   
c         if (lpr) write(l6,103) tit(it)
c  103    format(' anti-particle energies and gaps ',1x,
c     &          'in the canonical basis: ',a,/,1x,66(1h-))
c         do ib=1,nb
c            kap = kb(ib)
c            j   = iabs(kap)
c	        mul = 2*j
c            lf  = lfkap(kap)
c            lg  = lgkap(kap)
c            nf  = id(ib,1)
c            ng  = id(ib,2)
c            nh  = nf + ng
c            nhfb= nh + nh
c            i0f = ia(ib,1)
c            i0g = ia(ib,2)
c	        mf  = ib + (it-1)*nbx
c            kf  = kd(ib,it)
c            kg  = kd(ib,it+2)
c	        k0f = ka(ib,it)
c	        k0g = ka(ib,it+2)
c            k1 = kacan(ib,it+2)+1
c            k2 = kacan(ib,it+2)+kdcan(ib,it+2)
c            do k = k1,k2
c               e1 = eecan(k,it+2)
c               v1 = vvcan(k,it+2)
c               d1 = decan(k,it+2)
c               ipar = ipb(ib)
c
c               if (lpr) write(l6,104) 2*ijb(ib)-1,tp(ipar),
c     &                    tb(ib),e1,v1,d1
c  104          format(' j =',i2,'/2',a1,'  (',a8,')',
c     &                f9.2,2f10.3)
c           enddo  ! k        	        
c         enddo
         
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
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
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
c           write(6,'(3i3,a,/,(10f12.8))') mm,ma,me,
c     &      ' Eigenvalues degenerate:',(ea(i),i=ma,me)
c          write(l6,'(3i3,a,/,(10f12.8))') mm,ma,me,
c    &      ' Eigenvalues degenerate:',(ea(i),i=ma,me)

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
c           call aprint(1,1,6,1,1,mm,eb,' ',' ','H_x')
c	    call aprint(1,1,6,1,1,mm,eb(ma),' ',' ','H_x')
c           call aprint(1,1,1,na,mm,mm,zz,' ',' ','DD_x')
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
       subroutine centmas(lpr)
c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
      logical lpr
      
      call cmcd(lpr)
      call cmcn(lpr)
      
      return
      end

c=====================================================================c 
       subroutine cmcd(lpr)
c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpr
      character tb*5                                            ! blokap
      
      dimension l(2)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)   
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)
      common /masses/ amu,amsig,amome,amdel,amrho
      common /physco/ hbc,alphi,r0
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4)
      common /wavec/ wave(ngh2,nkx,2),dwave(ngh2,nkx,2)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
      
      if (icm.ne.2) return
      
      fac = one/(two*amu*amas)
      fac = fac*hbc
      
      do it=1,itx
         partd = zero
         do ib=1,nb
            kap = kb(ib)
            j   = iabs(kap)
	        mul = 2*j
            l(1)= lfkap(kap)
            l(2)= lgkap(kap)
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            k1  = kacan(ib,it)+1
            k2  = kacan(ib,it)+kdcan(ib,it)
            partds = zero
            do k=k1,k2
               vv = vvcan(k,it)*mul
               partdl = zero
               do ifg=1,2
                  ll = l(ifg)
                  comp1 = zero
                  comp2 = zero
                  do ih=1,ngh
                     wx = wdcor(ih)
                     wx1 = wx/rb(ih)**2
                     ihg = (ifg-1)*ngh+ih
                     comp1 = comp1 + dwave(ihg,k,it)**2*wx
                     comp2 = comp2 + wave(ihg,k,it)**2*wx1
                  enddo
                  partdl = partdl+comp1+ll*(ll+1)*comp2
               enddo
               partds = partds + vv*partdl              
            enddo !k
            partd = partd+partds
         enddo !ib
         ecmd(it) = -partd*fac
      enddo !it
      ecmd(3) = ecmd(1)+ecmd(2)
      return
      end
      
c=====================================================================c
  
       subroutine cmcn(lpr)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'

      logical lpr
      character tb*5                                            ! blokap
      
      dimension la(2),lb(2)
      common /baspar/ hom,hb0,b0
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)   
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)
      common /gfviv / iv(0:igfv)
      common /masses/ amu,amsig,amome,amdel,amrho
      common /physco/ hbc,alphi,r0
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /eeecan/ eecan(nkx,4),decan(nkx,4),vvcan(nkx,4),
     &                fgcan(nhx,nkx,4),mulcan(nkx,4)
      common /wavec/ wave(ngh2,nkx,2),dwave(ngh2,nkx,2)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
      
      if(icm .lt. 2) then
         ecm(1) = -0.75d0*hom
         ecm(2) = -0.75d0*hom
         ecm(3) = -0.75d0*hom
         return
      endif

      fac = one/(two*amu*amas)
      fac = fac*hbc
      
      do it=1,itx
         partn = zero
         do iba=1,nb
            kapa = kb(iba)
            ja   = iabs(kapa)
	        mula = 2*ja
            la(1)= lfkap(kapa)
            la(2)= lgkap(kapa)
            nfa  = id(iba,1)
            nga  = id(iba,2)
            nha  = nfa + nga
            k1a  = kacan(iba,it)+1
            k2a  = kacan(iba,it)+kdcan(iba,it)   
            do ibb = 1,nb
               kapb = kb(ibb)
               jb   = iabs(kapb)
	           mulb = 2*jb
               lb(1)= lfkap(kapb)
               lb(2)= lgkap(kapb)
               nfb  = id(ibb,1)
               ngb  = id(ibb,2)
               nhb  = nfb + ngb
               k1b  = kacan(ibb,it)+1
               k2b  = kacan(ibb,it)+kdcan(ibb,it) 
               do kka=k1a,k2a
                  do kkb=k1b,k2b
                     va = sqrt(vvcan(kka,it))
                     vb = sqrt(vvcan(kkb,it))
                     ua = sqrt(one-va**2)
                     ub = sqrt(one-vb**2)
                     fac2 = va*vb*(va*vb+ua*ub)
                     sum = zero
                     do ifg=1,2
                        lla=la(ifg)
                        llb=lb(ifg)
                        if(abs(lla-llb).eq.1) then
                           if( (ja-lla).eq.(jb-llb)) then
                              f1 = (lla+jb+two)*(lla+jb-one)
                              f2 = (2*lla+one)*(2*llb+one)
                              sixfcu = sqrt(f1/f2)
                           else
                              f1 = (jb-lla+one)*(lla-jb+two)
                              f2 = (2*lla+one)*(2*llb+one)
                              sixfcu = (2*jb-2*llb-1)*sqrt(f1/f2)                              
                           endif
                           if(lla .gt. llb) then
                              facdru = sqrt(dble(lla))
                              facctu = -dble(lla)
                           else
                              facdru = -sqrt(dble(llb))
                              facctu = dble(llb)
                           endif
                           rint = zero
                           do ih=1,ngh
                              ihg = (ifg-1)*ngh+ih
                              wx = wdcor(ih)
                              rrb = wave(ihg,kkb,it)
                              f1  = wave(ihg,kka,it)
                              f2  = dwave(ihg,kka,it)
                              rra = f2-(facctu-one)*f1/rb(ih)
                              rint = rint + rrb*rra*wx
                           enddo
                           sum = sum +sixfcu*facdru*rint
                        endif
                     enddo !ifg
                     partn = partn + fac2*sum**2
                  enddo ! kkb
               enddo !kka
            enddo !ibb
         enddo !iba
         ecmn(it) = fac*partn
         ecm(it) = ecmd(it)+ecmn(it)
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
      common /baspar/ hom,hb0,b0
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /procou/ ggc(0:ngh,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (icou.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN COULOM *******************************'
c
      do i = 0,ngh
         s = zero
         do k = 0,ngh
            s = s + ggc(i,k)*drvp(k)
         enddo
         cou(i) = s*b0
      enddo
c
c      if (lpr) call prigh(1,cou,one,'Coulom')
c
      if (lpr)
     &write(l6,*) ' ****** END COULOM *********************************'
      return
C-end-COULOM
      end
c======================================================================c

      subroutine greecou(lpr)

c======================================================================c
C
C     calculation of the meson and Coulomb-propagator
c
c     imes = 1:   sigma
c            2:   omega
c            3:   delta
c            4:   rho
c            0:   coulomb
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
c
      dimension gi(nox,nox),rr(nox,ngh)
c
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /procou/ ggc(0:ngh,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (icou.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GREECOU ******************************'
c
c
c---- Coulomb-propagator
      f0 = one/(6*alphi)
      do kh = 0,ngh
         f  = f0*wdcor(kh)
         do ih = 0,ngh
            r1 = xh(ih)
            r2 = xh(kh)
            rg = dmax1(r1,r2)
            rk = dmin1(r1,r2)
            ggc(ih,kh) =  f * ( 3*rg + rk**2/rg)
         enddo   ! ih
      enddo    ! ik
c      if (lpr) call aprint(1,1,6,ngh,ngh,ngh,ggc,' ',' ','VC')
c
      if (lpr)
     &write(l6,*) ' ****** END GREECOU ********************************'
c
      return
c-end-GREECOU
      end

c======================================================================c

      subroutine default()

c======================================================================c
c
c     Default for Relativistic Mean Field spherical
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tp*1,tis*1,tit*8,tl*1                   ! textex

      common /baspar/ hom,hb0,b0
      common /couplf/ ff(0:ngh,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /broyde2/ ibroyd
c

c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutron:','Proton: '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
c
c---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/
c======================================================================c
c---- signs and factorials
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



c======================================================================c
c     Coulomb-Field
c-----------------------------------------------------------------------
c     icou: Coulomb-field:  0   not at all 
c                           1   only direct term  
c                           2   plus exchange 
      icou = 1
c======================================================================c

C======================================================================c
c     pairing
c----------------------------------------------------------------------c
      do it = 1,2                   
	     spk0(it) = zero
	     del(it)  = zero
	     ala(it)  = -7.0
      enddo   ! it
c======================================================================c
c     iteration
c----------------------------------------------------------------------c
c
      maxi = 500               ! maximal number of iteration
      si   = one             ! actual error in the main iteration
      epsi = 1.d-6       ! accuracy for the main iteration
      iaut  = 1              ! automatic change of xmix: 0 (no) 1 (yes)
      inxt  = 3              ! number of iterations until next question
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
c     inink = 0: pairing potential read, 1: pairing potential monopol
      inink = 1
c
c     oscillator length b0 (is calcuated for b0 <= 0)
      b0 = -2.320
c
c


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
      lvpp = 13
c======================================================================c
c---- preparation of density dependence
c----------------------------------------------------------------------c
      do m = 1,4
         lmes(m) = 0
         ggmes(m) = zero
         do i = 0,ngh
            ff(i,m,1) = one
	    ff(i,m,2) = zero
         enddo   ! ngh
      enddo   ! m

      return
c-end-DEFAULT
      end
c======================================================================c

      blockdata block1

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tp*1,tis*1,tit*8,tl*1                   ! textex

      common /physco/ hbc,alphi,r0
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c

c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutron:','Proton: '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
c
c---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/

c-end-BLOCK1      
      end
C======================================================================c

      subroutine delta(it,lpr)

c======================================================================c
c
c     calculats the pairing field for ide 3,4 ...                       
c
c     ide = 0:     no pairing field
c           1:     constant gap
c           2:     constant G
c           3:     GOGNY
c 
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical  lpr
c
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
      dimension pn(0:n0fx)
c
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)
      common /vvvsep/ vn(mvx,0:n0fx)
      common /tmrpar/  gl(2),gal

      if (lpr)
     &write(l6,*) ' ****** BEGIN DELTA ********************************'
c
c----------------------------------------------------------------------c
c        Separable pairing
c----------------------------------------------------------------------c
	  do il=0,nnmax
		 pn(il) = 0.d0
		 do i12=1,mv
		    pn(il) = pn(il) + vn(i12,il)*aka(i12,it)
		 enddo
	  enddo	
      i34 = 0
	  do ib = 1,nb
	     i0 = ia(ib,1)
	     nf = id(ib,1)
	     ng = id(ib,2)
	     nh = nf + ng
	     m  = ib + (it-1)*nbx
         do n2=1,nf
            do n1=n2,nf
               i34=i34+1
               s= zero
               do il=0,nnmax
                  s = s+pn(il)*vn(i34,il)
               enddo
               de(n1+(n2-1)*nh,m)=-half*s*gl(it)
               de(n2+(n1-1)*nh,m)=-half*s*gl(it)
            enddo
         enddo   
	  enddo	! ib
c----------------------------------------------------------------------c
c     print out                   
c----------------------------------------------------------------------c
      if (lpr) then
         do ib = 1,nb
	        nf = id(ib,1) 
	        ng = id(ib,2) 
            nh = nf + ng
	        k0 = ia(ib,1)+1
	        m  = ib + (it-1)*nbx
c	        call aprint(1,3,1,nh,nf,nf,de(1,m),tt(k0),tt(k0),'DE++')
	     enddo   ! ib
         write(l6,*) ' ****** END DELTA *****************************'
      endif
c
      return
C-end-DELTA
      end

C=====================================================================c

      subroutine densit(lpr)

c=====================================================================c
C
C     density at the radius r = xh(ih)*b0 is given by
C     b0**(-3) * rv(ih) / ( 4*pi * r**2 * wh) in units of fm**(-3)
C
c---------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
c
      dimension lfg(2),rsfg(2)
      dimension drs(0:ngh,2),drv(0:ngh,2)
c
      common /baspar/ hom,hb0,b0
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /gfviv / iv(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      data eps/1.d-7/
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DENSIT *******************************'
c
      do it = 1,itx
         do ih = 0,ngh
            rs(ih,it)  = zero
            rv(ih,it)  = zero
            drs(ih,it) = zero
            drv(ih,it) = zero
         enddo
      enddo
c
c     loop over j-blocks
      il = 0
      do ib = 1,nb
         kappa = kb(ib)
	     nf = id(ib,1)
	     ng = id(ib,2)
	     nh = nf + ng
c---- loop over contributions from large and small components
         do ifg = 1,2
	        n0  = (ifg-1)*nf
	        nd  = id(ib,ifg)
	        i0  = ia(ib,ifg)
	        ivv = iv(ifg)
c
c---------- loop over oscillator basis sttes n2 and n1
            do n2 =  1,nd
            do n1 = n2,nd
               i12 = 2 - n2/n1
               do it = 1,itx
                  rsfg(it) = i12*rosh(n0+n1+(n0+n2-1)*nh,ib+(it-1)*nbx)
               enddo   
c
	           l  = lfgkap(kappa,ifg)
               ll = l*(l+1)
               nn = 2*(n1+n2+l)-1
c
c------------- loop over the meshpoints
               do ih = 0,ngh
                  s  = rnl(n1,l,ih)*rnl(n2,l,ih)
                  s1 = rnl1(n1,l,ih)*rnl1(n2,l,ih)
                  xx = xh(ih)*xh(ih)
		          ds = 2*(s*(xx+ll/xx-nn)+s1)
                  do it = 1,itx
		             fgr  = rsfg(it)*s
		             fgd  = rsfg(it)*ds
                     rs(ih,it)  = rs(ih,it)  - ivv*fgr
                     rv(ih,it)  = rv(ih,it)  + fgr                   
c                    Delta-rho 
                     drv(ih,it) = drv(ih,it) + fgd 
                     drs(ih,it) = drs(ih,it) - ivv*fgd 
                  enddo   ! it
               enddo   ! ih
   10       enddo   ! n1
            enddo   ! n2
         enddo   !   ifg  (large and small components)
      enddo   ! ib (loop over the blocks)
      
c---- check, whether integral over drv vanishes
      s = zero
      do it = 1,itx
         do ih = 1,ngh
            s = s + drv(ih,it)
         enddo
         if (lpr) write(l6,*) 'integral over dro',it,s
      enddo
c
c
c---- normalization and renormalization to particle number
      do it = 1,itx
         s  = zero
         sp = zero
         do ih = 1,ngh
            s  = s  + rv(ih,it)
         enddo
         if (lpr) then
	    write(6,100) ' norm of the vector density = ',it,s
         endif
         do ih = 0,ngh
            f           = one/wdcor(ih)
            rs(ih,it)   = f*rs(ih,it)
            rv(ih,it)   = f*rv(ih,it)
            drs(ih,it)  = f*drs(ih,it)/b0**2
            drv(ih,it)  = f*drv(ih,it)/b0**2
         enddo

      enddo   ! it
      if (itx.eq.1) then
         do ih = 0,ngh
            ro(ih,1)  = 2 * rs(ih,1)
            ro(ih,2)  = 2 * rv(ih,1)
            ro(ih,3)  = zero
            ro(ih,4)  = zero
	        dro(ih,1) = 2 * drs(ih,1) 
	        dro(ih,2) = 2 * drv(ih,1)
            dro(ih,3) = zero
            dro(ih,4) = zero 
	        drvp(ih)  = zero  
         enddo
      elseif (itx.eq.2) then
         do ih = 0,ngh
            ro(ih,1)  = + rs(ih,1) + rs(ih,2)
            ro(ih,2)  = + rv(ih,1) + rv(ih,2)
            ro(ih,3)  = - rs(ih,1) + rs(ih,2)
            ro(ih,4)  = - rv(ih,1) + rv(ih,2)
   	        dro(ih,1) = drs(ih,1) + drs(ih,2)
	        dro(ih,2) = drv(ih,1) + drv(ih,2)
            dro(ih,3) = - drs(ih,1) + drs(ih,2)
            dro(ih,4) = - drv(ih,1) + drv(ih,2)
	        drvp(ih)  = drv(ih,2)
         enddo
      else
	     stop 'in DENSIT: itx not properly defined'
      endif

  100 format(a,i3,2f15.8) 
c
      if (lpr)
     &write(l6,*) ' ****** END DENSIT *********************************'
      return
C-end-DENSIT
      end
c======================================================================c

      subroutine denssh(it,lpr)

c======================================================================c
c
c     calculates densities in oscillator basis 
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tt*8                                            ! bloqua
c
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(nhbx,nkx,4),equ(nkx,4)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN DENSSH *******************************'
c
c
c---- loop over the j-blocks
      sp = zero
      il = 0
      do ib = 1,nb
         kap = kb(ib)
	 nf  = id(ib,1)
	 ng  = id(ib,2)
	 nh  = nf + ng  
	 k1  = ka(ib,it) + 1
	 ke  = ka(ib,it) + kd(ib,it)
         k1a = ka(ib,it+2) + 1
         kea = ka(ib,it+2) + kd(ib,it+2)
	 mul = 2*iabs(kap)
	 m   = ib + (it-1)*nbx
         if (lpr) write(l6,'(/,a,1x,a)') tb(ib),tis(it)
c
c------- contributions of large components f*f to rho
         do n2 =  1,nh
         do n1 = n2,nh
            sr = zero
            do k = k1,ke
               sr = sr + fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
            enddo
            do k = k1a,kea
               sr = sr + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
            enddo
	    sr = sr*mul
            rosh(n1+(n2-1)*nh,m) = sr
            rosh(n2+(n1-1)*nh,m) = sr    
c
         enddo   ! n1
         enddo   ! n2
c
c------- contributions of large components f*f to kappa 
	 i0 = il
         do n2 =  1,nf
         do n1 = n2,nf
            i12 = (2 - n2/n1)
            il  = il + 1
            sk = zero
            do k = k1,ke
               sk = sk + fguv(nh+n1,k,it)*fguv(n2,k,it)
            enddo
            do k = k1a,kea
               sk = sk + fguv(nh+n1,k,it+2)*fguv(n2,k,it+2)
            enddo
            sk = mul*sk
            aka(il,it) =  i12*sk
	    if (n1.eq.n2) sp = sp + aka(il,it)
         enddo   ! n1
         enddo   ! n2

      enddo  ! ib     
      spk(it) = half*sp

      if (lpr)
     &write(l6,*) ' ****** END DENSSH *********************************'
      return
c-end-DENSSH
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
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /initia/ inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /vvvikf/ mv,ipos(nbx),nib(mvx),nni(2,mvx)
c
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
      if (inink.eq.0) then
c
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         read(laka) mv0
         if (mv0.ne.mv) stop 'in DINOUT: mv wrong'
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
c     
c---- initial constant pairing field
      elseif (inink.eq.1) then
         do it = 1,itx
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
      else
         stop 'in DINOUT: inink wrong'
      endif   ! inink
      if (lpr) write(l6,*) ' Pairing field has been calculated'
c
c==== writing of the pairing potential
      elseif (is.eq.2) then
c
c
         open(laka,file='dirhb.del',form='unformatted',status='unknown')
         write(laka) mv
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
c
      else
         stop 'in DINOUT: is wrong'
      endif   ! is
c
      if (lpr) then
         do it = 1,itx
            do ib = 1,nb
	           nf = id(ib,1) 
	           ng = id(ib,2)
	           nh = nf + ng
	           m  = ib + (it-1)*nbx
	           k0 = ia(ib,1)+1
	        enddo   ! ib
	     enddo   ! it
         write(l6,*) ' ****** END DINOUT ****************************'
      endif
      return
c-end-DINOUT
      end

c======================================================================c

      subroutine dirhb(it,lprh,lprl)

c======================================================================c
c
c     solves the Dirac-HFB-Equation in spherical oscillator basis
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
c
      logical lprh,lprl
c
      character*1 bbb
      character*8 tbb(nhbx)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character tk*8                                            ! blolev
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
      character nucnam*2                                        ! nucnuc
c
      dimension hb(nhbqx),e(nhbx),ez(nhbx)
c
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /deldel/ de(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /gamgam/ hh(nhhx,nb2x)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /nucnuc/ amas,nmas,npr(2),nucnam
      common /physco/ hqc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /waveuv/ fguv(nhbx,nkx,4),equ(nkx,4)
c
c
      data maxl/200/,epsl/1.d-8/,bbb/'-' /
c
      if (lprh.or.lprl)
     &write(l6,*) ' ****** BEGIN DIRHB **************************'

c======================================================================c
c     with pairing: Dirac-Bogoliubov equation
c======================================================================c
c------- loop over the different j-blocks
      dl    = 100.d0
      xh    = ala(it) + dl
      xl    = ala(it) - dl
      al    = ala(it)  
      sn    = zero
      do lit=1,maxl
         snold=sn
         sn=zero
         alx=al
         klp=0
         kla=0
         do ib = 1,nb
	    kap  = kb(ib)
	    j    = iabs(kap)
	    mul  = mb(ib)
	    lf   = lfkap(kap)
	    lg   = lgkap(kap)
            nf   = id(ib,1)
            ng   = id(ib,2)
            nh   = nf + ng
            nhb = nh + nh
            m    = ib + (it-1)*nbx
c
c---------- calculation of the Dirac-HFB-Matrix:
            do i2 = 1,nh
               do i1 = i2,nh
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
c---------- Diagonalization:
            if (lprh) then
	       i0f = ia(ib,1)  
               do i = 1,nh
                  tbb(i)    = tt(i0f+i)
               enddo
               do i = 1,nh
                  tbb(nh+i) = tbb(i)
               enddo
               write(l6,'(/,a)') tb(ib)
            endif
            call sdiag(nhb,nhb,hb,e,hb,ez,+1)
c---------- store eigenvalues and wave functions
c---------- particles
	    ka(ib,it) = klp
            do k = 1,nf
	       klp = klp + 1
	       ibk(klp,it) = ib
               write(tk(klp,it),100) k,tl(lf),j+j-1
  100       format(i2,a1,i2,2h/2)
               equ(klp,it) = e(nh+k)
               do i = 1,nhb
                  fguv(i,klp,it) = hb(i+(nh+k-1)*nhb)
               enddo
               v2 = zero
               do i = 1,nh
                  v2 = v2 + fguv(nh+i,klp,it)**2
               enddo
               sn=sn+v2*mul
            enddo
	    kd(ib,it) = klp - ka(ib,it)
c
c---------- antiparticles
	    ka(ib,it+2) = kla
            do k = 1,ng
	       kla = kla + 1
	       ibk(kla,it+2) = ib
            write(tk(kla,it+2),100) kla-ka(ib,it+1),tl(lg),j+j-1
               equ(kla,it+2) = e(ng-k+1) 
               do i = 1,nhb
                  fguv(i,kla,it+2) = hb(i+(ng-k)*nhb)
               enddo
               v2 = zero
               do i = 1,nh
                  v2 = v2 + fguv(nh+i,kla,it+2)**2
               enddo
               sn = sn + v2*mul 
            enddo
	    kd(ib,it+2) = kla - ka(ib,it+2)
c
      enddo   ! ib
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
      if (abs(al-alold).lt.epsl) goto 30
c
      if (lprl.or.lit.gt.10) then
          write(l6,113) lit,'. L-Iteration: ',bbb,alold,dn,al
          write(6,113)  lit,'. L-Iteration: ',bbb,alold,dn,al
  113     format(i4,a,a1,3f13.8)
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
  101    format(i4,a,i4,3f13.8)
      endif
      ala(it) = al
      nk(it)   = klp
      nk(it+2) = kla
      
      if (lprh.or.lprl)
     &write(l6,*) ' ****** END DIRHB *********************************'
      return
C-end-DIRHB
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
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
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
c     include 'diz.par'
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
c---- mathemathical constants
c     data zero/0.0d0/,one/1.d0/,two/2.d0/
c     data half/0.5d0/,third/0.333333333333333333d0/
c     data pi/3.141592653589793d0/
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
      if(ibroyd.eq.0) then
         if (inxt.eq.ii) then
            write(6,*) 
     &      'next stop? (0 right now, >0 fixed xmix, <0 autom. xmix)'
            read(*,*)  inx

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
             xmix = xmix * 1.04
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
c     include 'diz.par'
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
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
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

      double precision function racslj(K,L1,L2,J2,J1)

C=======================================================================
C
C     Calculates the Racah-coefficient     (  k   l1   l2 )
c                                          ( 1/2  j2   j1 )
c
C     for integer values        k = K, l1 = L1,     l2 = L2
C     and half integer values          j1 = J1-1/2, j2 = J2-1/2
C
C     Method of Edmonds
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
c     include 'diz.par'
      parameter (igfv = 100)
c
      common /gfviv / iv(0:igfv)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
C
      w(i,j,k,l,m,n) = sq(i)*sq(j)*sqi(k)*sqi(l)*sqi(m)*sqi(n) 
c
      racslj = 0.d0
c
      l12m = l1 - l2
      l12p = l1 + l2
      j12m = j1 - j2
      j12p = j1 + j2 - 1
c     check of triangular rule
      if ( (iabs(l12m).gt.k .or. k.gt.l12p) .or.
     &     (iabs(j12m).gt.k .or. k.gt.j12p) ) return
c
      if (j1.eq.l1+1) then
         if (j2.eq.l2+1) then
            racslj = -w(j12p+k+1,j12p-k,2*j1-1,j1,2*j2-1,j2)
         elseif (j2.eq.l2) then
            racslj =  w(k-l12m,k+j12m,2*j1-1,j1,l2,2*l2+1)
         endif
      elseif (j1.eq.l1) then
         if (j2.eq.l2+1) then
            racslj =  w(k+l12m,k-j12m,2*j2-1,j2,l1,2*l1+1)
         elseif (j2.eq.l2) then
            racslj =  w(l12p+k+1,l12p-k,l1,2*l1+1,l2,2*l2+1)
         endif
      endif      
      racslj = iv(l12p+k)*racslj/2
c
      return
c-end-RACSLJ
      end
c======================================================================c

      double precision function rtbrent(x1,x2,y1,y2,tol)      

c======================================================================c
c
c     Using Brent's method find the root of the function FUNC(X) 
c     known to lie between X1 and X2.
c     Y1 = FUNC(X1), Y2 = FUNC(X2) 
c     The root will be returned as RTBRENT with accuracy TOL
c
c     from: NUMERICAL RECIPIES, 9.3
c
c----------------------------------------------------------------------c 
      implicit real*8(a-h,o-z)
      parameter (itmax = 100, eps = 1.d-12)
c     maximum number of iterations, machine floating point precision
c      external func
      data zero/0.0/,one/1.d0/,half/0.5d0/
c
      a  = x1
      b  = x2
      fa = y1
      fb = y2
      if (fa*fb.gt.zero) stop ' in RTBRENT: root must be bracketed'
      fc = fb
      do 10 iter = 1,itmax
         if (fb*fc.gt.zero) then
c           rename a,b,c and adjust bounding interval
            c  = a              
            fc = fa
            d  = b - a
            e  = d
         endif
         if (abs(fc).lt.abs(fb)) then
            a  = b
            b  = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
         endif
c       
c        convergence check
         tol1 = 2*eps*abs(b)+half*tol
         xm = half*(c-b)
         if (abs(xm).le.tol1 .or. fb.eq.zero) then
            rtbrent = b
            return
         endif
c
         if (abs(e).ge.tol1. and. abs(fa).gt.abs(fb)) then
c           attempt inverse quadratic interpolation
            s = fb/fa
            if (a.eq.c) then
               p = 2*xm*s
               q = one - s
            else
               q = fa/fc
               r = fb/fc 
               p = s*(2*xm*q*(q-r) - (b-a)*(r-one))
               q = (q-one)*(r-one)*(s-one)
            endif
            if (p.gt.zero) q = -q
c           check whether in bounds
            p = abs(p)
            if (2*p.lt.dmin1(3*xm*q-abs(tol1*q),abs(e*q))) then
c              accept interpolation
               e = d
               d = p/q
            else
c              interpolation failed, use besection
               d = xm
               e = d
            endif
         else
c           bounds decreasing too slowly, use bisection
            d = xm
            e = d
         endif
c        move last best guess to a
         a  = b
         fa = fb
         if (abs(d).gt.tol1) then
c           evaluate new trial root
            b = b + d
         else
            b = b + sign(tol1,xm)
         endif
c         fb = pnum(b)
c        call numas(b,rs,av,pr,.false.)        
c        fb = pr
   10 continue
      stop ' in RTBRENT: exceeding maximum number of iterations'
      rtbrent = b
c
c-end-RTBRENT
      end
c======================================================================c

      subroutine brak0(func,x0,x1,x2,f1,f2,step)

c======================================================================c
c
c     subroutine for bracketing a root of a function
c
c     given a function monotonous groving f(x) = FUNC and 
c     given an initial point X0 this routine searches for an interval
c     such that ther is root of f(x) between x1 and x2
c     x1 < x2
c
c
c     INPUT:  X0   starting point
c             FUNC monotonously growing function of x
c             STEP inital step-size
c
c     OUTPUT: X1,X2 an intervall, which brakets a zero of f(x)
c             f1,f2 values of the function f(x) at x1 and x2
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      external func
c
      data maxit/100/
      x  = x0
      dx = step 
      s  = func(x)
      if (s.lt.0.0) then
	 x1 = x
	 f1 = s
	 do i = 1,maxit
	    x = x + dx 
            s = func(x)
	    if (s.gt.0.0) then
	       x2 = x
	       f2 = s
	       goto 20
            else
	       x1 = x
	       f1 = s
            endif
	    dx = 2*dx
         enddo   ! i
	 stop 'in BRAK0 no success in bracketing'
      else
	 x2 = x
	 f2 = s
	 do i = 1,maxit
	    x = x - dx 
            s = func(x)
	    if (s.lt.0.0) then
	       x1 = x
	       f1 = s
	       goto 20
            else 
	       x2 = x
	       f2 = s
            endif
	    dx = 2*dx
         enddo   ! i
	 stop 'in BRAK0 no success in bracketing'
      endif
   20 continue      
c     write(6,*) 'success in bracketing:     ',x1,x,x2,f1,s,f2
c 100 format(a,3f16.8,/,27x,3f16.8)
c
      return
c-end-BRAK0
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
C=======================================================================

      double precision function wiglll(l1,l2,l3)

C=======================================================================
C
C     Calculates the Wigner-coefficient    ( l1  l2  l3 )
c                                          (  0   0   0 )
c     for integer values of l1,l2,l3
c
C     Method of Edmonds
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
c     include 'diz.par'
      parameter (igfv = 100)
c
      common /gfvfak/ fak(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfviv / iv(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
c
      wiglll = 0.d0
c
      l = l1 + l2 + l3
      if (mod(l,2).ne.0) return
      l12p = l1 + l2
      l12m = l1 - l2
      lh   = l/2
c
c     check of triangular rule
      if (iabs(l12m).gt.l3 .or. l3.gt.l12p) return
c
      wiglll = iv(lh)*wf(l12p-l3)*wf(l12m+l3)*wfi(l+1)*wf(l3-l12m)* 
     &         fak(lh)*fi(lh-l1)*fi(lh-l2)*fi(lh-l3)
c
      return
c-end-WIGLLL
      end
c======================================================================c
      integer function lfkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lfkap = kappa
      else
         lfkap = - kappa - 1
      endif
      return
c-end-LFKAP
      end
c======================================================================c
      integer function lgkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lgkap = kappa - 1
      else
         lgkap = - kappa
      endif
      return
c-end-LGKAP
      end
c======================================================================c
      integer function lfgkap(kappa,is)
c======================================================================c
      if (is.eq.1) then
         if (kappa.gt.0) then
            lfgkap = kappa
         else
            lfgkap = - kappa - 1
         endif
      else
         if (kappa.gt.0) then
            lfgkap = kappa - 1
         else
            lfgkap = - kappa
         endif
      endif
      return
c-end-LFGKAP
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

c======================================================================c
c
      subroutine expect(lpr)
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
      dimension ekt(3),epart(3),xn(3),xs(3),r2(3),ept(3) 
      dimension emes(4)
c
      common /baspar/ hom,hb0,b0
      common /erwar / ea,rms
      common /fermi / ala(2),tz(2)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam 
      common /pair  / del(2),spk(2),spk0(2)
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /cenmas/ ecmd(3),ecmn(3),ecm(3)
      
c
      if (lpr) write(l6,*) '  ****** BEGIN EXPECT *******'
c
c
c======================================================================c
c---- particle number and radii
c======================================================================c
      do it = 1,3
         xn(it) = zero
         xs(it) = zero
         r2(it) = zero
      enddo
      do it = 1,itx
	 s0 = zero
	 s2 = zero
         do ih = 1,ngh
            r  = rb(ih)
            rr = r*r
            x  = rv(ih,it)*wdcor(ih)
            xn(it) = xn(it) + x
            r2(it) = r2(it) + x*rr
            x  = rs(ih,it)*wdcor(ih)
            xs(it) = xs(it)  + x
         enddo   ! ih
      enddo   ! it
      if (itx.eq.1) then
	     xn(2) = xn(1)
	     xs(2) = xs(1)
	     r2(2) = r2(1)
      endif
c
      do it = 1,2
         r2(it) = sqrt(r2(it)/xn(it))
      enddo
c
      xn(3) = xn(1) + xn(2)
      xs(3) = xs(1) + xs(2)
      r2(3) = sqrt((npr(1)*r2(1)**2+npr(2)*r2(2)**2)/amas)
      rc    = sqrt(r2(2)**2 + 0.64)
      rn    = r2(1)
      rms   = r2(3)
c
c
c======================================================================c
c---- single particle energies, kinetic energies and pairing energies
c======================================================================c
      do it = 1,itx
c------- kinetic energy 
         ekt(it) = ekin(it)
c------- particle energy
         epart(it) = epar(it)
c------- pairing energy
         ept(it) = epair(it)
      enddo   ! it
      if (itx.eq.1) then
	     ekt(2) = ekt(1)
	     ept(2) = ept(1)
	     epart(2) = epart(1)
      endif
      ekt(3)   = ekt(1) + ekt(2)
      ept(3)   = ept(1) + ept(2)
      epart(3) = epart(1) + epart(2)
c======================================================================c
c---- field energies
c======================================================================c
      call efield(emes,er,ecou)
      escsc  = emes(1)
      escve  = emes(2)
      evesc  = emes(3)
      eveve  = emes(4)
c
c
c======================================================================c
c---- Total energy
c======================================================================c
      etot0 = ekt(3) + escsc + escve + evesc +  eveve + ecou + ept(3) 
      etot  = etot0 + ecm(3)
      etot1 = epart(3)-escsc-escve-eveve-ecou-er+ept(3)
      etest = etot1-etot0
      ea    = etot/amas
c
c
c======================================================================c
c---- printout
c======================================================================c
      if (.not.lpr) return
c
      write(l6,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
c
c     particle number
      write(l6,'(a,6x,3f15.6)') ' particle number',xn
c
c     trace scalar density
      write(l6,'(a,1x,3f15.6)') ' trace scalar density',xs      
c
c     Lambda
      write(l6,100) ' lambda        ',ala

c
c     trace of kappa
      write(l6,100) ' spk           ',spk
c     
c     rms-Radius    
      write(l6,100) ' rms-Radius    ',r2
c
c     charge-Radius    
      write(l6,'(a,22x,f15.6)') ' charge-Radius ',rc
c
c
c     single-particle energy
      write(l6,100) ' Particle Energ',epart
      write(l6,'(a,30x,f15.6)') ' Selfconsistency Test.',etest
c
c     kinetic energy
      write(l6,*)
      write(l6,100) ' Kinetic Energy',ekt
c 
c     sigma energy 
      write(l6,101) ' E-scsc        ',escsc
c
c     omega energy  
      write(l6,101) ' E-scve        ',escve  
c
c     delta-energy
      write(l6,101) ' E-vesc        ',evesc
c 
c     rho-energy       
      write(l6,101) ' E-veve        ',eveve 
c 
c     rearrangement energy
      write(l6,101) ' E-rearrang.   ',er       
c
c     Coulomb energy (direct part)
      write(l6,101) ' Coulomb direct',ecou 
c
c     pairing energy
      write(l6,100) ' Pairing Energy',ept
c
c     total energy without center of mass correction
      write(l6,102) ' Sum without E-cm',etot0
c
c     center of mass correction
      if (icm.eq.0) then
         write(l6,101) ' E-cm  3/4*hom ',ecm(3)
      elseif (icm.eq.1) then
         write(l6,101) ' E-cm  3/4*hom ',ecm(3)
      elseif (icm.eq.2) then
         write(l6,100) 'DE-cm <P**2>/2M',ecmd
         write(l6,100) 'NE-cm <P**2>/2M',ecmn
         write(l6,100) ' E-cm <P**2>/2M',ecm
      else
         stop 'in ERWAR: icm not properly defined'
      endif
c
c     total energy
      write(l6,101) ' Total Energy  ',etot
c
      write(l6,101) ' E/A           ',ea 
      write(l6,*)
c
  100 format(a,7x,3f15.6)
  101 format(a,37x,f15.6)
  102 format(a,35x,f15.6)
  200 format(a,3f15.6)
  999 format(a,i4,a)
c
      write(l6,*) ' ****** END EXPECT *********************************'
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
      dimension h0(nhhx)
      common /baspar/ hom,hb0,b0
      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /physco/ hbc,alphi,r0
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /single/ sp(nfgx,nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
c
c
      ek = zero
      do ib = 1,nb
         kap = kb(ib)
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
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
         ek = ek + trabt(nh,nh,nh,nh,h0,rosh(1,m))
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
      character tb*5                                            ! blokap
      character tk*8                                            ! blolev
c
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /physco/ hbc,alphi,r0
      common /mathco/ zero,one,two,half,third,pi 
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
      common /gamgam/ hh(nhhx,nb2x)

c
      ep = zero
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         m   = ib + (it-1)*nbx
         ep = ep + trabt(nh,nh,nh,nh,hh(1,m),rosh(1,m))
      enddo   ! ib
      epar = ep
c
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
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /rhoshe/ rosh(nhhx,nb2x),aka(mvx,2)
      common /deldel/ de(nhhx,nb2x)

      s = zero
      il=0
      do ib = 1,nb
         kap = kb(ib)
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         m   = ib + (it-1)*nbx
         do n2=1,nf
            do n1=n2,nf
               il=il+1
               s=s+de(n1+(n2-1)*nh,m)*aka(il,it)!*(-1)**ib
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
      dimension emes(3)
c
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /coupld/ ddmes(4)
      common /couplf/ ff(0:ngh,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /fields/ phi(0:ngh,4)
      common /gaucor/ wdcor(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
c======================================================================c
c---- field energies
c======================================================================c
c
c     meson-fields
c
      if (ipc.eq.0) then
         er=zero
	 do m = 1,4
	    s = zero
            do i = 1,ngh
               s  = s + ggmes(m)*ff(i,m,1)*phi(i,m)*ro(i,m)*wdcor(i)
               er= er + ggmes(m)*ff(i,m,2)*phi(i,m)*ro(i,m)*ro(i,2)
     &                    *wdcor(i)
            enddo   ! i
	    emes(m) = half*hbc*s
         enddo   ! m
         er=er*hbc
c        rearrangement term

c     point coupling
      elseif (ipc.eq.1) then
         er=zero
	 do m = 1,4
	    s = zero
            do i = 1,ngh
               s = s + ggmes(m)*ff(i,m,1)*ro(i,m)**2*wdcor(i)
               er=er+ggmes(m)*ff(i,m,2)*ro(i,m)**2*ro(i,2)
     &                    *wdcor(i)
            enddo   ! i
            do i = 1,ngh
               s = s + ddmes(m)*ro(i,m)*dro(i,m)*wdcor(i)
            enddo   ! i
	    emes(m) = half*hbc*s
         enddo   ! m
         er=er*half*hbc
c        rearrangement term

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
         do i = 1,ngh
	        rv = half*(ro(i,2) + ro(i,4))
            ecou  = ecou + cou(i)*rv*wdcor(i)
         enddo   ! i
      endif   ! icou
      ecou  = hbc*ecou/2
c
c
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
      dimension so(ngh),ph(ngh)
      common /baspar/ hom,hb0,b0
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /couplf/ ff(0:ngh,4,2)
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /fields/ phi(0:ngh,4)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp   
c
      if (ipc.eq.1) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN FIELD ********************************'

      do imes=1,4
         do i = 1,ngh
            so(i) = + ff(i,imes,1)*ro(i,imes)
         enddo
         call gordon(imes,so,phi(1,imes))
      enddo   
 
      if (lpr)
     &write(l6,*) ' ****** END FIELD **********************************'
      return
C-end-FIELD
      end
c======================================================================c

      subroutine gordon(imes,so,phi)

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
      dimension so(ngh),phi(ngh)
C
      common /mathco/ zero,one,two,half,third,pi
      common /propag/ gg(ngh,ngh,4)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c     write(l6,*) ' ****** BEGIN GORDON *******************************'
c
c---- multiplication of source with Greens function
      do i = 1,ngh
         s = zero
         do k = 1,ngh
            s = s + gg(i,k,imes)*so(k)
         enddo
         phi(i) = s
      enddo
c
c     write(l6,*) ' ****** END GORDON *********************************'
      return
c-end-GORDON
      end
c======================================================================c

      subroutine greemes(lpr)

c======================================================================c
C
C     calculation of the meson and Coulomb-propagator
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
c
c
      dimension gi(nox,nox),rr(nox,ngh)
c
      common /baspar/ hom,hb0,b0
      common /bosqua/ no
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,ames(4)
      common /physco/ hbc,alphi,r0
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /propag/ gg(ngh,ngh,4)
      common /bospol/ rnb(1:nox,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GREEMES ******************************'
c
      if (ipc.eq.1) return
c
c     meson-propagators
      do imes = 1,4
         f = (one/((ames(imes)+1.d-10)*b0))**2
         do i = 1,no
         do k = 1,no
            gi(i,k) = zero
         enddo   ! k
         enddo   ! i
         do n = 1,no
            gi(n,n) = one + f*(2*n-half) 
            if (n.lt.no) then
               gi(n,n+1) = f*sq(n)*sqh(n)
               gi(n+1,n) = gi(n,n+1)
            endif
            do ih = 1,ngh
               rr(n,ih) = rnb(n,ih)
            enddo   ! ih
         enddo   ! n
         call lingd(nox,nox,no,ngh,gi,rr,d,ifl)
         f0 = one/(4*pi*b0**3)
         do kh = 1,ngh
	        f = f0*wdcor(kh)
            do ih = 1,ngh
               s = zero
               do n = 1,no
                  s = s + rnb(n,ih)*rr(n,kh)
               enddo
               gg(ih,kh,imes) = f * s
            enddo
         enddo      
      enddo   ! imes
c
      if (lpr)
     &write(l6,*) ' ****** END GREEMES ********************************'
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
c
      include 'dirhb.par'
c       
      logical lpr
      character parname*10                     ! common partyp

      common /dforce/ a_s,a_v,a_ts,a_tv,b_s,b_v,b_ts,b_tv,
     &                c_s,c_v,c_ts,c_tv,d_s,d_v,d_ts,d_tv,dsat
      common /masses/ amu,amsig,amome,amdel,amrho
      common /coupld/ ddsig,ddome,dddel,ddrho
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/  gl(2),gal
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN FORCES *******************************'
c---- forces
c=============================================================== 
       if (parname.eq.'DD-ME2') then
         amu    =  939.d0                 ! MeV
         amsig  =  550.12380d0            ! MeV
	     amome  =  783.d0                 ! MeV
	     amrho  =  763.d0                 ! MeV
         gsig   =   10.5396d0
         gome   =   13.0189d0
         gdel   =   zero
	     grho   =    3.6836d0
         b_s    =    1.0943d0
	     c_s    =    1.7057d0
         c_v    =    1.4620d0
         a_tv   =    0.5647d0
         dsat   =    0.152d0
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

         ipc    =  0
	 icm    =  2
	 idd    =  2
c=============================================================== 
       elseif (parname.eq.'DD-PC1') then
c---------------------------------------------------------------
c        G(x) = a + (b + c*x) * exp(-d*x)
c---------------------------------------------------------------

         dsat   =  0.152d0              ! fm^-3
         amu    =  939.d0               ! MeV
c
c        scalar-isoscalar
         a_s    = -10.0462d0           ! fm^-2
         b_s    =  -9.1504d0           ! fm^-2
         c_s    =  -6.4273d0           ! fm^-2
         d_s    =  +1.3724d0            
c
c        vector-isoscalar
         a_v    =  +5.9195d0           ! fm^-2
         b_v    =  +8.8637d0           ! fm^-2
         c_v    =   0.00000d0           ! fm^-2
         d_v    =  +0.6584d0            
c
c        scalar-isovector
         a_ts   =   zero
         b_ts   =   zero
         c_ts   =   zero
         d_ts   =   zero
c
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
         amsig = zero
         amome = zero
         amdel = zero
         amrho = zero

         icm    = 2
         idd    = 2
         ipc    = 1
      else
          stop 'This type of force is not defined'
      endif
c===============================================================      
      amu      = amu/hbc
      amsig    = amsig/hbc
      amome    = amome/hbc
      amdel    = zmdel/hbc
      amrho    = amrho/hbc

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
      
c----------------------------------------------------------------------
c---- Separable pairing force

      gl(1) = -728.d0
      gl(2) = -728.d0
      gal   = 0.415d0
c
c---- printout of force:
      if (lpr) call pripar
c
   10 if (lpr)
     &write(l6,*) ' ****** END FORCES *********************************'
c
      return
c-end-FORCES
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
      include 'dirhb.par'
c          
c
      character parname*10                     ! common partyp
c
      common /partyp/ parname
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /pair  / del(2),spk(2),spk0(2)
      common /coupld/ ddmes(4)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/  gl(2),gal
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
      logical lpr
      include 'dirhb.par'

      character tb*5                                            ! blokap      
c
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /single/ sp(nfgx,nbx)
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /gamgam/ hh(nhhx,nb2x)
c
      emcc2 = 2*amu*hbc
      f = hbc/b0
c
      do it=1,itx  
         do ib=1,nb
            kapp= kb(ib)
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            i0f = ia(ib,1)
            i0g = ia(ib,2)
            m   = ib + (it-1)*nbx
            lf  = lfkap(kapp)
            lg  = lgkap(kapp)
c 
            do i2 = 1,nf
               do i1 = 1,ng
                  hh(nf+i1+(i2-1)*nh,m) = f*sp(i1+(i2-1)*ng,ib)
               enddo
            enddo
            call pot(nf,nh,lf,vps(0,it),hh(1,m))
            call pot(ng,nh,lg,vms(0,it),hh(nf+1+nf*nh,m))
            do i = nf+1,nh
               hh(i+(i-1)*nh,m) = hh(i+(i-1)*nh,m) - emcc2
            enddo
c
c     symmetrize HH
            do i2 = 1,nh
               do i1 = i2+1,nh
                  hh(i2+(i1-1)*nh,m) = hh(i1+(i2-1)*nh,m)
               enddo
            enddo
         enddo !ib
      enddo  !it

      return
C-end-GAMMA
      end
c======================================================================c

      subroutine pot(n,nh,l,v,aa)

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      character tt*8                                            ! bloqua
c
      dimension aa(nh,nh),v(0:ngh)
c
      common /mathco/ zero,one,two,half,third,pi
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)

      do n2 =  1,n
      do n1 = n2,n
         s = zero
         do ih = 1,ngh
            s = s + v(ih)*rnl(n1,l,ih)*rnl(n2,l,ih)
         enddo  ! ih
         aa(n1,n2) = s
      enddo   ! n1
      enddo   ! n2
c
      return
c-end-POT
      end

c======================================================================c

      subroutine gaupol(lpr)

c======================================================================c
c
c     calculates the radial functions for the spherical oscillator
c
c     the wave function phi(nlj) of the spherical oscillator are: 
c
c     phi(r,Omega) = b^(-3/2) * R_nl(r) * Y_ljm(Omega) 
c     
c     R_nl(r) = N_nl * r^l  * L^(l+1/2)_(n-1)(x*x) * exp(-x*x/2)
c
c     N_nl    = sqrt(2 * (n-1)!/(n+l-1/2)!)     and    x=r/b
c
c     the contribution to the density from the shell j is
c
c     rho_j(r)= 1/(4*pi*b0**3) * (2*j+1) * R_nl(r)^2
c
c     the radial function at meshpoint xh(ih) is stored in RNL(n,l,ih)
c     in the following way: RNL is normalized in such way that the
c     norm integral reads
c
c     \int d^3r |phi(r)|^2 = 1 = \sum_i RNL(n,l,i)**2
c
c     this means, that RNL contains the following factors:
c
c     a)  the radial part of the wavefunction r * R_nl(r)
c     b)  the length units factor  b ** (3/2)
c     c)  the Gaussian weight sqrt( WH(i) ): 
c         \inf_0^inf f(x) dx = \sum_i f(x_i) * WH(i)
c
c     having RNL(n,l,i) we get the radial wavefunction:
c
c     R_nl(r) =  RNL(n,l,i) / ( x_i * sqrt(WH(i)) )  
c
c     and the density contribution from the shell j
c
c     rho_j(r) = (2*j+1) * RNL(n,l,i)**2 / ( 4 * pi x_i**2 * WH(i) * b**3)   
c
c----------------------------------------------------------------------c
c

c     RNL1 contains the radial derivatives in the following form:
c
c     d/dr R_nl(r) = 1/b * RNL1(n,l,i) / (x_i * sqrt(WH(i) )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tt*8                                            ! bloqua
c
      common /bosqua/ no
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /gfvsq / sq(0:igfv)
      common /gfvsqi/ sqi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvshi/ shi(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /bospol/ rnb(1:nox,0:ngh)
      common /sdimos/ nrm,nlm,nrbm
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GAUPOL ************************'
c
c     f = 2/pi**0.25
      f = 2*wgi(0)
c
      do ih = 0,ngh
         r  = xh(ih)
         rr = r*r
         ri = one/r 
         fe = f*exp(-half*rr)
c
c------------------------------------
c        the functions rnl,rnl1 contain already the measure
         u1 = fe*sqrt(wh(ih)*rr)
c------------------------------------
c
c------- basis for fermions
         do l = 0,nlm
            rnl(1,l,ih)  = u1
            rnl(2,l,ih)  = u1*(l+1.5d0-rr)*shi(l+1)
            u1           = u1*r*shi(l+1)
            rnl1(1,l,ih) =    (l-rr)*rnl(1,l,ih)*ri
            rnl1(2,l,ih) = ((2+l-rr)*rnl(2,l,ih) - 
     &                       2*sqh(l+1)*rnl(1,l,ih))*ri
c
            do n = 3,nrm
               rnl(n,l,ih)  = ((2*n+l-2.5d0-rr)*rnl(n-1,l,ih) -
     &           sq(n-2)*sqh(n-2+l)*rnl(n-2,l,ih))*sqi(n-1)*shi(n-1+l)
               rnl1(n,l,ih) = ((2*n+l-2-rr)*rnl(n,l,ih) -
     &           2*sq(n-1)*sqh(n-1+l)*rnl(n-1,l,ih))*ri
            enddo
         enddo
c
c------- basis for bosons
         rnb(1,ih)  = fe
         rnb(2,ih)  = fe*(1.5d0-rr)*shi(1)
         do n = 3,no
            rnb(n,ih) = ((2*n-2.5d0-rr)*rnb(n-1,ih) -
     &           sq(n-2)*sqh(n-2)*rnb(n-2,ih))*sqi(n-1)*shi(n-1)
         enddo
c
      enddo 
c
c
c---- Test of orthogonality
      if (lpr) then
         do 40 l = 0,nlm
            write(l6,'(/,80(1h*))')
            do 41 n = 1,nrm
               write(l6,'(a,2i3)') 
     &         ' Radial function and derivative for n,l =',n,l
               ix = 5
               write(l6,'(5f15.8)') (rnl(n,l,ih),ih=1,ix)
               write(l6,'(5f15.8)') (rnl1(n,l,ih),ih=1,ix)
   41       continue
            do 50 n2 = 1,nrm
            do 50 n1 = n2,nrm 
               s1 = zero
               s2 = zero
               s3 = zero
               sb = zero
               do ih = 0,ngh
                  rr = xh(ih)**2
                  s0 = rnl(n1,l,ih)*rnl(n2,l,ih)
                  s1 = s1 + s0
                  s2 = s2 + rr*s0
                  s3 = s3 + (rnl1(n1,l,ih)*rnl1(n2,l,ih)
     &                       + rnl1(n1,l,ih)*rnl(n2,l,ih)/xh(ih)
     &                       + rnl(n1,l,ih)*rnl1(n2,l,ih)/xh(ih)
     &                       + s0*(1+l*(l+1))/rr)
		  if (l.eq.0) then
		     sb = sb + rnb(n1,ih)*rnb(n2,ih)*rr*wh(ih)
                  endif
               enddo
               write(l6,'(a,2i3,4f12.8)') 
     &                  ' RNL(n,l) test ',n1,n2,s1,s2,s3,sb
   50       continue
   40    continue
         write(l6,'(/,80(1h*))')
      endif
c
      if (lpr)
     &write(l6,*) ' ****** END GAUPOL *************************'
      return
c-end-GAUPOL
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
c
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      dimension xt(2*ngh),pt(2*ngh)

      if (lpr)
     &write(l6,*) ' ****** BEGIN GAUSH *******************************'
      call gauher(xt,pt,2*ngh)
      do i=1,ngh
         xh(i)=xt(ngh+1-i)
         ph(i)=pt(ngh+1-i)
         wh(i) = ph(i)*dexp(xh(i)*xh(i))
      enddo
      xh(0)=1.d-10
      ph(0)=1.d-10
      wh(0)=1.d-10
c
c      write(l6,100) ngh
C
      if (.not.lpr) return
      write(l6,101)
c
      write(l6,*) ' xh'
      write(l6,105) (xh(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (xh(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (xh(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (xh(i),i=1,ngh)
c
      write(l6,*) ' wh'
      write(l6,105) (wh(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (wh(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (wh(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (wh(i),i=1,ngh)
c
      write(l6,*) ' ph'
      write(l6,105) (ph(i),i=1,ngh)
      if (ngh.eq.8)  write(l6,102) (ph(i),i=1,ngh)
      if (ngh.eq.12) write(l6,103) (ph(i),i=1,ngh)
      if (ngh.eq.16) write(l6,104) (ph(i),i=1,ngh)
c
  100 format('  GAUSH:  G-H-Integration  ngh =',i3)
  101 format(1x,36(1h-))
  102 format(2(3e19.11,/),2(e19.11))
  103 format(4(3e19.11,/))
  104 format(5(3e19.11,/),e19.11,/)
  105 format(3e19.11)
c
      if (lpr)
     &write(l6,*) ' ****** END GAUSH *********************************'
      return
c-end-GAUSH
      end

c======================================================================c
c
      subroutine gdd(lpr)
c
c----------------------------------------------------------------------c
c
c     calculates 
c        fmes(x,1)       density-dependent coupling constants
c     and
c        fmes(x,2)       their derivatives
c
c======================================================================c
      implicit real*8(a-h,o-z)
c
      include 'dirhb.par'
c     
      logical lpr

      common /couplf/ ff(0:ngh,4,2)
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /couplg/ ggmes(4),lmes(4)
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat 
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c     
      fun1(x,a,b,c,d) = a*(1+b*(x+d)**2)/(1+c*(x+d)**2)
      dun1(x,a,b,c,d) = 2*a*(b-c)*(x+d)/(1+c*(x+d)**2)**2
 
      fun2(x,a,b,c,d) = exp(-a*(x-one))
      dun2(x,a,b,c,d) = -a*exp(-a*(x-one))
      
      fun3(x,a,b,c,d) = a+(b+c*x)*exp(-d*x)
      dun3(x,a,b,c,d) = c*exp(-d*x)-d*(b+c*x)*exp(-d*x)

      if (lpr) then
      write(l6,*) '****** BEGIN GDD ***********************************'
      endif
      if (ipc.eq.0) then
         do i = 0,ngh
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
      elseif (ipc.eq.1) then
         do i = 0,ngh
            x=ro(i,2)/dsat
            do m = 1,4
               ff(i,m,1) = fun3(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun3(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
         enddo   ! i
      endif   ! ipc
c
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
c
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (is.eq.1.and.inin.ne.0) return
c
      if(lpr)
     &write(l6,*) ' ****** BEGIN INOUT ********************************'
c
c---- reading of meson fields from tape:
c-------------------------------------
      if (is.eq.1) then
         open(lwin,file='dirhb.wel',status='old',form='unformatted')
         read(lwin) ng0
         if (ng0.ne.ngh) stop 'in INOUT: ngh wrong'
         read(lwin) ala
 
c------- reading of the potentials 
         read(lwin) vms
         read(lwin) vps
         close(lwin)

         write(l6,*) ' potentials read from tape ','dirhb.wel'
 
c---- writing to tape:
      else
         open(lwou,file='dirhb.wel',status='unknown',form='unformatted')
         write(lwou) ngh
         write(lwou) ala
         write(lwou) vms
         write(lwou) vps
         close(lwou)

         if(lpr) write(l6,*) ' potentials written to tape dirhb.wel'
      endif
c
      if(lpr)
     &write(l6,*) ' ****** END INOUT **********************************'
      return
c-end-INOUT
      end      
c======================================================================c

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
      common /erwar / ea,rms
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /initia/ inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /broyde2/ ibroyd

c
      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
      write(l6,*) ' ****** BEGIN ITER *********************************'

      ii=0
      call gamma()
      call broyden(.false.)
      
      do ite = 1,maxi
         ii = ite
c
         if (lpr) then
            write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,'  mix =',xmix  
            if (l6.ne.6) 
     &      write(6,102)  ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,'  mix =',xmix
  102        format(i3,a,f12.7,2(a,f9.4),a,f5.2) 
         endif

c------- loop over neutrons and protons
         do it = 1,itx
c            
c------- diagonalization of the Dirac-Bogolibov equation	
            call dirhb(it,.false.,.false.)
c     
c---------- calculation of new densities in oscillator basis
            call denssh(it,.false.)
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
         call poten(.false.)

c
c------- potentials in r-space
         call gamma()
c
c------- pairing field
         do it = 1,itx
	    call delta(it,.false.)
            spk0(it) = spk(it)
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
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text2,si
      goto 40
c
   30 write(6,101) nucnam,nmas,text2,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ii,text3,si
c
  100 format(1x,68(1h*),/,1x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(a2,i4,a27,i4,a14,f17.10)
   40 write(l6,*) ' ****** END ITER ***********************************'
      return
c-end-ITER
      end

c======================================================================c
c 
c     PROGRAM DIRHB-Spherical  
c
c======================================================================c
c     Relativistic Hartree-Bogoliubov theory in a spherical basis
c     Main part of the code
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
c
c---- sets data
      call default
c
c---- reads in data     
      call reader
c
c---- force-parameters
      call forces(.true.)
c
c---- Gauss-Hermite mesh points
      call gaush(.false.)
c
c---- preparations
      call prep(.false.)
c
c---- initialization of the potentials
      call inout(1,.false.)
      call start(.false.)
c
c---- oscillator basis for single particle states
      call base(.false.)
c
c---- initialization of the pairing field
      call dinout(1,.false.)
c
c---- wavefunctions at Gauss-Meshpoints
      call gaupol(.false.)
c
c---- single-particle matix elements
      call singf(.false.)
c
c---- pairing matix elements
      call singd(.false.)
c      
c---- coulomb and meson propagators
      call greecou(.false.)
      call greemes(.false.)
c
c---- iteration
      call iter(.true.)
c
c---- transformation to the canonical basis
      call canon(.true.)
c
c---- center of the mass correction
      call centmas(.false.)
c
c---- results
      call resu(.true.)
c
c---- densities in the coordinate space
      call plot(.false.)
c---- punching of the potentials to the tape  dirhb.wel
      call inout(2,.false.)
c---- punching of the pairing field to the tape  dirhb.wel      
      call dinout(2,.false.)
      stop ' FINAL STOP OF DIRHBS'
c-end-DIRHBS
      end

c=====================================================================c

      subroutine plot(lpr)

c=====================================================================c
C
C     prepares plot of densities in coordinate space
C
c---------------------------------------------------------------------c
      include 'dirhb.par'
c
      implicit real*8 (a-h,o-z)
      logical lpr
c
      dimension pn(nox), ppn(nox),zp(0:ngh)
      
c
      parameter(yp0=1.d30)
      parameter(ypn=1.d30)
      common /baspar/ hom,hb0,b0
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /rhorho/ rs(0:ngh,2),rv(0:ngh,2)
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /potpot/ vps(0:ngh,1:2),vms(0:ngh,1:2)
      common /physco/ hbc,alphi,r0
      common /masses/ amu,amsig,amome,amdel,amrho
      common /gaucor/ wdcor(0:ngh)          
c
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN PLOT ********************************'
c
c     number of points for the plot
      mxpl  = 128 
c     plot step in (fm)
      stpl = 0.1
      am = amu*hbc


c---- plot for densities:

      open(lplo,file='dirhb.plo',status='unknown')

      call spline(rb,ro(0,2),zp,ngh,yp0,ypn)  
      r = 1.d-8
      do ist = 0,mxpl
         call splint(0,rb,ro(0,2),zp,ngh,r,y,y1,y2)
         write(lplo,100) r,y
         r = r + stpl
      enddo
      
100          format(f10.3,f15.6) 
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
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /coupld/ ddmes(4)
      common /couplf/ ff(0:ngh,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /dens  / ro(0:ngh,4),dro(0:ngh,4)
      common /fields/ phi(0:ngh,4)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN POTEN ********************************'
c
      do i = 0,ngh
c
c---- meson-fields
         if (ipc.eq.0) then
	    do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*phi(i,m)
	    enddo   ! m
c
c---- rearangement field
	    re = zero
	    do m = 1,4
               re = re + ggmes(m)*ff(i,m,2)*phi(i,m)*ro(i,m)
	    enddo   ! m
            glt(2) = glt(2) + re
c
c---- point-coupling models
	 elseif (ipc.eq.1) then
	    do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*ro(i,m)
	    enddo   ! m
	    do m = 1,4
               glt(m) = glt(m) + ddmes(m)*dro(i,m)
	    enddo   ! m
	    re = zero
	    do m = 1,4
               re = re + ggmes(m)*ff(i,m,2)*ro(i,m)**2
	    enddo   ! m
            glt(2) = glt(2) + half*re
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
c
      enddo   ! i
c
c      if (lpr) then
c         call prigh(0,sig,one,'X(FM) ')
c         call prigh(1,vps(1,1),one,'V+S  n')
c         call prigh(1,vms(1,1),one,'V-S  n')
c         call prigh(1,vps(1,2),one,'V+S  p')
c         call prigh(1,vms(1,2),one,'V-S  p')
c	     write(l6,*) 'si =',si
c      endif
c
      if (lpr)
     &write(l6,*) ' ****** END POTEN **********************************'
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
      character tp*1,tis*1,tit*8,tl*1                    ! common textex
      character nucnam*2                                 ! common nucnuc
      character tb*5
      logical lpr

c
      common /baspar/ hom,hb0,b0
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /gaucor/ wdcor(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      write(l6,*) ' ****** BEGIN PREP *********************************'
c
c---- basis parameters
      hb0 = hbc/(two*amu)
      hom = 41.0*amas**(-third)    
      if (icm.eq.1) hb0 = hb0*(one - one/amas)
      b0 = sqrt(two*hb0/hom)
c
      write(l6,*) ' '
      write(l6,100) ' hom =                         ',hom
      write(l6,100) ' hb0 =                         ',hb0
      write(l6,100) ' b0  =                         ',b0 
c
      do ih = 0,ngh
         rb(ih) = xh(ih)*b0
c        metric element for three-dimensional integration 
         wdcor(ih) = b0**3 * 4*pi * xh(ih)**2 * wh(ih)
      enddo
c
c---- printout pairing:
c      write(l6,*)
  100 format(a,4f11.6)
c      write(l6,100) ' Initial Gap   = ',del
c
   10 if (itx.eq.1) icou = 0
      if (icou.eq.0) write(l6,100) ' without Coulomb force'
      if (icou.eq.1) write(l6,100) ' with Coulomb force'
      if (icou.eq.2) write(l6,100) ' with Coulomb force with exchange'
      if(icm.eq.2) write(l6,*) 'With microscopic c.m. correction'
  101 format(a,i4)


      write(l6,*) ' ****** END PREP ***********************************'
      return
c-end PREP
      end

c======================================================================c

      subroutine reader

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c 
      character parname*10                                      ! partyp
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc

      CHARACTER*8 date
      CHARACTER*10 time
      CHARACTER*5 zone
      INTEGER*4 VALUES(8)
c
      common /basnnn/ n0f,n0b
      common /baspar/ hom,hb0,b0
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /wplot/ itpl,jpl,ippl,kpl      
c............................
caf ... reading from dirhb.dat
c............................
c
      if (lin.eq.0) return
      open(lin,file='dirhb.dat',status='old')
c
c
c---- Basisparameters:            
      read(lin,'(10x,2i5)') n0f,n0b
c
c---- Initialization of wavefunctions:
      read(lin,'(10x,2i5)') inin,inink
c
c---- Nucleus under consideration
      read(lin,'(a2,i4)') nucnam,nama

c---- Initial Gap Parameters
      read(lin,*)
      read(lin,102) del
c
c---- Parameterset of the Lagrangian
      read(lin,*)
      read(lin,'(12x,a10)') parname	
c      read(lin,'(/12x,i4)') itpl
c      read(lin,'(12x,i4)') jpl
c      read(lin,'(12x,i4)') ippl
c      read(lin,'(12x,i4)') kpl
c
c..................................
caf ... end of reading from DIS.dat
c..................................
c

      close(lin)
      call nucleus(2,npro,nucnam)
c
      nneu = nama - npro
      amas = nama 
      
c..............................................................
caf ... begin of output into .OUT file
c..............................................................
c
      if (l6.ne.6) open(l6,file='dirhb.out',status='unknown')
c
      call date_and_time( date, time, zone, values )
      write(l6,'(a)')
      write(l6,'(a)') '  ******************************************  '
      write(l6,'(a)') '  *           Fortran 77 Code              *  '
      write(l6,'(a)') '  *         Spherical H.O. basis            *  '
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
      write(6,'(a,16x,a10)') '  *',parname
      write(l6,'(a)') '  ******************************************  '

      write(6,'(a)')
      write(6,'(a)') '  ******************************************  '
      write(6,'(a)') '  *           Fortran 77 Code              *  '
      write(6,'(a)') '  *         Spherical H.O. basis            *  '
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

      write(l6,*) ' ****** BEGIN READER *******************************'
c
c---- Basisparameters:            
      write(l6,101) ' Number of oscillator shells : ',n0f,n0b
c
c---- Initialization of wavefunctions:
      write(l6,'(a,2i5)') ' Initialization inin,inink   : ',inin,inink  
c
c---- Nucleus under consideration
      write(l6,'(a,a,i4,i6,i4)') ' Nucleus: ',nucnam,nama,nneu,npro
c
c---- Initial Gap Parameters
      write(l6,103) ' Initial Gap Parameters      : ',del
c---- Parameterset of the Lagrangian
      write(l6,106) ' Parameter set               : ',parname
c.............................................................
caf ... end of output into .OUT file
c.............................................................
      tz(1) = nneu
      tz(2) = npro
c
c
  100 format(10x,2i5)
  101 format(a,2i5)
  102 format(10x,2f10.4) 
  103 format(a,2f10.4) 
  106 format(a,'   ',a10) 
  110 format(10x,4f10.4)
  111 format(9x,i1,4f10.4)

      write(l6,*) ' ****** END READER *********************************'
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

      common /blodir/ ka(nbx,4),kd(nbx,4)
      common /waveuv/ fguv(nhbx,nkx,4),equ(nkx,4)
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /fermi / ala(2),tz(2)
      common /erwfit/ etot,rc,rn,slp(2),nn(2),ll(2)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nmas,nneu,npro,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      data eqpmax/100.d0/

c
      if (.not.lpr) return    
c
      write(l6,*) ' ****** BEGIN RESU *********************************'

      call expect(.true.)

c---- quasi-particle energies
      do it = 1,itx
         write(l6,100) tit(it)
         write(l6,101) 'p/h','l j','[nr,l,j]','smax','Eqp','uu','vv'
  100    format(//,' quasi-particle properties ',a,/,1x,66('-'))
  101 format(5x,a,4x,a,5x,a,2x,a,7x,a,7x,a,7x,a)
  
         do ib = 1,nb
            kap = kb(ib)
            j = iabs(kap)
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
            k1 = ka(ib,it) + 1
            k2 = ka(ib,it) + kd(ib,it)
	        i0 = ia(ib,1)
c
            do k = k1,k2
c
c---------- search for main oscillator component
               smax = zero
               su = zero
               sv = zero
               do n=1,nh
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
               if (equ(k,it).lt.eqpmax) 
     &         write(l6,102) k,tph,tk(k,it),tt(i0+imax),smax, 
     &                    equ(k,it),su,sv
  102       format(i4,2x,a1,'  ',a8,'  ','(',a8,')',f5.2,3f10.3) 
   20       enddo  ! k 
         enddo   ! ib
      enddo   ! it
c
      write(l6,*) ' ****** END RESU ***********************************'
c
      return
c-end-RESU
      end

c======================================================================c

      subroutine singd(lpr)

c======================================================================c
c
c     calculates single particle matrix V_nn' for separable-pairing
c     in the spherical oscillator basis
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      dimension oscn(0:n0fx) 
c
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      common /baspar/ hom,hb0,b0
      common /physco/ hbc,alphi,r0
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),
     &                tt(ntx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /vvvsep/ vn(mvx,0:n0fx)
      common /tmrpar/  gl(2),gal
      common /gfviv / iv(0:igfv)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN SINGD ********************************'
c	
      call setwig
      
      call vn0(oscn)
      do i = 1,mvx
	     do nn = 0,n0fx
		    vn(i,nn) = zero
	     enddo    ! nn
	  enddo    ! i
      il = 0
	  do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         m  = ib + (it-1)*nbx
         kappa = kb(ib)
         j     = iabs(kappa)
         l     = lfkap(kappa)
         do n2 = 1,nf
            do n1 = n2,nf
               il = il + 1
               do nn = 0,n1-1+n2-1+l
                  nx = n1-1+n2-1+l-nn
                  vn(il,nn) = iv(l)*talmos(n1-1,l,n2-1,l,nn,
     &                        0,nx,0,0)*oscn(nx)/sqrt(two*l+one)
               enddo                  ! nn
               
            enddo                     ! n1
         enddo                     ! n2 
      enddo                     ! ib
      if (lpr)
     &write(l6,*) ' ****** END SINGD **********************************'
      return
c-end-SINGD
      end  
c======================================================================c

      subroutine setwig

c======================================================================c
c
c     computes and stores in an efficient way couplings for Moshinsky
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
c
      include 'dirhb.par'
c
      integer*2 locs
c
      common /gfvsqh/ sqh(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /gfvwgi/ wgi(0:igfv)
      common /wigwig/ wig(jdim),locs(0:lx,0:lx,0:lx)
c
      r2 = sqh(0)
      do 10 l1 = 0,lx
      do 10 l2 = 0,lx
      do 10 l3 = 0,lx
   10    locs(l1,l2,l3) = jdim

      ico = 0
      do 20 l1 = 0,lx
      do 20 l2 = l1,lx
         kmin = l2+mod(l1,2)
         ktop = min0(lx,l1+l2)
         if (kmin.gt.ktop) goto 20
         do l3 = kmin,ktop,2
            ico = ico+1
            if (ico.gt.jdim)  stop ' in SETWIG: jdim too small'
            ip = (l1+l2+l3)/2
            wig(ico) = wg(ip-l1)*wfi(ip-l1)*wg(ip-l2)*wfi(ip-l2)*
     &                 wg(ip-l3)*wfi(ip-l3)*wf(ip)*wgi(ip+1)*r2
            locs(l1,l2,l3) = ico
            locs(l1,l3,l2) = ico
            locs(l2,l3,l1) = ico
            locs(l2,l1,l3) = ico
            locs(l3,l1,l2) = ico
            locs(l3,l2,l1) = ico
         enddo
   20 continue
      wig(jdim) = 0.0d0
c
      return
c-end-SETWIG
      end
c======================================================================c

      subroutine vn0(oscn)

c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
c
      logical lpr
c
      dimension gll(4),oscn(0:n0fx)
      character tk*8                                            ! blolev
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
      common /baspar/ hom,hb0,b0
      common /blolev/ nk(4),ibk(nkx,4),tk(nkx,4)
      common /bloqub/ ijb(nbx),ilb(nbx,2),ipb(nbx),ikb(nbx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /wavefg/ fg(nhx,nkx,4)
      common /gfvfak/ fak(0:igfv) 
      common /gfviv / iv(0:igfv) 
      common /gfvwfi/ wfi(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /tmrpar/  gl(2),gal
c
	  a=sqrt(gal)
	  f1=one/sqrt(2*pi*a)
	  f1=f1**3
	  f2=a**3/b0**1.5
	  fac=f1*f2
c
      do i = 0,n0fx
         f1 = sqrt(2.d0)/pi/2.d0**0.75d0
         f2 = (one/b0)**1.5d0
         f3 = wg(i+1)
         f4 = wfi(i)
         f5 = (b0**2/(b0**2+a**2))**1.5d0
         f6 = ((b0**2-a**2)/(b0**2+a**2))**i
         s2 = f1*f2*f3*f4*f5*f6
	     oscn(i) = s2 
      enddo    ! i
      
      return      
c
      end
c======================================================================c
c
      function talmos(m1,k1,m2,k2,m3,k3,m4,k4,ilam)
c
c======================================================================c
C
c
C     Talmi - Moshinsky Bracket:( lambda n1 l1 n2 l2 n3 l3 n4 l4 )
C
c     < n1 l1, n2 l2, lambda | n3 l3, n4 l4, lambda >
C
c     see T.A. Brody & M. Moshinsky, "tablas de parentesis de
c     transformacion" (monografias del instituto de fisica, Mexico
c     1960)
c
c     zero is returned for cases that violate the triangle rules
c     energy consevation.
c
c     radial quantum number start from zero:   n1 = 0,1,2,....
c
c     the coefficients are calculated here in the phase convention
c     of Baranger. 
c     <n1 l1 n2 l2 n3 l3 n4 l4, lambda>(Brody-Moszkowsky)
c     = (-)^(b3+n3-lambda) * <n1 l1 n2 l2 n4 l4 n3 l3, lambda>(Baranger)
c       
c----------------------------------------------------------------------c
      implicit real *8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical testj
      integer*2 locs
c
      common /gfviv / iv(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvgmi/ gmi(0:igfv)
      common /gfvsqh/ sqh(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwg / wg(0:igfv)
      common /mathco/ zero,one,two,half,third,pi
      common /wigwig/ wig(jdim),locs(0:lx,0:lx,0:lx)
c
c     itestj is true if any triangle condition is violated
      testj(l1,l2,l3) = l1+l2.lt.l3 .or. iabs(l1-l2).gt.l3
c
      talmos = zero
c
c     check triangular rules and energy
      if (testj(k1,k2,ilam).or.testj(k3,k4,ilam)) return
      nn1 = 2*m1+k1
      nn2 = 2*m2+k2
      nn3 = 2*m3+k3
      nn4 = 2*m4+k4
      ichi = nn1+nn2
      if (ichi.ne.nn3+nn4) return
c
c---- reordering of indices
      if (min0(k1,k2).lt.min0(k3,k4)) then
         if (k1.gt.k2) then
            n3 = m2
            l3 = k2
            n4 = m1
            l4 = k1
            if (nn3.gt.nn4) then
               n1 = m4
               l1 = k4
               n2 = m3
               l2 = k3
               iph = iv(l1+l4)
            else
               n1 = m3
               l1 = k3
               n2 = m4
               l2 = k4
               iph = iv(l1+ilam)
            endif
         else
            n3 = m1
            l3 = k1
            n4 = m2
            l4 = k2
            if (nn3.gt.nn4) then
               n1 = m4
               l1 = k4
               n2 = m3
               l2 = k3
               iph = iv(l3+ilam)
            else
               n1 = m3
               l1 = k3
               n2 = m4
               l2 = k4
               iph = +1
            endif
         endif
      else
         if (k3.gt.k4) then
            n3 = m4
            l3 = k4
            n4 = m3
            l4 = k3
            if (nn1.gt.nn2) then
               n1 = m2
               l1 = k2
               n2 = m1
               l2 = k1
               iph = iv(l2+l3)
            else
               n1 = m1
               l1 = k1
               n2 = m2
               l2 = k2
               iph = iv(l1+ilam)
            endif
         else
            n3 = m3
            l3 = k3
            n4 = m4
            l4 = k4
            if (nn1.gt.nn2) then
               n1 = m2
               l1 = k2
               n2 = m1
               l2 = k1
               iph = iv(l3+ilam)
            else
               n1 = m1
               l1 = k1
               n2 = m2
               l2 = k2
               iph = +1
            endif
         endif
      endif
c
c---------------------------------------------------------------
      nn1 = 2*n1+l1
      nn2 = 2*n2+l2
      nn3 = 2*n3+l3
      nn4 = 2*n4+l4
      prout = sqh(l1)*wf(n1)*wg(n1+l1+1)*sqh(l2)*wf(n2)*wg(n2+l2+1)*
     &        sqh(l3)*wf(n3)*wg(n3+l3+1)*sqh(l4)*wf(n4)*wg(n4+l4+1)
      xf = sqh(0)**ichi
      iat = min0(nn1,nn3)
      if (l1.gt.lx.or.l2.gt.lx) goto 910
      if (iat.gt.lx)  goto 910
c
c---- main loop
      t = zero
      do 10 la = 0,iat
         lbl = iabs(la-l3)
         ibl = lbl+mod(la+lbl+l3,2)
         lbt = min0(nn2,nn3-la,la+l3)
         ibt = lbt-mod(la+lbt+l3,2)
         if (ibl.gt.ibt) goto 10
         if (ibt.gt.lx)  goto 910
c-------------
         do 20 lb = ibl,ibt,2
            lcl = iabs(la-l1)
            icl = lcl+mod(la+lcl+l1,2)
            lct = min0(nn1-la,nn4,la+l1)
            ict = lct-mod(la+lct+l1,2)
            if (icl.gt.ict) goto 20
            lo  = locs(la,lb,l3)
            c3  = wig(lo)
            if (ict.gt.lx) goto 910
c----------------
            do 30 lc = icl,ict,2
               ldl = max0(iabs(lb-l2),iabs(lc-l4))
               idl = ldl+mod(lb+ldl+l2,2)
               ldt = min0(nn2-lb,nn4-lc,lb+l2,lc+l4)
               idt = ldt-mod(lb+ldt+l2,2)
               if (idl.gt.idt) goto 30
               lo = locs(la,lc,l1)
               c1 = wig(lo)
               if (idt.gt.lx)  goto 910
c-------------------
               do 40 ld = idl,idt,2
                  mm1 = nn4-lc-ld
                  mm2 = nn2-lb-ld
                  ndt = min0(mm1/2,mm2/2)
                  if (ndt.lt.0) goto 40
                  iphd = iv(ld)
                  lo = locs(lb,ld,l2)
                  c2 = wig(lo)
                  lo = locs(lc,ld,l4)
                  c4 = wig(lo)
                  if (min0(ilam,la).gt.l3) then
                     qj = sninj0(l1,lc,la,l2,ld,lb,ilam,l4,l3)
                     qj = qj*iv(l1+l2+ilam)
                  else
                     qj = sninj0(la,lc,l1,lb,ld,l2,l3,l4,ilam)
                  endif
c----------------------
                  do 50 nd = 0,ndt
                     nnd = 2*nd+ld
                     na = (nn3-nn2+nnd-la)/2
                     if (na.lt.0) goto 50 
                     nb = (nn2-nnd-lb)/2
                     if (nb.lt.0) goto 50 
                     nc = (nn4-nnd-lc)/2
                     if (nc.lt.0) goto 50 
                     ll = (2*la+1)*(2*lb+1)*(2*lc+1)*(2*ld+1)
                     gg = fi(na)*gmi(na+la+1)*fi(nb)*gmi(nb+lb+1)*
     &                    fi(nc)*gmi(nc+lc+1)*fi(nd)*gmi(nd+ld+1)
                     t  = t + iphd*ll*c1*c2*c3*c4*qj*gg*prout
   50             continue
   40          continue
   30       continue
   20    continue
   10 continue
      talmos = iph*t*xf/pi
      if (dabs(talmos).lt.1.0d-10) talmos = zero
      return
c
 900  stop  'in TALMOS: n_j too large '
 910  stop ' in TALMOS: lx  too small '
c-end-TALMOS
      end
C=======================================================================

      function sninj0(j11,j12,j13,j21,j22,j23,j31,j32,j33)

C=======================================================================
C
C     calculates 9j-Symbol for integer j-values
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c 
      sninj0 = 0.0
c
c---- check for triangel rules (can be skipped if know before) 
      if (j12.lt.abs(j11-j13).or.j12.gt.j11+j13) return
      if (j22.lt.abs(j21-j23).or.j22.gt.j21+j23) return
      if (j32.lt.abs(j31-j33).or.j32.gt.j31+j33) return
      if (j21.lt.abs(j11-j31).or.j21.gt.j11+j31) return
      if (j22.lt.abs(j12-j32).or.j22.gt.j12+j32) return
      if (j23.lt.abs(j13-j33).or.j23.gt.j13+j33) return
c--------------------------------------------------------------
c
      k1 = max0(iabs(j21-j32),iabs(j11-j33),iabs(j12-j23))
      k2 = min0(j21+j32,j12+j23,j11+j33)
      do k = k1,k2
         sninj0 = sninj0 + (k+k+1)*
     &            racah0(j32,j21,k,j11,j33,j31)*
     &            racah0(j12,j22,j32,j21,k,j23)*
     &            racah0(k,j23,j12,j13,j11,j33)
      enddo
c
      return
c-end-SNINJ0
      end
C=======================================================================

      double precision function racah0(j1,j2,j3,l1,l2,l3)

C=======================================================================
C
C     Calculates 6j-symbol (notation of Edmonds) for integer ang.momenta
C
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      include 'dirhb.par'
c
      common /gfviv / iv(0:igfv)
      common /gfvfak/ fak(0:igfv)
      common /gfvfi / fi(0:igfv)
      common /gfvwf / wf(0:igfv)
      common /gfvwfi/ wfi(0:igfv)
c
      dsq(i,k,l) = wf(i+k-l)*wfi(i+k+l+1)*wf(i-k+l)*wf(k+l-i)
c
      racah0 = 0.
      i1 = j1+j2+j3
      i2 = j1+l2+l3
      i3 = l1+j2+l3
      i4 = l1+l2+j3
      i5 = j1+j2+l1+l2
      i6 = j2+j3+l2+l3
      i7 = j3+j1+l3+l1
      n1 = max0(i1,i2,i3,i4)
      n2 = min0(i5,i6,i7)
      if (n1.gt.n2) return
      do 10 n = n1,n2
   10 racah0 = racah0+iv(n)*fak(n+1)*
     &         fi(n-i1)*fi(n-i2)*fi(n-i3)*fi(n-i4)*
     &         fi(i5-n)*fi(i6-n)*fi(i7-n)
   20 racah0 = dsq(j1,j2,j3)*dsq(j1,l2,l3)*
     &         dsq(l1,j2,l3)*dsq(l1,l2,j3)*racah0
c
      return
c-end-RACAH0
      end

c======================================================================c

      subroutine singf(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions       
c     in the spherical oscillator basis
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
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /single/ sp(nfgx,nbx)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN SINGF ********************************'
c
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
c
c        SIGMA*P
c----------------
         call sigp(nf,ng,ib,sp(1,ib),lpr)
c
      enddo
C
      if (lpr)
     &write(l6,*) ' ****** END SINGF **********************************'
      return
c-end-SINGF
      end  
c=====================================================================c

      subroutine sigp(nf,ng,ib,aa,lpr)

c=====================================================================c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
c
      character tb*5                                            ! blokap
      character tt*8                                            ! bloqua
c
      dimension aa(ng,nf)
c
      common /baspar/ hom,hb0,b0
      common /blokap/ nb,kb(nbx),nrbl(nbx,4),mb(nbx),tb(nbx)
      common /bloosc/ ia(nbx,2),id(nbx,2)
      common /bloqua/ nt,nr(ntx),nl(ntx),nj(ntx),kk(ntx),np(ntx),tt(ntx)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /radosc/ rnl(1:nrx,0:nlx,0:ngh),rnl1(1:nrx,0:nlx,0:ngh)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      kappa = kb(ib)
      lf    = lfkap(kappa)
      lg    = lgkap(kappa)
c
      kkk = - kappa - 1
      do n2 = 1,nf
      do n1 = 1,ng
         s = zero
         do ih = 1,ngh
            s = s + rnl(n1,lg,ih) * 
     &              ( - rnl1(n2,lf,ih) + kkk*rnl(n2,lf,ih)/xh(ih))    
         enddo
         aa(n1,n2) = s
      enddo
      enddo

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
      if (yp0.gt. 0.999d30) then
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
c
c     initializes potentials
c     inin = 0:   reads fields from tape lwin
c            1:   saxon-woods
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      include 'dirhb.par'
c
      logical lpr
      character*2 nucnam
c
      dimension vso(2),r0v(2),av(2),rso(2),aso(2)
c
      common /baspar/ hom,hb0,b0
      common /coulmb/ cou(0:ngh),drvp(0:ngh)
      common /gaussh/ xh(0:ngh),wh(0:ngh),ph(0:ngh),rb(0:ngh)
      common /initia/ inin,inink
      common /masses/ amu,amsig,amome,amdel,amrho
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nama,npr(2),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /potpot/ vps(0:ngh,2),vms(0:ngh,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c=======================================================================
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
      data v0/-71.28/,akv/0.4616/
      data r0v/1.2334,1.2496/,av/0.615,0.6124/
      data vso/11.1175,8.9698/
      data rso/1.1443,1.1401/,aso/0.6476,0.6469/
c---- potentials are read in INOUT
      if (inin.eq.0) return
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN START ********************************'
c
c
      if (lpr) then
         write(l6,'(a,f10.4)') ' v0     = ',v0
         write(l6,'(a,f10.4)') ' kappa  = ',akv
         write(l6,'(a,2f10.4)') ' lambda = ',vso
         write(l6,'(a,2f10.4)') ' r0     = ',r0v
         write(l6,'(a,2f10.4)') ' a      = ',av
         write(l6,'(a,2f10.4)') ' r0-so  = ',rso
         write(l6,'(a,2f10.4)') ' a-so   = ',aso
      endif
      do ih = 0,ngh
         r = rb(ih)
c
c------- Woods-Saxon potential
         do it = 1,itx
            ita = 3-it
            rav = r0v(it)*amas**third
            rao = rso(it)*amas**third
            vp  = v0*(one - akv*(npr(it)-npr(ita))/amas)
            vls = vp * vso(it)
c
            argv = (r - rav)/av(it)
            if (argv.le.65.d0) then
               u = vp/(one + exp(argv))
            else
               u = zero
            endif
            argo = (r - rao)/aso(it)
            if (argo.le.65.d0) then
               w = -vls/(one + exp(argo))
            else
               w = zero
            endif
            vps(ih,it) = u 
            vms(ih,it) = w
         enddo   ! it
	     if (itx.eq.1) then
	        vms(ih,2) = vms(ih,1)
	        vps(ih,2) = vps(ih,1)
         endif   ! itx=1
c
c------- Coulomb potential
         cou(ih) = zero
         if (icou.ne.0) then
            rc = r0v(2)*amas**third
            if (r.lt.rc) then
               c = half*(3/rc - r*r/rc**3)
            else
               c = one/r
            endif
            cou(ih)   = c*npr(2)/alphi
	        vps(ih,2) = vps(ih,2) + cou(ih)*hbc
	        vms(ih,2) = vms(ih,2) + cou(ih)*hbc
	     endif   ! icou.ne.0
      enddo   ! ih
      if (lpr)
     &write(l6,'(/,a)') ' Initial potentials of Saxon-Woods shape '

      if (lpr)
     &write(l6,*) ' ****** END START **********************************'
      return
c-end START
      end

