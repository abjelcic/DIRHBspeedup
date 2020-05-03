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
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)   
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
      subroutine centmas
c======================================================================c
      implicit real*8 (a-h,o-z)
      include 'dirhb.par'
      
      call cmcd
      call cmcn
      
      return
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
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)         

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
      common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(2)         

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

c======================================================================c
c
c      PROGRAM DIRRHB-triaxial
c
c======================================================================c
c
c     Relativistic Hartree-Bogoliubov theory in a triaxial basis
c     Main part of the code
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
c
c
c---- sets data
      call default()
c
c---- reads in data     
      call reader(.true.)
c
c---- force-parameters
      call forces(.true.)
c
c---- Gauss-Hermite mesh points
      call gaush(.false.)
c
c---- oscillator basis for single particle states
      call base(.false.)
c
c---- preparations
      call prep(.true.)
c
c---- initialization of the potentials
      call inout(1,.false.)
      call dinout(1,.false.)
      
      call start(.false.)
c
c---- wavefunctions at Gauss-Meshpoints
      call gaupol(.false.)
c
c---- single-particle matix elements
      call singf(.false.)
c      
c---- pairing matrix elements
      call singd(.false.)      
c
c---- meson propagators
      call greemes(.false.)
      call greecou(.false.)    

c
c---- iteration
      call iter(.true.)
      
c---- transformation to the canonical basis
      call canon(.true.)
      
c---- center-of-mass correction      
      call centmas  
c
c---- results
      call resu(.true.)
      call plot(.false.)
c
c---- punching of potentials to tape  
      call inout(2,.false.)
      call dinout(2,.false.)
c
c
      stop ' FINAL STOP OF DIRHBT'
c-end-DIRHBT
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

