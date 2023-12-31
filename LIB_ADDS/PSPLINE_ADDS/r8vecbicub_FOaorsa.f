      subroutine r8vecbicub_FOaorsa(ictA,ict,ivec,xvec,yvec,ivd,fvalA,
     >     fvalBXRe,fvalBYRe,fvalBZRe,fvalBXIm,fvalBYIm,fvalBZIm,
     >     fvalEXRe,fvalEYRe,fvalEZRe,fvalEXIm,fvalEYIm,fvalEZIm,
     >     nx,xpkg,ny,ypkg,fsplA,fsplBXRe,fsplBYRe,fsplBZRe,
     >     fsplBXIm,fsplBYIm,fsplBZIm,fsplEXRe,fsplEYRe,fsplEZRe,
     >     fsplEXIm,fsplEYIm,fsplEZIm,inf2,iwarn,ier)
c
c  vectorized spline evaluation routine -- 2d *compact* spline
c  1.  call vectorized zone lookup routine
c  2.  call vectorized spline evaluation routine
c
c--------------------------
c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iwarn1,iwarn2
!============
      integer ictA(6)          ! selector:
      integer ict(6)          ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for d2f/dx2 (don't evaluate if ict(4)=0)
c        ict(5)=1 for d2f/dy2 (don't evaluate if ict(5)=0)
c        ict(6)=1 for d2f/dxdy (don't evaluate if ict(6)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c
c  list of (x,y) pairs:
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
      REAL*8 yvec(ivec)                   ! y-locations at which to evaluate
c
      integer ivd                       ! 1st dimension of output array
c
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      REAL*8 fvalA(ivd,*)       ! output array
      REAL*8 fvalBXRe(ivd,*)    ! output array
      REAL*8 fvalBYRe(ivd,*)    ! output array
      REAL*8 fvalBZRe(ivd,*)    ! output array
      REAL*8 fvalBXIm(ivd,*)    ! output array
      REAL*8 fvalBYIm(ivd,*)    ! output array
      REAL*8 fvalBZIm(ivd,*)    ! output array
      REAL*8 fvalEXRe(ivd,*)    ! output array
      REAL*8 fvalEYRe(ivd,*)    ! output array
      REAL*8 fvalEZRe(ivd,*)    ! output array
      REAL*8 fvalEXIm(ivd,*)    ! output array
      REAL*8 fvalEYIm(ivd,*)    ! output array
      REAL*8 fvalEZIm(ivd,*)    ! output array
c
c  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
c  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
c   --etc--
c
c input:
      integer nx,ny                     ! dimension of spline grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      integer inf2                      ! fspl 3rd array dimension, .ge.nx
      REAL*8 fsplA(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBXRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBYRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBZRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBXIm(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBYIm(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplBZIm(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEXRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEYRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEZRe(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEXIm(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEYIm(0:3,inf2,ny) ! (compact) spline coefficients
      REAL*8 fsplEZIm(0:3,inf2,ny) ! (compact) spline coefficients
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer ix(ivec)                  ! zone indices {j}
      REAL*8 dxn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hx(ivec)                     ! h(j) vector
      REAL*8 hxi(ivec)                    ! 1/h(j) vector
c
      integer iy(ivec)                  ! zone indices {j}
      REAL*8 dyn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hy(ivec)                     ! h(j) vector
      REAL*8 hyi(ivec)                    ! 1/h(j) vector
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(nx.lt.2) then
         write(6,*) ' ?vecbicub:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ny.lt.2) then
         write(6,*) ' ?vecbicub:  ny.lt.2:  ny = ',ny
         ier=1
      endif
c
      if(ivec.le.0) then
         write(6,*) ' ?vecbicub:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecbicub:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      ix=0
      iy=0
      call r8xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call r8xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      iwarn=iwarn1+iwarn2
c
c  vectorized evaluation
c
      call r8fvbicub(ictA,ivec,ivd,fvalA,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplA,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBXRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBXRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBYRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBYRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBZRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBZRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBXIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBXIm,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBYIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBYIm,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalBZIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplBZIm,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEXRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEXRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEYRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEYRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEZRe,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEZRe,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEXIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEXIm,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEYIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEYIm,inf2,ny)
      call r8fvbicub(ict,ivec,ivd,fvalEZIm,ix,iy,dxn,dyn,
     >     hx,hxi,hy,hyi,fsplEZIm,inf2,ny)
c
      return
      end
