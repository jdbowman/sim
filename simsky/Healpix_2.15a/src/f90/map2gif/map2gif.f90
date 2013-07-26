!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
PROGRAM MAP2GIF
   !-------------------------------------------------------------------
   ! This program is a simplified kluge of mollview.pro, a basic tool to
   ! view a Mollweide projection of maps made in the Healpix pixelisation
   ! scheme. A map in FITS format is read and converted to a GIF image.
   !
   ! PURPOSE:
   !    basic tool to view a Mollweide projection of maps made
   !    in Healpix pixelisation
   !
   ! CALLING SEQUENCE:
   !    MAP2GIF [-inp input_file] [-out output_file]'
   !            [-max usermax] [-min usermin]'
   !            [-add offset] [-mul factor]'
   !            [-pro projection] [-res gridresn]'
   !            [-lon lon0] [-lat lat0]'
   !            [-col color_table] [-sig signal]'
   !            [-xsz xsize] [-log logflg]'
   !            [-ash ashflg]'
   !            [-bar barflg] [-ttl title]'
   !            [-hlp]'
   !
   ! INPUTS:
   !    input_file = string containing the name of a FITS file of
   !                 the healpix map in an extension or in the image field
   !
   !    output_file= string containing the name of the output GIF file(s)
   !
   !    color_table= color table number, in [1,5]
   !                 if unset the color table will be 5 (Blue-Red)
   !
   !    signal     = integer defining which field of the FITS file should be plotted
   !                 for an IQU map    1 = I; 2 = Q; 3 = U
   !                 for a polarized WMAP map  4=5= N_obs, 6 = N_QQ, ...
   !
   !    xsize      = x-dimension of image array to be stored in GIF file
   !                 if unset the x-size will be ysize = xsize/2
   !
   !     offset    = displays   map = (data+offset) * factor, default = 0
   !
   !     factor    = displays   map = (data+offset) * factor, default = 1
   !
   !     logflg    = if .true.,  displays    log(map)
   !                   default = .false. (ie. use actual temperatures)
   !
   !     ashflg    = if .true.,  displays   sinh^(-1)(map) (inverse Hyperbolic sinus)
   !                   default = .false. (ie. use actual temperatures)
   !
   !     title     = character string to be used as plot title
   !
   !     usermax   = sets the upper bound of plotted map
   !
   !     usermin   = sets the lower bound of plotted map
   !
   !     lon0      = longitude (deg) of central point
   !
   !     lat0      = latitude (deg) of central point
   !
   !     barflg    =   if set, add color bar
   !
   !     projection =  either "MOL" or "mol" (Mollweide, default) 
   !                       or "GNO" or "gno" (gnomic)
   !
   !     gridresn  = resolution in (u,v) plan in arcmin (only for Gnomic projection)
   !
   !
   ! PROCEDURES USED:
   !    parseOptions, getnumext_fits, getsize_fits, input_map, proj_mollw, proj_gnomic
   !
   ! EXTERNAL LIBRARY:
   !    this code uses the CFITSIO library that can be found at
   !    http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
   !
   ! RELATED LITERATURE:
   !    see the web site http://www.tac.dk/~healpix
   !
   ! MODIFICATION HISTORY:
   !    October 1997,  Eric Hivon,      TAC (original IDL code)
   !    December 1998, A.J. Banday,     MPA (simple conversion to F90)
   !    December 1998, M.S. Bartelmann, MPA (cleanup, color tables, ...)
   !    January 1999,  A.J. Banday,     MPA (added gnomic routine)
   !    ??             M. Reinecke,     MPA
   !    November 2004, Eric Hivon,      IPAC (upgraded for WMAP polarized maps)
   !    Jan 2005,      EH               IPAC (edited to reduce stack consomption)
   !
   !-------------------------------------------------------------------

   USE healpix_types
   USE misc_utils
   USE fitstools
   USE gifmod

   IMPLICIT NONE

   INTEGER(I4B) :: nsmax
   INTEGER(I4B) :: npixtot, nmaps
   INTEGER(I4B) :: ordering
   INTEGER(I4B) :: i

   REAL(SP),     DIMENSION(:,:),   ALLOCATABLE :: map_IQU
   REAL(SP),     DIMENSION(:),     ALLOCATABLE :: map
   REAL(SP),     DIMENSION(:,:),   ALLOCATABLE :: image

   INTEGER(I4B), ALLOCATABLE, DIMENSION(:,:) :: imgint

   INTEGER(I4B), POINTER, DIMENSION(:,:) :: imgbar
   INTEGER(I4B), POINTER, DIMENSION(:,:) :: imgttl

   LOGICAL(LGT), DIMENSION(:,:),   ALLOCATABLE :: mask
   LOGICAL(LGT), DIMENSION(:),     ALLOCATABLE :: observed

   REAL(SP)     :: mindata, maxdata, Tmin, Tmax

   INTEGER(I4B) :: status

   CHARACTER(LEN=7), PARAMETER :: code = 'MAP2GIF'
   character(len=*), parameter :: VERSION = HEALPIX_VERSION

   CHARACTER(LEN=filenamelen) :: input_file = ''
   CHARACTER(LEN=filenamelen) :: output_file = ''
   CHARACTER(LEN=filenamelen) :: title = ''
   CHARACTER(LEN=3)   :: projection = 'MOL'
   REAL(SP)           :: usermax = 1.1e30
   REAL(SP)           :: usermin = -1.1e30
   REAL(SP)           :: lon0 = 0.0
   REAL(SP)           :: lat0 = 0.0
   REAL(SP)           :: gridresn = 2.0
   INTEGER(I4B)       :: color_table = 5
   INTEGER(I4B)       :: signal = 1         ! I/intensity measurement
   INTEGER(I4B)       :: xsize  = 800
   LOGICAL(LGT)       :: logflg = .false.
   LOGICAL(LGT)       :: bar = .false.
   LOGICAL(LGT)       :: ashflag = .false.
   REAL(SP)           :: offset = 0.0
   REAL(SP)           :: factor = 1.0

   integer(i4b)       :: n_ext, i_ext, nmaps_sum, icolumn
   !-------------------------------------------------------------------

   !--- acquire input parameters ---
   call parseOptions()

   !--- read in the FITS file header, find out the pixel number of the map
   !    and its ordering ---
   nmaps_sum = 0
   icolumn   = signal
   n_ext   = getnumext_fits(input_file)
   do i_ext = 0, n_ext-1
      npixtot = getsize_fits(input_file, nmaps = nmaps, ordering=ordering, extno=i_ext)
      nmaps_sum = nmaps_sum + nmaps
      if (signal <= nmaps_sum) exit
      icolumn = signal - nmaps_sum
   enddo

   !--- are requested signal and map compatible? ---
   if (signal > nmaps_sum) THEN
      print *,' requested signal incompatible with file'
      call fatal_error
   end if

   nsmax = NINT(SQRT( npixtot / 12.) )
   call assert(npixtot == 12*nsmax*nsmax,'wrong number of pixels')

   !--- allocate space for arrays ---
   ALLOCATE(map(0:npixtot-1),stat = status)
   call assert_alloc(status,'map2gif','map') 
   ALLOCATE(map_IQU(0:npixtot-1,1:nmaps),stat = status)
   call assert_alloc(status,'map2gif','map_IQU') 

   !--- read in the FITS file ---
   call input_map(input_file, map_IQU, npixtot, nmaps, fmissval=-1.6375e30_sp, extno=i_ext)

   !--- select correct signal to plot ---
   map(:) = map_IQU(:,icolumn)

   deallocate(map_IQU,stat = status)
   allocate(observed(0:npixtot-1),stat = status)
   call assert_alloc(status,'map2gif','observed') 

   !--- Test for unobserved pixels ---
   observed = map > -1.6375e30
   !!where (observed) map = (map+offset)*factor ! << requires large stack size for large maps
   do i = 0, npixtot-1
      if (observed(i)) map(i) = (map(i)+offset)*factor
   enddo
   if (ashflag) then
      do i=0,npixtot-1
        if (observed(i)) then
          if (map(i)>=0) then
            map(i) = log(map(i)+sqrt(map(i)**2+1))
          else
            map(i) = -log(-map(i)+sqrt(map(i)**2+1))
          endif
        endif
      end do
   endif
   IF( logflg )THEN
      observed = map > 0.
!       where(observed)
!         map=log10(map)
!       elsewhere
!         map=-1.6375e30
!       endwhere
      do i=0, npixtot - 1
         if (observed(i)) then
            map(i) = log10(map(i))
         else
            map(i) = -1.6375e30
         endif
      enddo
   END IF
   mindata = minval(map, mask = observed)
   maxdata = maxval(map, mask = observed)

   if (mindata==maxdata) then
     print *,"Warning: map contains completely uniform values!"
     maxdata=mindata+1e-5
   endif

   !--- sets MIN and MAX using linear scaling ---
   if (ABS((maxdata+mindata)/(maxdata-mindata)) < 5.e-2) then
      !--- center of color scale is 0
      Tmax = MAXVAL(ABS((/mindata,maxdata/)))
      Tmin = -Tmax
   else
      Tmax = maxdata
      Tmin = mindata
   end if

   !--- allow for user min/max settings ---
   if (usermax <  1.1e30) Tmax = usermax
   if (usermin > -1.1e30) Tmin = usermin

   !--- adjust 'offscale' temperatures ---
   !! WHERE(observed) map = max(min(map,Tmax),Tmin) ! << requires large stack size for large maps
   do i = 0, npixtot - 1
      if (observed(i)) map(i) = max(min(map(i), Tmax), Tmin)
   enddo
   deallocate(observed)

   !--- create the projection ---
   select case (projection)
   case ('MOL', 'mol')
      !--- allocate space for image and mask ---
      ALLOCATE(image(0:xsize-1,0:xsize/2-1),stat = status)
      call assert_alloc(status,'map2gif','image') 
      ALLOCATE(mask(0:xsize-1,0:xsize/2-1),stat = status)
      call assert_alloc(status,'map2gif','mask') 

      !--- create the mollweide projection ---
      call proj_mollw(&
           map,      &
           nsmax,    &
           ordering, &
           xsize,    &
           lon0,     &
           image,    &
           mask      &
           )
   case ('GNO', 'gno')
      !--- allocate space for image and mask ---
      ALLOCATE(image(0:xsize-1,0:xsize-1),stat = status)
      call assert_alloc(status,'map2gif','image') 
      ALLOCATE(mask(0:xsize-1,0:xsize-1),stat = status)
      call assert_alloc(status,'map2gif','mask') 

      !--- create the gnomic projection ---
      call proj_gnomic(&
           map,        &
           nsmax,      &
           ordering,   &
           xsize,      &
           lon0,       &
           lat0,       &
           gridresn,   &
           image,      &
           mask        &
           )
   case default
      print '(" unknown projection type: ",a3)', projection
      call fatal_error
   end select

   !--- load the color table ---
   call setcol( color_table )

   !--- scale image ---
   allocate(imgint(size(image,1),size(image,2)))
   call imgscl(&
        image,  &
        imgint, &
        mask    &
        )
   deallocate(image)

   !--- add color bar ---
   if (bar) then
      call addbar(&
           imgint, &
           imgbar  &
           )
      deallocate(imgint)
      allocate(imgint(size(imgbar,1),size(imgbar,2)))
      imgint = imgbar
      deallocate(imgbar)
   end if

   !--- add title ---
   if (title /= '') then
      call addttl(&
           imgint, &
           imgttl, &
           title   &
           )
      deallocate(imgint)
      allocate(imgint(size(imgttl,1),size(imgttl,2)))
      imgint = imgttl
      deallocate(imgttl)
   end if

   !--- create the GIF image ---
   call gifmap(&
        imgint,           &
        trim(output_file) &
        )

   !--- deallocate memory for arrays ---
   DEALLOCATE(map,stat = status)
   DEALLOCATE(mask,stat = status)

   !--- Exit, stage left ...

contains

subroutine euler_matrix(       &
                        a1,    &
                        a2,    &
                        a3,    &
                        matrix &
                                )

   ! --------------------------------------------------------------
   !
   ! this subroutine computes the Euler matrix corresponding
   ! to aeronautics conventionn following the prescription
   !    rotation a1 around original Z
   !    rotation a2 around interm   Y
   !    rotation a3 around final    X
   ! see H. Goldestein, Classical Mechanics (2nd Ed.) p. 147 for
   ! discussion
   ! INPUTS:
   !        a1, a2, a3 = Euler angles, in radians
   !        all the angles are measured counterclockwise
   !
   ! MODIFICATION HISTORY:
   !    October 1997, Eric Hivon, TAC (original IDL code)
   !    January 1999, A.J. Banday MPA (conversion to F90)
   !    March 22 1999, E.H. Caltech   (modification to match IDL version,
   !    correction of a bug on the matrix product)
   !
   ! --------------------------------------------------------------

   USE healpix_types

   IMPLICIT NONE

   REAL(DP)                 :: a1,a2,a3, c1,c2,c3, s1,s2,s3
   REAL(DP)                 :: ze,un
   REAL(DP), DIMENSION(3,3) :: matrix
   REAL(DP), DIMENSION(3,3) :: m1,m2,m3

   ! --------------------------------------------------------------

   !-- convert to sin and cosine values --
   c1 = COS(a1)
   s1 = SIN(a1)
   c2 = COS(a2)
   s2 = SIN(a2)
   c3 = COS(a3)
   s3 = SIN(a3)

   ze = 0.0_dp
   un = 1.0_dp

   !-- form full rotation matrix --
   !-- FOR REFERENCE:
   ! m1 = [[ c1, s1,  0],[-s1, c1,  0],[  0,  0,  1]] ! around orig  z
   ! m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ! around orig  y
   ! m3 = [[  1,  0,  0],[  0, c3, s3],[  0,-s3, c3]] ! around orig  x
   ! A ## B = matrix product AB
   ! matrix = m3 # ( m2 # m1)  = m1 m2 m3 ! here was the bug

   m1(:,1) = (/ c1, s1, ze /)
   m1(:,2) = (/-s1, c1, ze /)
   m1(:,3) = (/ ze, ze, un /)

   m2(:,1) = (/ c2, ze, s2 /)
   m2(:,2) = (/ ze, un, ze /)
   m2(:,3) = (/-s2, ze, c2 /)

   m3(:,1) = (/ un, ze, ze /)
   m3(:,2) = (/ ze, c3, s3 /)
   m3(:,3) = (/ ze,-s3, c3 /)

   matrix = matmul(m1,matmul(m2,m3))

end subroutine euler_matrix

SUBROUTINE PROJ_MOLLW(&
     map,      &
     nside,    &
     ordering, &
     xsize,    &
     lon0,     &
     image,    &
     mask      &
     )
   ! --------------------------------------------------------------
   ! This subroutine takes an input HEALPIX sky map array and
   ! generates the (u,v) position on a mollweide projection.
   !
   ! INPUTS:
   !       map     = array containing pixel temperatures
   !       nside   = defines number of pixels on the sky
   !       ordering= whether input map is in RING or NESTED format
   !       xsize   = x-dimension of output image array
   !       lon0    = longitude (in degrees) of the center of the plot
   !       image   = output mollweide projection image
   !       mask    = logical array defining usable pixels in output image array
   !
   ! PROCEDURES USED:
   !            ang_pix, ang_nest_pix
   !
   !  RELATED LITERATURE:
   !  see the web site http://www.tac.dk/~healpix
   !
   ! MODIFICATION HISTORY:
   !    October 1997, Eric Hivon, TAC (original IDL code)
   !       December 1998, A.J. Banday MPA (conversion to F90)
   !
   ! --------------------------------------------------------------

   USE healpix_types
   USE pix_tools

   IMPLICIT NONE

   INTEGER(I4B) :: I,J
   INTEGER(I4B) :: id_pix
   INTEGER(I4B) :: ordering
   INTEGER(I4B) :: nside

   INTEGER(I4B) :: xsize, ysize

   REAL(SP), DIMENSION(0:12*nside**2-1) :: map

   REAL(SP),     DIMENSION(0:xsize-1,0:xsize/2-1) :: image
   LOGICAL(LGT), DIMENSION(0:xsize-1,0:xsize/2-1) :: mask

   REAL(SP) :: lon0
   REAL(DP) :: lon0rad
   REAL(DP) :: DtoR = PI/180.0_dp

   INTEGER(I4B) :: xc, dx, yc, dy

   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: u, v
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: out_lat, out_lon

   INTEGER(I4B) :: status

   ! --------------------------------------------------------------

   ysize = xsize/2

   !--- generates the (u,v) position on the mollweide map ---
   image(:,:) = 1.                         ! white background ???

   xc = (xsize-1)/2
   dx = xc
   yc = (ysize-1)/2
   dy = yc

   !--- allocate memory for arrays ---
   ALLOCATE(u(0:xsize-1,0:ysize-1),stat = status)
   call assert_alloc(status,'map2gif','u')
   do i = 0,xsize-1
      u(i,:) = real(i,dp)
   end do
   ALLOCATE(v(0:xsize-1,0:ysize-1),stat = status)
   call assert_alloc(status,'map2gif','v')
   do i = 0,ysize-1
      v(:,i) = real(i,dp)
   end do
   ALLOCATE(out_lat(0:xsize-1,0:ysize-1),stat = status)
   call assert_alloc(status,'map2gif','out_lat')
   ALLOCATE(out_lon(0:xsize-1,0:ysize-1),stat = status)
   call assert_alloc(status,'map2gif','out_lon')

   u =  2._dp*(u - real(xc,dp))/(real(dx,dp)/1.02_dp)   ! 1.02 : fudge factor in [-2,2]*1.02
   v =    (v - real(yc,dp))/(real(dy,dp)/1.02_dp)   ! in [-1,1] * 1.02

   !--- for each point on the mollweide map looks for the corresponding
   !    position (lon, lat) on the sphere ---
   mask = (u**2/4._dp + v**2) <= 1._dp
   out_lat = 1.e12_dp ! values for points out of the sphere
   out_lon = 1.e12_dp
   lon0rad = real(lon0,dp) * DtoR
   WHERE( mask )
      out_lat = PI/2._dp - (ASIN ( 2._dp/PI * (&
           & ASIN(v) + v*SQRT( (1._dp-v)*(1._dp+v) )) ))
      ! colat in [0,pi]
      out_lon = -lon0rad - PI/2._dp * u/MAX(SQRT( (1._dp-v)*(1._dp+v) ),1.e-6_dp)
      ! lon in [-pi,pi], the minus sign for astro convention
   END WHERE

   WHERE(out_lon < 0._dp)
      out_lon = out_lon  + 2._dp*PI
      ! lon in RAD in [0,2pi]
   END WHERE

   !--- converts the position on the sphere into pixel number and project the
   !    corresponding data value on the map ---
   DO I = 0,xsize-1
      DO J = 0,ysize-1
         IF(ABS(out_lat(I,J)) .LE. 1.e5_dp)THEN ! keep only meaningful points
            if (ordering .eq. 1) then
               call ang2pix_ring( nside, out_lat(I,J), out_lon(I,J), id_pix )
               ! id_pix in [0,12*nside^2-1]
            else
               call ang2pix_nest( nside, out_lat(I,J), out_lon(I,J), id_pix )
               ! id_pix in [0,12*nside^2-1]
            end if
            image(I,J) = map(id_pix)
         END IF
      END DO
   END DO

   !--- deallocate memory for arrays ---
   deallocate( u )
   deallocate( v )
   deallocate( out_lat )
   deallocate( out_lon )

END SUBROUTINE PROJ_MOLLW

SUBROUTINE PROJ_GNOMIC(&
     map,      &
     nside,    &
     ordering, &
     xsize,    &
     lon0,     &
     lat0,     &
     resgrid,  &
     image,    &
     mask      &
     )

   ! --------------------------------------------------------------
   ! This subroutine takes an input HEALPIX sky map array and
   ! generates the gnomic (tangential) projection on a rectangular grid.
   !
   ! INPUTS:
   !       map     = array containing pixel temperatures
   !       nside   = defines number of pixels on the sky
   !       ordering= whether input map is in RING or NESTED format
   !       xsize   = x-dimension of output image array
   !       lon0    = longitude (in degrees) of the center of the plot
   !       lat0    = latitude (in degrees) of the center of the plot
   !       gridresn= resolution of the grid elements in ARCMINUTES
   !       image   = output gnomic projection image
   !       mask    = logical array defining usable pixels in output
   !                 image array
   !
   ! PROCEDURES USED:
   !            ang_pix, ang_nest_pix, euler_matrix
   !
   !  RELATED LITERATURE:
   !  see the web site http://www.tac.dk/~healpix
   !
   ! MODIFICATION HISTORY:
   !    October 1997, Eric Hivon, TAC (original IDL code)
   !    January 1999, A.J. Banday MPA (conversion to F90)
   !    March   1999, EH Caltech (consistent with new IDL version,
   !                       fixed a bug on rotation)
   !
   ! --------------------------------------------------------------

   USE healpix_types
   USE pix_tools

   IMPLICIT NONE

   INTEGER(I4B) :: I,J
   INTEGER(I4B) :: id_pix
   INTEGER(I4B) :: ordering
   INTEGER(I4B) :: nside

   INTEGER(I4B) :: xsize

   REAL(SP), DIMENSION(0:12*nside**2-1) :: map

   REAL(SP),     DIMENSION(0:xsize-1,0:xsize-1) :: image
   LOGICAL(LGT), DIMENSION(0:xsize-1,0:xsize-1) :: mask

   REAL(SP) :: lon0, lat0, resgrid
   REAL(DP) :: lon0rad, lat0rad, resgrid_rad
   REAL(DP) :: DtoR = PI/180.0_dp

   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xg, yg, zg, rg
   REAL(DP), DIMENSION(3,3)              :: mm
   REAL(DP), DIMENSION(3)                :: nvs

   INTEGER(I4B) :: status

   ! --------------------------------------------------------------

   !-- some basic variable manipulations --
   resgrid_rad = (real(resgrid,dp)/60.0_dp) * DtoR !/ real(xsize-1,dp)
   lon0rad = real(lon0,dp) * DtoR
   lat0rad = real(lat0,dp) * DtoR
!   theta0  = halfpi - lat0rad
!   phi0    = lon0rad

   !-- generates the projection around the chosen contact point
   !   on the planar grid  (1,Y,Z) --
   image(:,:) = 1.0
   mask(:,:)  = .true.

   !--- allocate memory for arrays ---
   ALLOCATE(xg(0:xsize-1,0:xsize-1),stat = status)
   call assert_alloc(status,'map2gif','xg')
   xg = 1.0_dp
   ALLOCATE(yg(0:xsize-1,0:xsize-1),stat = status)
   call assert_alloc(status,'map2gif','yg')
   do i = 0,xsize-1
      yg(i,:) = -(real(i,dp) - 0.5_dp*real(xsize,dp) + 0.5_dp) *&
           & resgrid_rad  ! minus sign = astro convention
   end do
   ALLOCATE(zg(0:xsize-1,0:xsize-1),stat = status)
   call assert_alloc(status,'map2gif','zg')
   do i = 0,xsize-1
      zg(:,i) = (real(i,dp) - 0.5_dp*real(xsize,dp) + 0.5_dp) *&
           & resgrid_rad
   end do
   ALLOCATE(rg(0:xsize-1,0:xsize-1),stat = status)
   call assert_alloc(status,'map2gif','rg')
   rg = sqrt(1.0_dp + yg**2 + zg**2)

   !-- the contact point is rotated by theta0 around y (=colatitude)
   !   and phi0 around z (=azimut) --
!   call euler_matrix(0.0_dp,theta0,phi0,mm) ! /y option used in IDL here
   call euler_matrix(lon0rad,lat0rad,0.0_dp,mm) ! /z option used in IDL here

   xg = 1.0_dp / rg
   yg = yg / rg
   zg = zg / rg


   !--- converts the position on the sphere into pixel number and project the
   !    corresponding data value on the map ---
   DO I = 0,xsize-1
      DO J = 0,xsize-1
         !-- rotated positions --
         nvs = matmul(mm,(/ xg(i,j), yg(i,j), zg(i,j) /))
         if (ordering .eq. 2) then
            call vec2pix_nest( nside, nvs, id_pix )
            ! id_pix in [0,12*nside^2-1]
         else ! default
            call vec2pix_ring( nside, nvs, id_pix )
            ! id_pix in [0,12*nside^2-1]
         end if
         image(I,J) = map(id_pix)
      END DO
   END DO

!   print*,mm
!   print*,matmul(mm,(/ 1, 0, 0 /))
!   print*,matmul(mm,(/ 0, 1, 0 /))
!   print*,matmul(mm,(/ 0, 0, 1 /))

!   call euler_matrix(1.0_dp,0.0_dp,0.0_dp,mm) ! /z option used in IDL here
!   print*,mm
!   print*,matmul(mm,(/ 1, 0, 0 /))
!   print*,matmul(mm,(/ 0, 1, 0 /))
!   print*,matmul(mm,(/ 0, 0, 1 /))


   !--- deallocate memory for arrays ---
   deallocate( xg )
   deallocate( yg )
   deallocate( zg )

END SUBROUTINE PROJ_GNOMIC

   subroutine parseOptions ()

     USE extension, ONLY : getArgument, nArguments
      implicit none

      integer(I4B)       :: i,n
      character(len=4)   :: opt
      character(len=filenamelen) :: arg

      n = nArguments()

      if (n < 4) then
         print *, 'usage: MAP2GIF [-inp input_file] [-out output_file]'
         print *, '               [-max usermax] [-min usermin]'
         print *, '               [-add offset] [-mul factor]'
         print *, '               [-pro projection] [-res gridresn]'
         print *, '               [-lon lon0] [-lat lat0]'
         print *, '               [-col color_table] [-sig signal]'
         print *, '               [-xsz xsize] [-log logflg]'
         print *, '               [-ash ashflg]'
         print *, '               [-bar add color bar] [-ttl title]'
         print *, '               [-hlp]'
         call fatal_error
      end if


      do i = 1,n,2
         call getArgument(i,opt)
         if (opt == '-hlp') then
            print 1, input_file,output_file,color_table,usermax,usermin,&
                 & signal,xsize,logflg,bar,title
            stop
         end if
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            call fatal_error
         end if
         call getArgument(i+1,arg)
         select case (opt)
         case ('-inp')
            input_file = trim(arg)
         case ('-out')
            output_file = trim(arg)
         case ('-col')
            read (arg,*) color_table
         case('-sig')
            read (arg,*) signal
         case ('-xsz')
            read (arg,*) xsize
         case ('-log')
            read (arg,*) logflg
         case ('-ash')
            read (arg,*) ashflag
         case ('-bar')
            read (arg,*) bar
         case ('-ttl')
            title = trim(arg)
         case ('-max')
            read (arg,*) usermax
         case ('-min')
            read (arg,*) usermin
         case ('-add')
            read (arg,*) offset
         case ('-mul')
            read (arg,*) factor
         case ('-pro')
            projection = trim(arg)
         case ('-lat')
            read (arg,*) lat0
         case ('-lon')
            read (arg,*) lon0
         case ('-res')
            read (arg,*) gridresn
         case default
            print '("unknown option ",a4," ignored")', opt
         end select
      end do


      call assert_present(input_file)
      if (output_file(1:1) /= '!') then
         call assert_not_present(output_file)
      else
         output_file = output_file(2:len(output_file))
      endif

      if (usermax == usermin) then
         print *, 'it is silly to set the maximum and the minimum'
         print *, 'to the same value, is it not?'
         call fatal_error
      end if

      if (ashflag .and. logflg) then
         print *, '-log and -ash cannot both be true!'
         call fatal_error
      end if

1     format(/,&
           & " -inp [] ",A,/,&
           & " -out [] ",A/,&
           & " -col [] ",I1,/,&
           & " -max [] ",e10.2,/,&
           & " -min [] ",e10.2,/,&
           & " -add [] ",e10.2,/,&
           & " -mul [] ",e10.2,/,&
           & " -pro [] ",A,/,&
           & " -lon [] ",e10.2,/,&
           & " -lat [] ",e10.2,/,&
           & " -res [] ",e10.2,/,&
           & " -sig [] ",I1,/,&
           & " -xsz [] ",I6,/,&
           & " -log [] ",L1,/,&
           & " -ash [] ",L1,/,&
           & " -bar [] ",L1,/,&
           & " -ttl [] ",A,/)

   end subroutine parseOptions

END PROGRAM MAP2GIF
