!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to Z pencil

subroutine transpose_y_to_z_real_s(src, dst)

  implicit none

  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst

  call transpose_y_to_z_real(src, dst, decomp_main)

end subroutine transpose_y_to_z_real_s

subroutine transpose_y_to_z_real(src, dst, decomp)

  implicit none

  ! Arguments
  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_Neighbor_alltoallw( &
      src,decomp%xcnts_yz(decomp%xranks_yz),decomp%xdispls_yz(decomp%xranks_yz),decomp%xtypes_yzr(decomp%xranks_yz), &
      dst,decomp%zcnts_yz(decomp%zranks_yz),decomp%zdispls_yz(decomp%zranks_yz),decomp%ztypes_yzr(decomp%zranks_yz), &
      decomp%ytozNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
  call MPI_Alltoallw(src,decomp%xcnts_yz,decomp%xdispls_yz,decomp%xtypes_yzr, &
                     dst,decomp%zcnts_yz,decomp%zdispls_yz,decomp%ztypes_yzr,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

end subroutine transpose_y_to_z_real

subroutine transpose_y_to_z_complex_s(src, dst)

  implicit none

  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst

  call transpose_y_to_z_complex(src, dst, decomp_main)

end subroutine transpose_y_to_z_complex_s

subroutine transpose_y_to_z_complex(src, dst, decomp)

  implicit none

  ! Arguments
  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_Neighbor_alltoallw( &
      src,decomp%xcnts_yz(decomp%xranks_yz),decomp%xdispls_yz(decomp%xranks_yz),decomp%xtypes_yzc(decomp%xranks_yz), &
      dst,decomp%zcnts_yz(decomp%zranks_yz),decomp%zdispls_yz(decomp%zranks_yz),decomp%ztypes_yzc(decomp%zranks_yz), &
      decomp%ytozNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
  call MPI_Alltoallw(src,decomp%xcnts_yz,decomp%xdispls_yz,decomp%xtypes_yzc, &
                     dst,decomp%zcnts_yz,decomp%zdispls_yz,decomp%ztypes_yzc,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

end subroutine transpose_y_to_z_complex

#ifdef OCC

subroutine transpose_y_to_z_real_start_s(src, dst, req)

  implicit none

  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
  integer, intent(out) :: req

  call transpose_y_to_z_real_start(src, dst, req, decomp_main)

end subroutine transpose_y_to_z_real_start_s

subroutine transpose_y_to_z_real_start(src, dst, req, decomp)

  implicit none

  ! Arguments
  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
  integer, intent(out) :: req
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_INeighbor_alltoallw( &
      src,decomp%xcnts_yz(decomp%xranks_yz),decomp%xdispls_yz(decomp%xranks_yz),decomp%xtypes_yzr(decomp%xranks_yz), &
      dst,decomp%zcnts_yz(decomp%zranks_yz),decomp%zdispls_yz(decomp%zranks_yz),decomp%ztypes_yzr(decomp%zranks_yz), &
      decomp%ytozNeighborComm,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
  call MPI_IAlltoallw(src,decomp%xcnts_yz,decomp%xdispls_yz,decomp%xtypes_yzr, &
                      dst,decomp%zcnts_yz,decomp%zdispls_yz,decomp%ztypes_yzr,MPI_COMM_WORLD,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_IALLTOALLW")
#endif

end subroutine transpose_y_to_z_real_start

subroutine transpose_y_to_z_complex_start_s(src, dst, req)

  implicit none

  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
  integer, intent(out) :: req

  call transpose_y_to_z_complex_start(src, dst, req, decomp_main)

end subroutine transpose_y_to_z_complex_start_s

subroutine transpose_y_to_z_complex_start(src, dst, req, decomp)

  implicit none

  ! Arguments
  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
  integer, intent(out) :: req
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_INeighbor_alltoallw( &
      src,decomp%xcnts_yz(decomp%xranks_yz),decomp%xdispls_yz(decomp%xranks_yz),decomp%xtypes_yzc(decomp%xranks_yz), &
      dst,decomp%zcnts_yz(decomp%zranks_yz),decomp%zdispls_yz(decomp%zranks_yz),decomp%ztypes_yzc(decomp%zranks_yz), &
      decomp%ytozNeighborComm,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_INEIGHBOR_ALLTOALLW")
#else
  call MPI_IAlltoallw(src,decomp%xcnts_yz,decomp%xdispls_yz,decomp%xtypes_yzc, &
                     dst,decomp%zcnts_yz,decomp%zdispls_yz,decomp%ztypes_yzc,MPI_COMM_WORLD,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_IALLTOALLW")
#endif

end subroutine transpose_y_to_z_complex_start

#endif
