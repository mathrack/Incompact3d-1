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

! This file contains the routines that transpose data from Z to X pencil

subroutine transpose_z_to_x_real_s(src, dst)

  implicit none

  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst

  call transpose_z_to_x_real(src, dst, decomp_main)

end subroutine transpose_z_to_x_real_s

subroutine transpose_z_to_x_real(src, dst, decomp)

  implicit none

  ! Arguments
  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks_xz),decomp%zdispls_xz(decomp%zranks_xz),decomp%ztypes_xzr(decomp%zranks_xz), &
      dst,decomp%xcnts_xz(decomp%xranks_xz),decomp%xdispls_xz(decomp%xranks_xz),decomp%xtypes_xzr(decomp%xranks_xz), &
      decomp%ztoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
  call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzr, &
                     dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzr,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

end subroutine transpose_z_to_x_real

subroutine transpose_z_to_x_complex_s(src, dst)

  implicit none

  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst

  call transpose_z_to_x_complex(src, dst, decomp_main)

end subroutine transpose_z_to_x_complex_s

subroutine transpose_z_to_x_complex(src, dst, decomp)

  implicit none

  ! Arguments
  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks_xz),decomp%zdispls_xz(decomp%zranks_xz),decomp%ztypes_xzc(decomp%zranks_xz), &
      dst,decomp%xcnts_xz(decomp%xranks_xz),decomp%xdispls_xz(decomp%xranks_xz),decomp%xtypes_xzc(decomp%xranks_xz), &
      decomp%ztoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
  call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzc, &
                     dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzc,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

end subroutine transpose_z_to_x_complex

#ifdef OCC

subroutine transpose_z_to_x_real_start_s(src, dst, req)

  implicit none

  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
#ifdef MPI3
  type(mpi_request), intent(out) :: req
#else
  integer, intent(out) :: req
#endif

  call transpose_z_to_x_real_start(src, dst, req, decomp_main)

end subroutine transpose_z_to_x_real_start_s

subroutine transpose_z_to_x_real_start(src, dst, req, decomp)

  implicit none

  ! Arguments
  real(mytype), dimension(:,:,:), intent(IN) :: src
  real(mytype), dimension(:,:,:), intent(OUT) :: dst
#ifdef MPI3
  type(mpi_request), intent(out) :: req
#else
  integer, intent(out) :: req
#endif
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_INeighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks_xz),decomp%zdispls_xz(decomp%zranks_xz),decomp%ztypes_xzr(decomp%zranks_xz), &
      dst,decomp%xcnts_xz(decomp%xranks_xz),decomp%xdispls_xz(decomp%xranks_xz),decomp%xtypes_xzr(decomp%xranks_xz), &
      decomp%ztoxNeighborComm,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_INEIGHBOR_ALLTOALLW")
#else
  call MPI_IAlltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzr, &
                      dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzr,MPI_COMM_WORLD,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_IALLTOALLW")
#endif

end subroutine transpose_z_to_x_real_start

subroutine transpose_z_to_x_complex_start_s(src, dst, req)

  implicit none

  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
#ifdef MPI3
  type(mpi_request), intent(out) :: req
#else
  integer, intent(out) :: req
#endif

  call transpose_z_to_x_complex_start(src, dst, req, decomp_main)

end subroutine transpose_z_to_x_complex_start_s

subroutine transpose_z_to_x_complex_start(src, dst, req, decomp)

  implicit none

  ! Arguments
  complex(mytype), dimension(:,:,:), intent(IN) :: src
  complex(mytype), dimension(:,:,:), intent(OUT) :: dst
#ifdef MPI3
  type(mpi_request), intent(out) :: req
#else
  integer, intent(out) :: req
#endif
  TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
  call MPI_INeighbor_alltoallw( &
      src,decomp%zcnts_xz(decomp%zranks_xz),decomp%zdispls_xz(decomp%zranks_xz),decomp%ztypes_xzc(decomp%zranks_xz), &
      dst,decomp%xcnts_xz(decomp%xranks_xz),decomp%xdispls_xz(decomp%xranks_xz),decomp%xtypes_xzc(decomp%xranks_xz), &
      decomp%ztoxNeighborComm,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_INEIGHBOR_ALLTOALLW")
#else
  call MPI_IAlltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzc, &
                      dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzc,MPI_COMM_WORLD,req,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_IALLTOALLW")
#endif

end subroutine transpose_z_to_x_complex_start

#endif
