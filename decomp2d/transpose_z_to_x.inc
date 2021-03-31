! -*- mode: f90 -*-
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
      src,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xzr(decomp%zranks), &
      dst,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xzr(decomp%xranks), &
      decomp%ztoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
    call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzr, &
                       dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzr,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

  return
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
      src,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xzc(decomp%zranks), &
      dst,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xzc(decomp%xranks), &
      decomp%ztoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
    call MPI_Alltoallw(src,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xzc, &
                       dst,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xzc,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

  return
end subroutine transpose_z_to_x_complex
