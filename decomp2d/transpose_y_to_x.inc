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

! This file contains the routines that transpose data from Y to X pencil

subroutine transpose_y_to_x_real_s(src, dst)

implicit none

real(mytype), dimension(:,:,:), intent(IN) :: src
real(mytype), dimension(:,:,:), intent(OUT) :: dst

call transpose_y_to_x_real(src, dst, decomp_main)

end subroutine transpose_y_to_x_real_s

subroutine transpose_y_to_x_real(src, dst, decomp)

implicit none

! Arguments
real(mytype), dimension(:,:,:), intent(IN) :: src
real(mytype), dimension(:,:,:), intent(OUT) :: dst
TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
    call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xy(decomp%zranks_xy),decomp%zdispls_xy(decomp%zranks_xy),decomp%ztypes_xyr(decomp%zranks_xy), &
      dst,decomp%xcnts_xy(decomp%xranks_xy),decomp%xdispls_xy(decomp%xranks_xy),decomp%xtypes_xyr(decomp%xranks_xy), &
      decomp%ytoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
    call MPI_Alltoallw(src,decomp%zcnts_xy,decomp%zdispls_xy,decomp%ztypes_xyr, &
                       dst,decomp%xcnts_xy,decomp%xdispls_xy,decomp%xtypes_xyr,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

return
end subroutine transpose_y_to_x_real

subroutine transpose_y_to_x_complex_s(src, dst)

implicit none

complex(mytype), dimension(:,:,:), intent(IN) :: src
complex(mytype), dimension(:,:,:), intent(OUT) :: dst

call transpose_y_to_x_complex(src, dst, decomp_main)

end subroutine transpose_y_to_x_complex_s

subroutine transpose_y_to_x_complex(src, dst, decomp)

implicit none

! Arguments
complex(mytype), dimension(:,:,:), intent(IN) :: src
complex(mytype), dimension(:,:,:), intent(OUT) :: dst
TYPE(DECOMP_INFO), intent(IN) :: decomp

  ! Local variable
  integer :: ierror

#ifdef MPI3
    call MPI_Neighbor_alltoallw( &
      src,decomp%zcnts_xy(decomp%zranks_xy),decomp%zdispls_xy(decomp%zranks_xy),decomp%ztypes_xyc(decomp%zranks_xy), &
      dst,decomp%xcnts_xy(decomp%xranks_xy),decomp%xdispls_xy(decomp%xranks_xy),decomp%xtypes_xyc(decomp%xranks_xy), &
      decomp%ytoxNeighborComm,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_NEIGHBOR_ALLTOALLW")
#else
    call MPI_Alltoallw(src,decomp%zcnts_xy,decomp%zdispls_xy,decomp%ztypes_xyc, &
                       dst,decomp%xcnts_xy,decomp%xdispls_xy,decomp%xtypes_xyc,MPI_COMM_WORLD,ierror)
  if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLTOALLW")
#endif

return
end subroutine transpose_y_to_x_complex
