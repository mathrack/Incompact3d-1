!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module user_sim

  USE decomp_2d
  USE variables
  USE param
  use MPI
  use var, only : nzmsize
  use iso_fortran_env, only : output_unit

  IMPLICIT NONE

  ! Temperature on the left wall
  real(mytype), parameter :: temp_l = 0.5_mytype

  ! Temperature difference between both walls
  real(mytype), parameter :: deltaT = 1._mytype

  ! IO unit and file name for bulk quantities
  integer, save :: io_bulk = output_unit
  character(len=*), parameter :: bulk_file = "cavity_bulk.dat"

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_user, boundary_conditions_user, postprocess_user, &
            visu_user_init, visu_user

contains

  subroutine init_user (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d_io
    USE MPI

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    ! Local variables
    integer :: i

    ! Open IO unit for bulk quantities on master rank
    if (nrank.eq.0) then
       open(newunit=io_bulk, file=bulk_file, form='formatted')
    endif

    ! Velocity is zero
    ux1 = zero
    uy1 = zero
    uz1 = zero

    ! Linear temperature profile
    if (numscalar.ge.1) then
       do i = 1, xsize(1)
          phi1(i,:,:,:) = temp_l - deltaT * (i-1) / real(xsize(1)-1, kind=mytype)
       enddo
    endif

  end subroutine init_user

  subroutine boundary_conditions_user (ux,uy,uz,phi,ep)

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    ! Local arguments
    integer :: i

    ! Velocity
    IF (nclx1.EQ.2) THEN
    ENDIF
    IF (nclxn.EQ.2) THEN
    ENDIF
    IF (ncly1.EQ.2.and.xstart(2).eq.1) THEN
    ENDIF
    IF (nclyn.EQ.2.and.xend(2).eq.ny) THEN
    ENDIF

    ! Scalar
    if (numscalar.ge.1) then
       if (nclxS1.eq.2) then
          phi(1,:,:,:) = temp_l
       endif
       if (nclxSn.eq.2) then
          phi(xsize(1),:,:,:) = temp_l - deltaT
       endif
       if (nclyS1.eq.2.and.xstart(2).eq.1) then
          do i = 1, xsize(1)
             phi(i,1,:,:) = temp_l - deltaT * (i-1) / real(xsize(1)-1, kind=mytype)
          enddo
       endif
       if (nclySn.eq.2.and.xend(2).eq.ny) then
          do i = 1, xsize(1)
             phi(i,xsize(2),:,:) = temp_l - deltaT * (i-1) / real(xsize(1)-1, kind=mytype)
          enddo
       endif
    endif

  end subroutine boundary_conditions_user

  subroutine postprocess_user(ux1,uy1,uz1,phi1,ep1)

    implicit none

    ! Arguments
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    ! Local variables
    real(mytype), dimension(3) :: phiavg, phivar
    integer :: ivar, i, j, k, ierror

    ! Nothing to do if no scalar
    if (numscalar.le.0) return

    ! Init
    phiavg = zero
    phivar = zero

    ! Monitor the space-averaged and RMS of velocity and temperature
    ! First, compute the local sum
    ivar = 1
    do k=1,xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             phiavg(ivar) = phiavg(ivar) + ux1(i,j,k)
             phivar(ivar) = phivar(ivar) + ux1(i,j,k)**2
          enddo
       enddo
    enddo
    ivar = ivar + 1
    do k=1,xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             phiavg(ivar) = phiavg(ivar) + uy1(i,j,k)
             phivar(ivar) = phivar(ivar) + uy1(i,j,k)**2
          enddo
       enddo
    enddo
    ivar = ivar + 1
    do k=1,xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             phiavg(ivar) = phiavg(ivar) + phi1(i,j,k,1)
             phivar(ivar) = phivar(ivar) + phi1(i,j,k,1)**2
          enddo
       enddo
    enddo
    ! Parallel sum if needed
    if (nproc>1) then
       call MPI_ALLREDUCE(MPI_IN_PLACE, &
                          (/phiavg,phivar/), &
                          6, &
                          real_type, &
                          MPI_SUM, &
                          MPI_COMM_WORLD, &
                          ierror)
       if (ierror/=0) call decomp_2d_abort(ierror, "MPI_REDUCE in postprocess_user")
    endif
    ! Rescale
    phiavg = phiavg / real(nx*ny, kind=mytype)
    phivar = phivar / real(nx*ny, kind=mytype)
    ! Rank 0 save the values in the file
    ! Replace E14.6 => E24.16 to print all the digits
    if (io_bulk/=output_unit) write(io_bulk,'(6(E14.6))') phiavg, phivar-phiavg**2

  end subroutine postprocess_user

  !############################################################################
  !!
  !!  SUBROUTINE: visu_user
  !!      AUTHOR: CF
  !! DESCRIPTION: Performs case-specific visualization
  !!
  !############################################################################
  subroutine visu_user(ux1, uy1, uz1, pp3, phi1, ep1, num)

    implicit none

    ! Arguments
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

  end subroutine visu_user

  ! Register fields
  subroutine visu_user_init(visu_initialised)

    implicit none

    logical, intent(out) :: visu_initialised

    visu_initialised = .true.

  end subroutine visu_user_init

end module user_sim
