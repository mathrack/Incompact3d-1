!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: moving_cylinder.f90
!!! DESCRIPTION: This module provides the position / velocity of a moving cylinder
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module moving_cylinder

  use MPI
  use, intrinsic :: iso_fortran_env, only: real32
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi, only : nrank, decomp_2d_abort
  use param, only : twopi
  use ibm_param, only : cex, cey, ubcx, ubcy, cyl_oscil, cyl_period, cyl_amp

  implicit none

  logical, save :: pos_vel_from_expe
  logical, parameter :: debug = .true.
  real(mytype), dimension(:,:), allocatable, save :: cyl_data

  private ! All functions/subroutines private by default
  public :: get_cyl_xpos, get_cyl_ypos, get_cyl_xvel, get_cyl_yvel, &
            moving_cylinder_init, moving_cylinder_fin

contains

  !
  ! This is called at the beginning of the simulation to initialize the module
  !
  ! If present, the argument is the name of the text file containing time, position and velocity
  !
  subroutine moving_cylinder_init(filename)

    implicit none

    ! Argument
    character(len=*), intent(in), optional :: filename

    ! Local variables
    integer :: iounit, i, n, code

    ! Define the flag for experimental data
    if (present(filename)) then
      pos_vel_from_expe = .true.
    else
      pos_vel_from_expe = .false.
    end if
    if (debug) write(*,*) "DEBUG. Moving cylinder, CPU ", nrank, " expe ", pos_vel_from_expe

    ! Allocate memory and read data if needed
    if (pos_vel_from_expe) then

       ! Open the file, read the number of lines, broadcast
       if (nrank.eq.0) then
          open(newunit=iounit, &
               file=trim(filename), &
               action='read')
          read(iounit,*) n
       end if
       call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
       if (debug) write(*,*) "DEBUG. Moving cylinder, CPU ", nrank, " number of lines ", n ! This can be removed later

       ! Safety check
       if (n<=0) then
          pos_vel_from_expe = .false.
          return
       end if

       ! Allocate memory
       allocate(cyl_data(3,n))
       cyl_data = 0._mytype

       ! Read the file, broadcast
       if (nrank.eq.0) then
          do i = 1, n
             read(iounit,*) cyl_data(:,i)
          end do
          close(iounit)
       end if
       call MPI_BCAST(cyl_data, 3*n, real_type, 0, MPI_COMM_WORLD,code)
       if (debug) write(*,*) "DEBUG. Moving cylinder, CPU ", nrank, " data ", real(cyl_data(:,1:2), kind=real32) ! This can be removed later

    end if

  end subroutine moving_cylinder_init

  !
  ! This is called at the end of the simulation to finalize the module
  !
  subroutine moving_cylinder_fin()

    implicit none

    if (allocated(cyl_data)) deallocate(cyl_data)

  end subroutine moving_cylinder_fin

  !
  ! Return the location of the moving cylinder in the direction x
  !
  function get_cyl_xpos(time)

    implicit none

    ! Argument and output value
    real(mytype), intent(in) :: time
    real(mytype) :: get_cyl_xpos

    get_cyl_xpos = cex
    if (time > 0._mytype) then
      if (cyl_oscil) then
        ! Nothing to do
        return
      else
        get_cyl_xpos = get_cyl_xpos + ubcx * time
      end if
    end if

  end function get_cyl_xpos

  ! 
  ! Return the location of the moving cylinder in the direction y
  !                                                                                        
  function get_cyl_ypos(time)
    
    implicit none                          
    
    ! Argument and output value
    real(mytype), intent(in) :: time
    real(mytype) :: get_cyl_ypos                                                           

    if (pos_vel_from_expe) then
       get_cyl_ypos = cey + interpolate(cyl_data, 1, time, 2)
       if (debug) write(*,*) "DEBUG. Moving cylinder, CPU ", nrank, " position ", real(get_cyl_ypos, kind=real32) ! This can be removed later
       return
    end if

    get_cyl_ypos = cey
    if (time >= 0._mytype) then
      if (cyl_oscil) then
        get_cyl_ypos = get_cyl_ypos + cyl_amp * sin(twopi * (time / cyl_period + 1._mytype/4._mytype))
      else
        get_cyl_ypos = get_cyl_ypos + ubcy * time                                          
      end if
    end if
       
  end function get_cyl_ypos

  !
  ! Return the velocity of the moving cylinder in the direction x
  !
  function get_cyl_xvel(time)

    implicit none

    ! Argument and output value
    real(mytype), intent(in) :: time
    real(mytype) :: get_cyl_xvel

    if (cyl_oscil) then
      get_cyl_xvel = 0._mytype
    else
      get_cyl_xvel = ubcx
    end if

  end function get_cyl_xvel

  !
  ! Return the velocity of the moving cylinder in the direction y
  !
  function get_cyl_yvel(time)

    implicit none

    ! Argument and output value
    real(mytype), intent(in) :: time
    real(mytype) :: get_cyl_yvel

    if (pos_vel_from_expe) then
       get_cyl_yvel = interpolate(cyl_data, 1, time, 3)
       if (debug) write(*,*) "DEBUG. Moving cylinder, CPU ", nrank, " velocity ", real(get_cyl_yvel, kind=real32) ! This can be removed later
       return
    end if

    if (cyl_oscil) then
      get_cyl_yvel = cyl_amp * (twopi / cyl_period) * cos(twopi * (time / cyl_period + 1._mytype/4._mytype))
    else
      get_cyl_yvel = ubcy
    end if

  end function get_cyl_yvel

  !
  ! Interpolate the input data
  !
  function interpolate(data, ix, x, iy)

    implicit none

    ! Arguments and output value
    real(mytype), intent(in) :: data(:,:), x
    integer, intent(in) :: ix, iy
    real(mytype) :: interpolate

    ! Local value
    integer :: i1

    if (x <= data(ix,1)) then
      interpolate = data(iy,1)
      return
    end if

    if (x >= data(ix,size(data,2))) then
      interpolate = data(iy,size(data,2))
      return
    end if

    ! Assume constant time step
    i1 = 1 + int((x-data(ix,1))/(data(ix,2)-data(ix,1)))

    ! Safety
    i1 = min(1, i1)
    i1 = max(i1, size(data,2)-1)

    ! Linear interpolation
    interpolate = data(iy,i1) + (data(iy,i1+1)-data(iy,i1)) * (x-data(ix,i1)) / (data(ix,i1+1)-data(ix,i1))

  end function interpolate

end module moving_cylinder
