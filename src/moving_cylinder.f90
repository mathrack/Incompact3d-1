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

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank, decomp_2d_abort
  use ibm_param, only : cex, cey, ubcx, ubcy, cyl_oscil, cyl_period, cyl_amp

  implicit none

  private ! All functions/subroutines private by default
  public :: get_cyl_xpos, get_cyl_ypos, get_cyl_xvel, get_cyl_yvel

contains

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
                                                                                           
    get_cyl_ypos = cey
    if (time > 0._mytype) then
      if (cyl_oscil) then
        get_cyl_ypos = get_cyl_ypos + cyl_amp * sin(cyl_period * time)                     
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

    if (cyl_oscil) then
      get_cyl_yvel = cyl_amp * cyl_period * cos(cyl_period * time)
    else
      get_cyl_yvel = ubcy
    end if

  end function get_cyl_yvel

end module moving_cylinder
