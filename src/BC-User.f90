!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module user_sim

   USE decomp_2d
   use decomp_2d_io
   USE variables
   USE param
   use ibm_param
   use MPI
   use var, only: di1, di2, ux2, uy2, phi2, ta1, ta2, nzmsize
   use iso_fortran_env, only: output_unit
   use hist, only: hist_type

   IMPLICIT NONE

   ! Flags to control monitoring
   logical, save :: init_not_done = .true.
   logical, parameter :: monitor_bulk = .true.
   logical, parameter :: monitor_minmax = .false.
   logical, parameter :: use_hist = .true.

   ! Temperature on the left wall
   real(mytype), parameter :: temp_l = 0.5_mytype

   ! Temperature difference between both walls
   real(mytype), parameter :: deltaT = 1._mytype

   ! IO unit and file name for bulk quantities
   integer, save :: io_bulk = output_unit
   character(len=*), parameter :: bulk_file = "cavity_bulk.dat"

   ! Monitor minmax
   real(mytype), save, dimension(:, :, :, :), allocatable :: mnmx_u, mnmx_v, mnmx_t, &
                                                             mnmx_udx, mnmx_vdx, mnmx_tdx, &
                                                             mnmx_udy, mnmx_vdy, mnmx_tdy

   ! Histograms
   type(hist_type), save, dimension(:, :), allocatable :: hst_u, hst_v, hst_t, &
                                                          hst_udx, hst_vdx, hst_tdx, &
                                                          hst_udy, hst_vdy, hst_tdy

   PRIVATE ! All functions/subroutines private by default
   PUBLIC :: init_user, boundary_conditions_user, postprocess_user, &
             visu_user_init, visu_user, fin_user

contains

   subroutine init_user(ux1, uy1, uz1, ep1, phi1)

      USE decomp_2d_io
      USE MPI

      implicit none

      ! Arguments
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      ! Local variables
      integer :: i, j

      ! Do not initialize twice
      init_not_done = .false.

      ! Open IO unit for bulk quantities on master rank
      if (monitor_bulk .and. nrank == 0) then
         open (newunit=io_bulk, file=bulk_file, form='formatted')
         write (io_bulk, *) "u      v       T       u'²     v'²     T'²"
      end if

      ! This does not apply in case of restart
      if (irestart==0) then

         ! Velocity is zero
         ux1 = zero
         uy1 = zero
         uz1 = zero

         ! Linear temperature profile
         if (numscalar >= 1) then
            do i = 1, xsize(1)
               phi1(i, :, :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if

      endif

      ! Monitor minmax
      if (monitor_minmax) then
         allocate (mnmx_u(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_v(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_t(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_udx(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_vdx(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_tdx(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_udy(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_vdy(xsize(1), xsize(2), xsize(3), 2))
         allocate (mnmx_tdy(xsize(1), xsize(2), xsize(3), 2))
      end if

      ! Add histogram everywhere for u, v, T and the space derivative
      if (use_hist) then
         allocate (hst_u(xsize(1), xsize(2)))
         allocate (hst_v(xsize(1), xsize(2)))
         allocate (hst_t(xsize(1), xsize(2)))
         allocate (hst_udx(xsize(1), xsize(2)))
         allocate (hst_vdx(xsize(1), xsize(2)))
         allocate (hst_tdx(xsize(1), xsize(2)))
         allocate (hst_udy(xsize(1), xsize(2)))
         allocate (hst_vdy(xsize(1), xsize(2)))
         allocate (hst_tdy(xsize(1), xsize(2)))
         call decomp_2d_read_one(1,ta1,".","min_u","min_u",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_u","max_u",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_u(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_v","min_v",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_v","max_v",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_v(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_t","min_t",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_t","max_t",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_t(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_udx","min_udx",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_udx","max_udx",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_udx(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_vdx","min_vdx",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_vdx","max_vdx",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_vdx(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_tdx","min_tdx",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_tdx","max_tdx",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_tdx(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_udy","min_udy",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_udy","max_udy",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_udy(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_vdy","min_vdy",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_vdy","max_vdy",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_vdy(i, j)%init(ta1(i,j,1), di1(i,j,1))
            enddo
         enddo
         call decomp_2d_read_one(1,ta1,".","min_tdy","min_tdy",decomp_main)
         call decomp_2d_read_one(1,di1,".","max_tdy","max_tdy",decomp_main)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_tdy(i, j)%init(ta1(i,j,1), di1(i,j,1))
            end do
         end do
      end if

   end subroutine init_user

   subroutine boundary_conditions_user(ux, uy, uz, phi, ep)

      implicit none

      ! Arguments
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, ep
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi

      ! Local arguments
      integer :: i

      ! Velocity
      IF (nclx1 == 2) THEN
      END IF
      IF (nclxn == 2) THEN
      END IF
      IF (ncly1 == 2 .and. xstart(2) == 1) THEN
      END IF
      IF (nclyn == 2 .and. xend(2) == ny) THEN
      END IF

      ! Scalar
      if (numscalar >= 1) then
         if (nclxS1 == 2) then
            phi(1, :, :, :) = temp_l
         end if
         if (nclxSn == 2) then
            phi(xsize(1), :, :, :) = temp_l - deltaT
         end if
         if (nclyS1 == 2 .and. xstart(2) == 1) then
            do i = 1, xsize(1)
               phi(i, 1, :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
         if (nclySn == 2 .and. xend(2) == ny) then
            do i = 1, xsize(1)
               phi(i, xsize(2), :, :) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
      end if

   end subroutine boundary_conditions_user

   subroutine postprocess_user(ux1, uy1, uz1, phi1, ep1)

      implicit none

      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      ! Local variables
      real(mytype), dimension(3) :: phiavg, phivar
      integer :: ivar, i, j, k, ierror

      ! Nothing to do if no scalar
      if (numscalar <= 0) return

      ! If first time step and restart
      if (init_not_done) call init_user(ux1, uy1, uz1, ep1, phi1)

      ! Monitor bulk quantities
      if (monitor_bulk) then

         ! Init
         phiavg = zero
         phivar = zero

         ! Monitor the space-averaged and RMS of velocity and temperature
         ! First, compute the local sum
         ivar = 1
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               do i = 1, xsize(1)
                  phiavg(ivar) = phiavg(ivar) + ux1(i, j, k)
                  phivar(ivar) = phivar(ivar) + ux1(i, j, k)**2
               end do
            end do
         end do
         ivar = ivar + 1
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               do i = 1, xsize(1)
                  phiavg(ivar) = phiavg(ivar) + uy1(i, j, k)
                  phivar(ivar) = phivar(ivar) + uy1(i, j, k)**2
               end do
            end do
         end do
         ivar = ivar + 1
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               do i = 1, xsize(1)
                  phiavg(ivar) = phiavg(ivar) + phi1(i, j, k, 1)
                  phivar(ivar) = phivar(ivar) + phi1(i, j, k, 1)**2
               end do
            end do
         end do
         ! Parallel sum if needed
         if (nproc > 1) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, &
                               (/phiavg, phivar/), &
                               6, &
                               real_type, &
                               MPI_SUM, &
                               MPI_COMM_WORLD, &
                               ierror)
            if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_REDUCE in postprocess_user")
         end if
         ! Rescale
         phiavg = phiavg/real(nx*ny, kind=mytype)
         phivar = phivar/real(nx*ny, kind=mytype)
         ! Rank 0 save the values in the file
         ! Replace E14.6 => E24.16 to print all the digits
         if (io_bulk /= output_unit) write (io_bulk, '(6(E14.6))') phiavg, phivar - phiavg**2

      end if

      ! Monitor min/max
      if (monitor_minmax) then

         ! u, v, T
         call update_minmax(mnmx_u, ux1)
         call update_minmax(mnmx_v, uy1)
         call update_minmax(mnmx_t, phi1(:,:,:,1))

         ! x-derivative
         call derx(ta1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call update_minmax(mnmx_udx, ta1)
         call derx(ta1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call update_minmax(mnmx_vdx, ta1)
         call derxS(ta1, phi1(:, :, :, 1), di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
         call update_minmax(mnmx_tdx, ta1)

         ! y-derivative
         call transpose_x_to_y(ux1, ux2, decomp_main)
         call dery(ta2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call update_minmax(mnmx_udy, ta1)
         call transpose_x_to_y(uy1, uy2, decomp_main)
         call dery(ta2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call update_minmax(mnmx_vdy, ta1)
         call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1), decomp_main)
         call deryS(ta2, phi2(:, :, :, 1), di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call update_minmax(mnmx_tdy, ta1)

      end if

      ! Histogram
      if (use_hist) then

         ! u, v, T
         call hst_update(hst_u, ux1)
         call hst_update(hst_v, uy1)
         call hst_update(hst_t, phi1(:,:,:,1))

         ! x-derivative
         call derx(ta1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0, ubcx)
         call hst_update(hst_udx, ta1)
         call derx(ta1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1, ubcy)
         call hst_update(hst_vdx, ta1)
         call derxS(ta1, phi1(:, :, :, 1), di1, sx, ffxpS, fsxpS, fwxpS, xsize(1), xsize(2), xsize(3), 1, zero)
         call hst_update(hst_tdx, ta1)

         ! y-derivative
         call transpose_x_to_y(ux1, ux2, decomp_main)
         call dery(ta2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1, ubcx)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call hst_update(hst_udy, ta1)
         call transpose_x_to_y(uy1, uy2, decomp_main)
         call dery(ta2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0, ubcy)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call hst_update(hst_vdy, ta1)
         call transpose_x_to_y(phi1(:, :, :, 1), phi2(:, :, :, 1), decomp_main)
         call deryS(ta2, phi2(:, :, :, 1), di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
         call transpose_y_to_x(ta2, ta1, decomp_main)
         call hst_update(hst_tdy, ta1)

      endif

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
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
      character(len=32), intent(in) :: num

   end subroutine visu_user

   ! Register fields
   subroutine visu_user_init(visu_initialised)

      implicit none

      logical, intent(out) :: visu_initialised

      visu_initialised = .true.

   end subroutine visu_user_init

   ! Print stuff and release memory
   subroutine fin_user()

      implicit none

      integer :: i, j, io_unit
      character(len=30) :: filename

      if (monitor_bulk .and. nrank == 0) close (io_bulk)

      if (monitor_minmax) then
         ! Write planes
         call decomp_2d_write_one(1, mnmx_u(:, :, :, 1), "min_u", decomp_main)
         call decomp_2d_write_one(1, mnmx_u(:, :, :, 2), "max_u", decomp_main)
         call decomp_2d_write_one(1, mnmx_v(:, :, :, 1), "min_v", decomp_main)
         call decomp_2d_write_one(1, mnmx_v(:, :, :, 2), "max_v", decomp_main)
         call decomp_2d_write_one(1, mnmx_t(:, :, :, 1), "min_t", decomp_main)
         call decomp_2d_write_one(1, mnmx_t(:, :, :, 2), "max_t", decomp_main)
         call decomp_2d_write_one(1, mnmx_udx(:, :, :, 1), "min_udx", decomp_main)
         call decomp_2d_write_one(1, mnmx_udx(:, :, :, 2), "max_udx", decomp_main)
         call decomp_2d_write_one(1, mnmx_vdx(:, :, :, 1), "min_vdx", decomp_main)
         call decomp_2d_write_one(1, mnmx_vdx(:, :, :, 2), "max_vdx", decomp_main)
         call decomp_2d_write_one(1, mnmx_tdx(:, :, :, 1), "min_tdx", decomp_main)
         call decomp_2d_write_one(1, mnmx_tdx(:, :, :, 2), "max_tdx", decomp_main)
         call decomp_2d_write_one(1, mnmx_udy(:, :, :, 1), "min_udy", decomp_main)
         call decomp_2d_write_one(1, mnmx_udy(:, :, :, 2), "max_udy", decomp_main)
         call decomp_2d_write_one(1, mnmx_vdy(:, :, :, 1), "min_vdy", decomp_main)
         call decomp_2d_write_one(1, mnmx_vdy(:, :, :, 2), "max_vdy", decomp_main)
         call decomp_2d_write_one(1, mnmx_tdy(:, :, :, 1), "min_tdy", decomp_main)
         call decomp_2d_write_one(1, mnmx_tdy(:, :, :, 2), "max_tdy", decomp_main)

         ! Free memory
         deallocate (mnmx_u)
         deallocate (mnmx_v)
         deallocate (mnmx_t)
         deallocate (mnmx_udx)
         deallocate (mnmx_vdx)
         deallocate (mnmx_tdx)
         deallocate (mnmx_udy)
         deallocate (mnmx_vdy)
         deallocate (mnmx_tdy)

      end if

      if (use_hist) then
         ! Each rank print stuff to a dedicated IO unit
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               write(filename, "('out/histogram_',I3.3,'_',I3.3,'.txt')") i, j+xstart(2)-1
               open (newunit=io_unit, file=trim(filename))
               write (io_unit, *) ""
               write (io_unit, *) "u"
               call hst_u(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "v"
               call hst_v(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "t"
               call hst_t(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "udx"
               call hst_udx(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "vdx"
               call hst_vdx(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "tdx"
               call hst_tdx(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "udy"
               call hst_udy(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "vdy"
               call hst_vdy(i, j)%print(io_unit)
               write (io_unit, *) ""
               write (io_unit, *) "tdy"
               call hst_tdy(i, j)%print(io_unit)
               close (io_unit)
            end do
         end do

         ! Release memory
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               call hst_u(i, j)%fin()
               call hst_v(i, j)%fin()
               call hst_t(i, j)%fin()
               call hst_udx(i, j)%fin()
               call hst_vdx(i, j)%fin()
               call hst_tdx(i, j)%fin()
               call hst_udy(i, j)%fin()
               call hst_vdy(i, j)%fin()
               call hst_tdy(i, j)%fin()
            end do
         end do
         deallocate (hst_u)
         deallocate (hst_v)
         deallocate (hst_t)
         deallocate (hst_udx)
         deallocate (hst_vdx)
         deallocate (hst_tdx)
         deallocate (hst_udy)
         deallocate (hst_vdy)
         deallocate (hst_tdy)
      end if

   end subroutine fin_user

   !
   ! Update map of min-max
   !
   subroutine update_minmax(minmax, var)

      implicit none

      ! Arguments
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), 2), intent(inout) :: minmax
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: var

      ! Local variables
      integer :: i, j, k

      if (itime == ifirst) then
         do k = 1, xsize(3)
         do j = 1, xsize(2)
         do i = 1, xsize(1)
            minmax(i, j, k, 1) = var(i, j, k)
            minmax(i, j, k, 2) = var(i, j, k)
         end do
         end do
         end do
         return
      end if

      do k = 1, xsize(3)
      do j = 1, xsize(2)
      do i = 1, xsize(1)
         minmax(i, j, k, 1) = min(minmax(i, j, k, 1), var(i, j, k))
         minmax(i, j, k, 2) = max(minmax(i, j, k, 2), var(i, j, k))
      end do
      end do
      end do

   end subroutine update_minmax

   subroutine hst_update(hst, var)

      implicit none

      ! Arguments
      type(hist_type), dimension(xsize(1), xsize(2)), intent(inout) :: hst
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: var

      ! Local variables
      integer :: i, j

      do j = 1, xsize(2)
         do i = 1, xsize(1)
            call hst(i,j)%update(var(i,j,1))
         enddo
      enddo

   end subroutine hst_update

end module user_sim
