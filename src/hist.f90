
module hist

   use decomp_2d
   use iso_fortran_env, only: output_unit

   implicit none

   integer, parameter :: nsamp = 100

   ! Object with the histogram
   type :: hist_type

      ! Number of events in each bin
      integer, dimension(:), allocatable :: num
      ! Min and max values for the bins
      real(mytype) :: hist_min, hist_max, step
      ! Average and max value out of bounds
      real(mytype) :: hist_minavg, hist_maxavg
      real(mytype) :: hist_minval, hist_maxval

   contains

      procedure, public :: init => hist_type_init
      procedure, public :: fin => hist_type_fin
      procedure, public :: print => hist_type_print
      generic, public :: update => hist_type_update_scalar!, update_vector, ...
      procedure, private :: hist_type_update_scalar!, update_vector, ...

   end type hist_type

   private
   public :: hist_type

contains

   !
   ! Initialize the object with given min/max bounds
   !
   subroutine hist_type_init(obj, obj_min, obj_max)

      implicit none

      class(hist_type) :: obj
      real(mytype), intent(in) :: obj_min, obj_max

      ! Safety check
      if (abs(obj_min - obj_max) < nsamp*tiny(obj_min)) return

      ! Save bounds
      obj%hist_min = min(obj_min, obj_max)
      obj%hist_max = max(obj_min, obj_max)
      obj%step = (obj%hist_max - obj%hist_min)/nsamp

      ! Safety check
      if (allocated(obj%num)) deallocate (obj%num)

      ! Allocate memory
      allocate (obj%num(nsamp + 2))

      ! Init values
      obj%num = 0
      obj%hist_minavg = 0._mytype
      obj%hist_maxavg = 0._mytype
      obj%hist_minval = HUGE(obj%hist_minval)
      obj%hist_maxval = -HUGE(obj%hist_minval)

   end subroutine hist_type_init

   !
   ! Release the memory
   !
   subroutine hist_type_fin(obj)

      implicit none

      class(hist_type) :: obj

      if (allocated(obj%num)) deallocate (obj%num)

   end subroutine hist_type_fin

   !
   ! Update the histogram with the given scalar
   !
   subroutine hist_type_update_scalar(obj, input)

      implicit none

      ! Arguments
      class(hist_type) :: obj
      real(mytype), intent(in) :: input

      ! Local variables
      integer :: ibin

      ! Safety check
      if (.not. allocated(obj%num)) return

      ! Find the size bin
      ibin = 1 + ceiling((input - obj%hist_min)/obj%step)

      ! Enforce bounds
      if (ibin < 1) then
         ibin = 1
      else if (ibin > size(obj%num)) then
         ibin = size(obj%num)
      end if

      ! Update num
      obj%num(ibin) = obj%num(ibin) + 1

      ! Special case, out of bounds
      if (ibin == 1) then
         obj%hist_minavg = obj%hist_minavg + (input - obj%hist_minavg)/obj%num(ibin)
         obj%hist_minval = min(obj%hist_minval, input)
      else if (ibin == size(obj%num)) then
         obj%hist_maxavg = obj%hist_maxavg + (input - obj%hist_maxavg)/obj%num(ibin)
         obj%hist_maxval = max(obj%hist_maxval, input)
      end if

   end subroutine hist_type_update_scalar

   !
   ! Print the histogram to a given IO unit or stdout
   !
   subroutine hist_type_print(obj, given_io_unit)

      implicit none

      ! Arguments
      class(hist_type) :: obj
      integer, intent(in), optional :: given_io_unit

      ! Local variables
      integer :: ibin, io_unit

      ! Safety check
      if (.not. allocated(obj%num)) return

      ! Given IO unit or stdout
      if (present(given_io_unit)) then
         io_unit = given_io_unit
      else
         io_unit = output_unit
      end if

      ! Output each size bin
      !
      ! Text output if stdout
      !
      ! Binary output otherwise
      !
      if (present(io_unit)) then
         write (io_unit) size(obj%num)
         if (obj%num(1) == 0) then
            write (io_unit) obj%hist_min, obj%hist_min, obj%num(1), obj%hist_minavg
         else
            write (io_unit) obj%hist_minval, obj%hist_min, obj%num(1), obj%hist_minavg
         end if
         do ibin = 2, size(obj%num) - 1
            write (io_unit) obj%hist_min + (ibin - 2)*obj%step, obj%hist_min + (ibin - 1)*obj%step, obj%num(ibin)
         end do
         if (obj%num(size(obj%num)) == 0) then
            write (io_unit, *) obj%hist_max, obj%hist_max, obj%num(size(obj%num)), obj%hist_maxavg
         else
            write (io_unit, *) obj%hist_max, obj%hist_maxval, obj%num(size(obj%num)), obj%hist_maxavg
         end if
      else
         if (obj%num(1) == 0) then
            write (io_unit, *) "No sample below min"
         else
            write (io_unit, *) "Bin ]", obj%hist_minval, ", ", obj%hist_min, "] ", &
               obj%num(1), obj%hist_minavg
         end if
         do ibin = 2, size(obj%num) - 1
            write (io_unit, *) "Bin ]", obj%hist_min + (ibin - 2)*obj%step, ", ", &
               obj%hist_min + (ibin - 1)*obj%step, "] ", obj%num(ibin)
         end do
         if (obj%num(size(obj%num)) == 0) then
            write (io_unit, *) "No sample above max"
         else
            write (io_unit, *) "Bin ]", obj%hist_max, ", ", obj%hist_maxval, "[ ", &
               obj%num(size(obj%num)), obj%hist_maxavg
         end if
      end if

   end subroutine hist_type_print

end module hist
