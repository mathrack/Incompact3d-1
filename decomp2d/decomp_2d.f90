!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the main 2D pencil decomposition module

module decomp_2d

  use MPI

  implicit none

  private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
  integer, parameter, public :: mytype = KIND(0.0D0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: real2_type = MPI_2DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef SAVE_SINGLE
  integer, parameter, public :: mytype_single = KIND(0.0)
  integer, parameter, public :: real_type_single = MPI_REAL
#else
  integer, parameter, public :: mytype_single = KIND(0.0D0)
  integer, parameter, public :: real_type_single = MPI_DOUBLE_PRECISION
#endif
#else
  integer, parameter, public :: mytype = KIND(0.0)
  integer, parameter, public :: real_type = MPI_REAL
  integer, parameter, public :: real2_type = MPI_2REAL
  integer, parameter, public :: complex_type = MPI_COMPLEX
  integer, parameter, public :: mytype_single = KIND(0.0)
  integer, parameter, public :: real_type_single = MPI_REAL
#endif

  integer, save, public :: mytype_bytes

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank  ! local MPI rank 
  integer, save, public :: nproc  ! total number of processors

  ! parameters for 2D Cartesian topology 
  integer, save, dimension(2) :: dims, coord
  logical, save, dimension(2) :: periodic
  integer, save, public :: DECOMP_2D_COMM_CART_X, &
       DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z 
  integer, save :: DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL

  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(3,6) :: neighbour 

  ! flags for periodic condition in three dimensions
  logical, save :: periodic_x, periodic_y, periodic_z

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist, &
          x1st, y1st, y2st, z2st, &
          x1en, y1en, y2en, z2en

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp

     ! buffer counts for MPI_ALLTOALL: either for evenly distributed data
     ! or for padded-alltoall
     integer :: x1count, y1count, y2count, z2count


     ! buffer counts, displacements and types for MPI_Alltoallw to transform
     ! directly between x- and z-pencils
     integer, allocatable, dimension(:) :: xcnts_xz, xtypes_xzr, xtypes_xzc
     integer, allocatable, dimension(:) :: zcnts_xz, ztypes_xzr, ztypes_xzc
     ! directly between x- and y-pencils
     integer, allocatable, dimension(:) :: xcnts_xy, xtypes_xyr, xtypes_xyc
     integer, allocatable, dimension(:) :: zcnts_xy, ztypes_xyr, ztypes_xyc
     ! directly between y- and z-pencils
     integer, allocatable, dimension(:) :: xcnts_yz, xtypes_yzr, xtypes_yzc
     integer, allocatable, dimension(:) :: zcnts_yz, ztypes_yzr, ztypes_yzc

#ifdef MPI3                               
     ! use MPI_ADDRESS_KIND for MPI_Neighbor_alltoallw call
     ! x <=> z transpose
     integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: xdispls_xz, zdispls_xz
     integer :: xtozNeighborComm, ztoxNeighborComm
     integer, allocatable, dimension(:) :: xranks_xz, zranks_xz
     ! x <=> y transpose
     integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: xdispls_xy, zdispls_xy
     integer :: xtoyNeighborComm, ytoxNeighborComm
     integer, allocatable, dimension(:) :: xranks_xy, zranks_xy
     ! y <=> z transpose
     integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: xdispls_yz, zdispls_yz
     integer :: ytozNeighborComm, ztoyNeighborComm
     integer, allocatable, dimension(:) :: xranks_yz, zranks_yz
#else                                     
     ! use default integer for MPI_Alltoallw call
     ! x <=> z transpose
     integer, allocatable, dimension(:) :: xdispls_xz, zdispls_xz
     ! x <=> y transpose
     integer, allocatable, dimension(:) :: xdispls_xy, zdispls_xy
     ! y <=> z transpose
     integer, allocatable, dimension(:) :: xdispls_yz, zdispls_yz
#endif

     ! evenly distributed data
     logical :: even

  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save :: decomp_main
  TYPE(DECOMP_INFO), save, public :: phG,ph1,ph2,ph3,ph4

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To define smaller arrays using every several mesh points
  integer, save, dimension(3), public :: xszS,yszS,zszS,xstS,ystS,zstS,xenS,yenS,zenS
  integer, save, dimension(3), public :: xszV,yszV,zszV,xstV,ystV,zstV,xenV,yenV,zenV
  integer, save, dimension(3), public :: xszP,yszP,zszP,xstP,ystP,zstP,xenP,yenP,zenP
  logical, save :: coarse_mesh_starts_from_1
  integer, save :: iskipS, jskipS, kskipS
  integer, save :: iskipV, jskipV, kskipV
  integer, save :: iskipP, jskipP, kskipP


  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_x_to_z, transpose_z_to_x, &
       transpose_z_to_y, transpose_y_to_x, &
#ifdef OCC
       transpose_x_to_y_start, transpose_y_to_z_start, &
       transpose_z_to_y_start, transpose_y_to_x_start, &
       transpose_x_to_y_wait, transpose_y_to_z_wait, &
       transpose_z_to_y_wait, transpose_y_to_x_wait, &
       transpose_test, &
#endif
       decomp_info_init, decomp_info_finalize, partition, &
       init_coarser_mesh_statS,fine_to_coarseS,&
       init_coarser_mesh_statV,fine_to_coarseV,&
       init_coarser_mesh_statP,fine_to_coarseP,&
       alloc_x, alloc_y, alloc_z, &
       update_halo, decomp_2d_abort, &
       get_decomp_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are routines to perform global data transpositions
  ! 
  !   Four combinations are available, enough to cover all situations
  !    - transpose_x_to_y (X-pencil --> Y-pencil)
  !    - transpose_y_to_z (Y-pencil --> Z-pencil)
  !    - transpose_z_to_y (Z-pencil --> Y-pencil)
  !    - transpose_y_to_x (Y-pencil --> X-pencil)
  !
  !   Generic interface provided here to support multiple data types
  !    - real and complex types supported through generic interface
  !    - single/double precision supported through pre-processing
  !       * see 'mytype' variable at the beginning
  !    - an optional argument can be supplied to transpose data whose 
  !      global size is not the default nx*ny*nz 
  !       * as the case in fft r2c/c2r interface 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface transpose_x_to_y
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_real_s
     module procedure transpose_x_to_y_complex
     module procedure transpose_x_to_y_complex_s
  end interface transpose_x_to_y

  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_real_s
     module procedure transpose_y_to_z_complex
     module procedure transpose_y_to_z_complex_s
  end interface transpose_y_to_z

  interface transpose_x_to_z
     module procedure transpose_x_to_z_real
     module procedure transpose_x_to_z_real_s
     module procedure transpose_x_to_z_complex
     module procedure transpose_x_to_z_complex_s
  end interface transpose_x_to_z

  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_real_s
     module procedure transpose_z_to_y_complex
     module procedure transpose_z_to_y_complex_s
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_real_s
     module procedure transpose_y_to_x_complex
     module procedure transpose_y_to_x_complex_s
  end interface transpose_y_to_x

  interface transpose_z_to_x
     module procedure transpose_z_to_x_real
     module procedure transpose_z_to_x_real_s
     module procedure transpose_z_to_x_complex
     module procedure transpose_z_to_x_complex_s
  end interface transpose_z_to_x

#ifdef OCC
  interface transpose_x_to_y_start
     module procedure transpose_x_to_y_real_start
     module procedure transpose_x_to_y_complex_start
  end interface transpose_x_to_y_start

  interface transpose_y_to_z_start
     module procedure transpose_y_to_z_real_start
     module procedure transpose_y_to_z_complex_start
  end interface transpose_y_to_z_start

  interface transpose_z_to_y_start
     module procedure transpose_z_to_y_real_start
     module procedure transpose_z_to_y_complex_start
  end interface transpose_z_to_y_start

  interface transpose_y_to_x_start
     module procedure transpose_y_to_x_real_start
     module procedure transpose_y_to_x_complex_start
  end interface transpose_y_to_x_start

  interface transpose_x_to_y_wait
     module procedure transpose_x_to_y_real_wait
     module procedure transpose_x_to_y_complex_wait
  end interface transpose_x_to_y_wait

  interface transpose_y_to_z_wait
     module procedure transpose_y_to_z_real_wait
     module procedure transpose_y_to_z_complex_wait
  end interface transpose_y_to_z_wait

  interface transpose_z_to_y_wait
     module procedure transpose_z_to_y_real_wait
     module procedure transpose_z_to_y_complex_wait
  end interface transpose_z_to_y_wait

  interface transpose_y_to_x_wait
     module procedure transpose_y_to_x_real_wait
     module procedure transpose_y_to_x_complex_wait
  end interface transpose_y_to_x_wait
#endif

  interface update_halo
     module procedure update_halo_real
     module procedure update_halo_complex
  end interface update_halo

  interface alloc_x
     module procedure alloc_x_real
     module procedure alloc_x_complex
  end interface alloc_x

  interface alloc_y
     module procedure alloc_y_real
     module procedure alloc_y_complex
  end interface alloc_y

  interface alloc_z
     module procedure alloc_z_real
     module procedure alloc_z_complex
  end interface alloc_z

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col
    logical, dimension(3), intent(IN), optional :: periodic_bc

    integer :: errorcode, ierror, row, col

    nx_global = nx
    ny_global = ny
    nz_global = nz

    if (present(periodic_bc)) then
       periodic_x = periodic_bc(1)
       periodic_y = periodic_bc(2)
       periodic_z = periodic_bc(3)
    else
       periodic_x = .false.
       periodic_y = .false.
       periodic_z = .false.
    end if

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call decomp_2d_abort(errorcode, &
               'Invalid 2D processor grid - nproc /= p_row*p_col')
       else
          row = p_row
          col = p_col
       end if
    end if

    ! Create 2D Catersian topology
    ! Note that in order to support periodic B.C. in the halo-cell code,
    ! need to create multiple topology objects: DECOMP_2D_COMM_CART_?,
    ! corresponding to three pencil orientations. They contain almost
    ! identical topological information but allow different combinations
    ! of periodic conditions.
    dims(1) = row
    dims(2) = col
    periodic(1) = periodic_y
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., &  ! do not reorder rank
         DECOMP_2D_COMM_CART_X, ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Y, ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_y
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Z, ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_CREATE")

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_COORDS")

    ! derive communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
         DECOMP_2D_COMM_COL,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_SUB")
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
         DECOMP_2D_COMM_ROW,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_SUB")

    ! gather information for halo-cell support code
    call init_neighbour

    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)

    ! make a copy of the decomposition information associated with the
    ! default global size in these global variables so applications can
    ! use them to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_SIZE")

    return
  end subroutine decomp_2d_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize

    implicit none

    call decomp_info_finalize(decomp_main)

    return
  end subroutine decomp_2d_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_decomp_info(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(OUT) :: decomp

    decomp = decomp_main

    return
  end subroutine get_decomp_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode

    ! verify the global size can actually be distributed as pencils
    if (nx<dims(1) .or. ny<dims(1) .or. ny<dims(2) .or. nz<dims(2)) then
       errorcode = 6
       call decomp_2d_abort(errorcode, &
            'Invalid 2D processor grid. ' // &
            'Make sure that min(nx,ny) >= p_row and ' // &
            'min(ny,nz) >= p_col')
    end if

    if (mod(nx,dims(1))==0 .and. mod(ny,dims(1))==0 .and. &
         mod(ny,dims(2))==0 .and. mod(nz,dims(2))==0) then
       decomp%even = .true.
    else
       decomp%even = .false.
    end if

    ! distribute mesh points
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
             decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1))
    allocate(decomp%x1st(0:dims(1)-1),decomp%x1en(0:dims(1)-1), &
             decomp%y1st(0:dims(1)-1),decomp%y1en(0:dims(1)-1))
    allocate(decomp%y2st(0:dims(2)-1),decomp%y2en(0:dims(2)-1), &
             decomp%z2st(0:dims(2)-1),decomp%z2en(0:dims(2)-1))
    call get_dist(nx,ny,nz,decomp)

    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), &
         decomp%zst, decomp%zen, decomp%zsz)

    ! prepare send/receive buffer displacement and count for ALLTOALL(V)
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
             decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1))
    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
             decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1))
    ! allocate arrays for MPI_ALLTOALLW calls
    ! x <=> z transpose
    allocate(decomp%xcnts_xz(nproc),decomp%zcnts_xz(nproc))
    allocate(decomp%xtypes_xzr(nproc),decomp%ztypes_xzr(nproc))
    allocate(decomp%xtypes_xzc(nproc),decomp%ztypes_xzc(nproc))
    allocate(decomp%xdispls_xz(nproc),decomp%zdispls_xz(nproc))
    ! x <=> y transpose
    allocate(decomp%xcnts_xy(nproc),decomp%zcnts_xy(nproc))
    allocate(decomp%xtypes_xyr(nproc),decomp%ztypes_xyr(nproc))
    allocate(decomp%xtypes_xyc(nproc),decomp%ztypes_xyc(nproc))
    allocate(decomp%xdispls_xy(nproc),decomp%zdispls_xy(nproc))
    ! y <=> z transpose
    allocate(decomp%xcnts_yz(nproc),decomp%zcnts_yz(nproc))
    allocate(decomp%xtypes_yzr(nproc),decomp%ztypes_yzr(nproc))
    allocate(decomp%xtypes_yzc(nproc),decomp%ztypes_yzc(nproc))
    allocate(decomp%xdispls_yz(nproc),decomp%zdispls_yz(nproc))
    !
    call prepare_buffer(decomp)

    return
  end subroutine decomp_info_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    ! Arguments
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! Local variables
    integer :: i, ierror

    deallocate(decomp%x1dist,decomp%y1dist)
    deallocate(decomp%y2dist,decomp%z2dist)
    deallocate(decomp%x1st,decomp%x1en,decomp%y1st,decomp%y1en)
    deallocate(decomp%y2st,decomp%y2en,decomp%z2st,decomp%z2en)
    deallocate(decomp%x1cnts,decomp%y1cnts)
    deallocate(decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%x1disp,decomp%y1disp)
    deallocate(decomp%y2disp,decomp%z2disp)

    ! *cnts_*
    deallocate(decomp%xcnts_xz,decomp%zcnts_xz)
    deallocate(decomp%xcnts_xy,decomp%zcnts_xy)
    deallocate(decomp%xcnts_yz,decomp%zcnts_yz)
    ! *types_*(1:nproc)
    do i = 1, nproc
      ! XZ
      if (decomp%xtypes_xzr(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_xzr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_xzr(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%ztypes_xzr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%xtypes_xzc(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_xzc(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_xzc(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%ztypes_xzc(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      ! XY
      if (decomp%xtypes_xyr(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_xyr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_xyr(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%ztypes_xyr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%xtypes_xyc(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_xyc(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_xyc(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%ztypes_xyc(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      ! YZ
      if (decomp%xtypes_yzr(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_yzr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_yzr(i).ne.MPI_INTEGER) then                       
        call MPI_Type_free(decomp%ztypes_yzr(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%xtypes_yzc(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_yzc(i),ierror)                      
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
      if (decomp%ztypes_yzc(i).ne.MPI_INTEGER) then                        
        call MPI_Type_free(decomp%ztypes_yzc(i),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_FREE")
      endif
    enddo
    ! *types_*r
    deallocate(decomp%xtypes_xzr,decomp%ztypes_xzr)
    deallocate(decomp%xtypes_xyr,decomp%ztypes_xyr)
    deallocate(decomp%xtypes_yzr,decomp%ztypes_yzr)
    ! *types_*c
    deallocate(decomp%xtypes_xzc,decomp%ztypes_xzc)
    deallocate(decomp%xtypes_xyc,decomp%ztypes_xyc)
    deallocate(decomp%xtypes_yzc,decomp%ztypes_yzc)
    ! *displs_*
    deallocate(decomp%xdispls_xz,decomp%zdispls_xz)
    deallocate(decomp%xdispls_xy,decomp%zdispls_xy)
    deallocate(decomp%xdispls_yz,decomp%zdispls_yz)
#ifdef MPI3
    ! *ranks_*
    deallocate(decomp%xranks_xz,decomp%zranks_xz)
    deallocate(decomp%xranks_xy,decomp%zranks_xy)
    deallocate(decomp%xranks_yz,decomp%zranks_yz)
    ! x <=> z
    call MPI_COMM_FREE(decomp%xtozNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(decomp%ztoxNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
    ! x <=> y
    call MPI_COMM_FREE(decomp%xtoyNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(decomp%ytoxNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
    ! y <=> z
    call MPI_COMM_FREE(decomp%ytozNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(decomp%ztoyNeighborComm,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
#endif
#ifdef SHM
    deallocate(decomp%x1disp_o,decomp%y1disp_o,decomp%y2disp_o, &
         decomp%z2disp_o)
    deallocate(decomp%x1cnts_s,decomp%y1cnts_s,decomp%y2cnts_s, &
         decomp%z2cnts_s)
    deallocate(decomp%x1disp_s,decomp%y1disp_s,decomp%y2disp_s, &
         decomp%z2disp_s)
#endif

    return
  end subroutine decomp_info_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for statistic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statS(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipS = i_skip
    jskipS = j_skip
    kskipS = k_skip

    skip(1)=iskipS
    skip(2)=jskipS
    skip(3)=kskipS

    do i=1,3
       if (from1) then
          xstS(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstS(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = xend(i)/skip(i)
       end if
       xszS(i) = xenS(i)-xstS(i)+1
    end do

    do i=1,3
       if (from1) then
          ystS(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystS(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = yend(i)/skip(i)
       end if
       yszS(i) = yenS(i)-ystS(i)+1
    end do

    do i=1,3
       if (from1) then
          zstS(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstS(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = zend(i)/skip(i)
       end if
       zszS(i) = zenS(i)-zstS(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for visualization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statV(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipV = i_skip
    jskipV = j_skip
    kskipV = k_skip

    skip(1)=iskipV
    skip(2)=jskipV
    skip(3)=kskipV

    do i=1,3
       if (from1) then
          xstV(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstV(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = xend(i)/skip(i)
       end if
       xszV(i) = xenV(i)-xstV(i)+1
    end do

    do i=1,3
       if (from1) then
          ystV(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystV(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = yend(i)/skip(i)
       end if
       yszV(i) = yenV(i)-ystV(i)+1
    end do

    do i=1,3
       if (from1) then
          zstV(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstV(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = zend(i)/skip(i)
       end if
       zszV(i) = zenV(i)-zstV(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for probe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statP(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipP = i_skip
    jskipP = j_skip
    kskipP = k_skip

    skip(1)=iskipP
    skip(2)=jskipP
    skip(3)=kskipP

    do i=1,3
       if (from1) then
          xstP(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstP(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = xend(i)/skip(i)
       end if
       xszP(i) = xenP(i)-xstP(i)+1
    end do

    do i=1,3
       if (from1) then
          ystP(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystP(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = yend(i)/skip(i)
       end if
       yszP(i) = yenP(i)-ystP(i)+1
    end do

    do i=1,3
       if (from1) then
          zstP(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstP(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = zend(i)/skip(i)
       end if
       zszP(i) = zenP(i)-zstP(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statP

  ! Copy data from a fine-resolution array to a coarse one for statistic
  subroutine fine_to_coarseS(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystS(1):yenS(1),ystS(2):yenS(2),ystS(3):yenS(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstS(1):zenS(1),zstS(2):zenS(2),zstS(3):zenS(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseS

  ! Copy data from a fine-resolution array to a coarse one for visualization
  subroutine fine_to_coarseV(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystV(1):yenV(1),ystV(2):yenV(2),ystV(3):yenV(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstV(1):zenV(1),zstV(2):zenV(2),zstV(3):zenV(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseV

  ! Copy data from a fine-resolution array to a coarse one for probe
  subroutine fine_to_coarseP(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystP(1):yenP(1),ystP(2):yenP(2),ystP(3):yenP(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstP(1):zenP(1),zstP(2):zenP(2),zstP(3):zenP(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find sub-domain information held by current processor
  !   INPUT: 
  !     nx, ny, nz - global data dimension
  !     pdim(3)    - number of processor grid in each dimension, 
  !                  valid values: 1 - distibute locally; 
  !                                2 - distribute across p_row; 
  !                                3 - distribute across p_col
  !   OUTPUT:
  !     lstart(3)  - starting index
  !     lend(3)    - ending index
  !     lsize(3)   - size of the sub-block (redundant) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3

       if (i==1) then
          gsize = nx
       else if (i==2) then
          gsize = ny
       else if (i==3) then
          gsize = nz
       end if

       if (pdim(i) == 1) then        ! all local
          lstart(i) = 1
          lend(i)   = gsize
          lsize(i)  = gsize
       elseif (pdim(i) == 2) then    ! distribute across dims(1)
          allocate(st(0:dims(1)-1))
          allocate(en(0:dims(1)-1))
          allocate(sz(0:dims(1)-1))
          call distribute(gsize,dims(1),st,en,sz)
          lstart(i) = st(coord(1))
          lend(i)   = en(coord(1))
          lsize(i)  = sz(coord(1))
          deallocate(st,en,sz)
       elseif (pdim(i) == 3) then    ! distribute across dims(2)
          allocate(st(0:dims(2)-1))
          allocate(en(0:dims(2)-1))
          allocate(sz(0:dims(2)-1))
          call distribute(gsize,dims(2),st,en,sz)
          lstart(i) = st(coord(2))
          lend(i)   = en(coord(2))
          lsize(i)  = sz(coord(2))
          deallocate(st,en,sz)
       end if

    end do
    return   

  end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine distribute(data1,proc,st,en,sz)

    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nu

    size1=data1/proc
    nu = data1 - size1 * proc
    nl = proc - nu
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
       st(i) = st(i-1) + size1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,proc-1
       st(i) = en(i-1) + 1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1

    return
  end subroutine distribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist(nx,ny,nz,decomp)

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    call distribute(nx,dims(1),decomp%x1st,decomp%x1en,decomp%x1dist)
    call distribute(ny,dims(1),decomp%y1st,decomp%y1en,decomp%y1dist)

    call distribute(ny,dims(2),decomp%y2st,decomp%y2en,decomp%y2dist)
    call distribute(nz,dims(2),decomp%z2st,decomp%z2en,decomp%z2dist)

    return
  end subroutine get_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_buffer(decomp)

    implicit none

    ! Arguments
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! Local variables
    integer :: i, k, rk, rank_x, rank_z, subsize, offset, ierror
    integer, dimension(2) :: coord
#ifdef MPI3
    integer, dimension(nproc) :: xranks, zranks, xweights, zweights
    integer :: index_src, index_dest
#endif

    ! MPI_ALLTOALLV buffer information

    do i=0, dims(1)-1
       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
       if (i==0) then
          decomp%x1disp(i) = 0  ! displacement is 0-based index
          decomp%y1disp(i) = 0
       else
          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
       end if
    end do

    do i=0, dims(2)-1
       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
       if (i==0) then
          decomp%y2disp(i) = 0  ! displacement is 0-based index
          decomp%z2disp(i) = 0
       else
          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
       end if
    end do

    ! MPI_ALLTOALL buffer information

    ! For evenly distributed data, following is an easier implementation.
    ! But it should be covered by the more general formulation below.
    !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
    !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
    !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
    !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    ! for X <=> Y transposes
    decomp%x1count = decomp%x1dist(dims(1)-1) * &
                     decomp%y1dist(dims(1)-1) * decomp%xsz(3)
    decomp%y1count = decomp%x1count
    ! for Y <=> Z transposes
    decomp%y2count = decomp%y2dist(dims(2)-1) * &
                     decomp%z2dist(dims(2)-1) * decomp%zsz(1)
    decomp%z2count = decomp%y2count

    ! Information for MPI_Alltoallw
    ! x <=> z transpose
    decomp%xdispls_xz(:)=0
    decomp%zdispls_xz(:)=0
    decomp%xcnts_xz(:)=0
    decomp%zcnts_xz(:)=0
    decomp%xtypes_xzr(:)=MPI_INTEGER
    decomp%ztypes_xzr(:)=MPI_INTEGER
    decomp%xtypes_xzc(:)=MPI_INTEGER
    decomp%ztypes_xzc(:)=MPI_INTEGER
    ! x <=> y transpose
    decomp%xdispls_xy(:)=0
    decomp%zdispls_xy(:)=0
    decomp%xcnts_xy(:)=0                                                                                 
    decomp%zcnts_xy(:)=0
    decomp%xtypes_xyr(:)=MPI_INTEGER
    decomp%ztypes_xyr(:)=MPI_INTEGER
    decomp%xtypes_xyc(:)=MPI_INTEGER                                                                     
    decomp%ztypes_xyc(:)=MPI_INTEGER
    ! y <=> z transpose
    decomp%xdispls_yz(:)=0
    decomp%zdispls_yz(:)=0
    decomp%xcnts_yz(:)=0
    decomp%zcnts_yz(:)=0
    decomp%xtypes_yzr(:)=MPI_INTEGER
    decomp%ztypes_yzr(:)=MPI_INTEGER
    decomp%xtypes_yzc(:)=MPI_INTEGER
    decomp%ztypes_yzc(:)=MPI_INTEGER

    ! Init local variables (x <=> z transpose)
#ifdef MPI3
    xranks(:) = 0
    zranks(:) = 0
    xweights(:) = 0
    zweights(:) = 0
    index_src=0
    index_dest=0
#endif

    do k=0,dims(1)-1
    do i=0,dims(2)-1

      ! Get rank_x and rank_z
      call MPI_Cart_rank(DECOMP_2D_COMM_CART_X,(/k,i/),rank_x,ierror)
      if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_RANK")
      call MPI_Cart_rank(DECOMP_2D_COMM_CART_Z,(/k,i/),rank_z,ierror)
      if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_RANK")
      ! Safety check
      if (rank_x.eq.rank_z) then
        rk = rank_x
      else
        print *, "Rank ", nrank, " : error in prepare_buffer."
        print *, "Rank ", nrank, " : [rank_x, rank_z] = [",rank_x, rank_z,"]"
        call decomp_2d_abort(13, "prepare_buffer: incompatible ranks.")
      endif

      !
      ! Local data
      !   X pencil : [1,nx]x[xst(2),xen(2)]x[xst(3),xen(3)]
      !   Y pencil : [yst(1),yen(1)]x[1,ny]x[yst(3),yen(3)]
      !   Z pencil : [zst(1),zen(1)]x[zst(2),zen(2)]x[1,nz]
      ! Remote data on CPU rk = rank_x = rank_z located at (k,i)
      !   X pencil : [1,nx]x[y1st(rk),y1en(rk)]x[z2st(rk),z2en(rk)]
      !   Y pencil : [x1st(rk),x1en(rk)]x[1,ny]x[z2st(rk),z2en(rk)]
      !   Z pencil : [x1st(rk),x1en(rk)]x[y2st(rk),y2en(rk)]x[1,nz]
      !

      !
      ! Transpose X <=> Z
      ! No checks on X or Z dimension
      ! Transform from Z into X pencils, so these always overlap
      !
      ! First define the MPI subarray for Z pencil
      if (decomp%zst(2).le.decomp%y1en(k) .and. &
          decomp%zen(2).ge.decomp%y1st(k)) then

        ! Safety check
        if (decomp%ztypes_xzr(rk+1).ne.MPI_INTEGER .or. &
            decomp%ztypes_xzc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer." 
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%zcnts_xz(rk+1)=1
        subsize=min(decomp%zen(2),decomp%y1en(k))-max(decomp%zst(2),decomp%y1st(k))+1
        offset =max(decomp%zst(2),decomp%y1st(k))-decomp%zst(2)

#ifdef MPI3
        index_src=index_src+1
        zranks(index_src)=rk
        zweights(index_src)=decomp%zsz(1)*subsize*decomp%z2dist(i)
#endif

        call MPI_Type_create_subarray(3,decomp%zsz, &
               (/decomp%zsz(1),subsize,decomp%z2dist(i)/), &
               (/0,offset,decomp%z2st(i)-decomp%zst(3)/), &
               MPI_ORDER_FORTRAN,real_type,decomp%ztypes_xzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_xzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%zsz, &
               (/decomp%zsz(1),subsize,decomp%z2dist(i)/), &
               (/0,offset,decomp%z2st(i)-decomp%zst(3)/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%ztypes_xzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_xzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

      ! Then define the MPI subarray for X pencil
      if (decomp%xst(2).le.decomp%y2en(i) .and. &
          decomp%xen(2).ge.decomp%y2st(i)) then

        ! Safety check
        if (decomp%xtypes_xzr(rk+1).ne.MPI_INTEGER .or. &
            decomp%xtypes_xzc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer."
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%xcnts_xz(rk+1)=1
        subsize=min(decomp%xen(2),decomp%y2en(i))-max(decomp%xst(2),decomp%y2st(i))+1
        offset =max(decomp%xst(2),decomp%y2st(i))-decomp%xst(2)

#ifdef MPI3
        index_dest=index_dest+1
        xranks(index_dest)=rk
        xweights(index_dest)=decomp%x1dist(k)*subsize*decomp%xsz(3)
#endif

        call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),subsize,decomp%xsz(3)/), &
               (/decomp%x1st(k)-decomp%xst(1),offset,0/), &
               MPI_ORDER_FORTRAN,real_type,decomp%xtypes_xzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_xzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),subsize,decomp%xsz(3)/), &
               (/decomp%x1st(k)-decomp%xst(1),offset,0/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%xtypes_xzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_xzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

    enddo
    enddo

#ifdef MPI3
    allocate(decomp%xranks_xz(index_dest))
    allocate(decomp%zranks_xz(index_src))

    decomp%xranks_xz=xranks(1:index_dest)+1
    decomp%zranks_xz=zranks(1:index_src)+1

    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      MPI_INFO_NULL,.true.,decomp%xtozNeighborComm,ierror)
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      MPI_INFO_NULL,.true.,decomp%ztoxNeighborComm,ierror)
#endif

    ! Init local variables (x <=> y transpose)
#ifdef MPI3
    xranks(:) = 0
    zranks(:) = 0
    xweights(:) = 0
    zweights(:) = 0
    index_src=0
    index_dest=0
#endif

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_COORDS")
    i = coord(2)
    do k=0,dims(1)-1

      ! Get rank_x and rank_z
      call MPI_Cart_rank(DECOMP_2D_COMM_CART_X,(/k,i/),rk,ierror)
      if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_RANK")

      !
      ! Local data
      !   X pencil : [1,nx]x[xst(2),xen(2)]x[xst(3),xen(3)]
      !   Y pencil : [yst(1),yen(1)]x[1,ny]x[yst(3),yen(3)]
      !   Z pencil : [zst(1),zen(1)]x[zst(2),zen(2)]x[1,nz]
      ! Remote data on CPU rk located at (k,i)
      !   X pencil : [1,nx]x[y1st(rk),y1en(rk)]x[z2st(rk),z2en(rk)]
      !   Y pencil : [x1st(rk),x1en(rk)]x[1,ny]x[z2st(rk),z2en(rk)]
      !   Z pencil : [x1st(rk),x1en(rk)]x[y2st(rk),y2en(rk)]x[1,nz]
      !

      !
      ! Transpose X <=> Y
      ! No checks on X or Y dimension
      ! Transform from Y into X pencils, so these always overlap
      !
      ! First define the MPI subarray for Y pencil
      if (decomp%yst(3).le.decomp%z2en(i) .and. &
          decomp%yen(3).ge.decomp%z2st(i)) then

        ! Safety check
        if (decomp%ztypes_xyr(rk+1).ne.MPI_INTEGER .or. &
            decomp%ztypes_xyc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer."
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%zcnts_xy(rk+1)=1
        subsize=min(decomp%yen(3),decomp%z2en(i))-max(decomp%yst(3),decomp%z2st(i))+1
        offset =max(decomp%yst(3),decomp%z2st(i))-decomp%yst(3)

#ifdef MPI3
        index_src=index_src+1
        zranks(index_src)=rk
        zweights(index_src)=decomp%ysz(1)*decomp%y1dist(k)*subsize
#endif

        call MPI_Type_create_subarray(3,decomp%ysz, &
               (/decomp%ysz(1),decomp%y1dist(k),subsize/), &
               (/0,decomp%y1st(k)-decomp%yst(2),offset/), &
               MPI_ORDER_FORTRAN,real_type,decomp%ztypes_xyr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_xyr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%ysz, &
               (/decomp%ysz(1),decomp%y1dist(k),subsize/), &
               (/0,decomp%y1st(k)-decomp%yst(2),offset/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%ztypes_xyc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_xyc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

      ! Then define the MPI subarray for X pencil
      if (decomp%xst(3).le.decomp%z2en(i) .and. &
          decomp%xen(3).ge.decomp%z2st(i)) then

        ! Safety check
        if (decomp%xtypes_xyr(rk+1).ne.MPI_INTEGER .or. &
            decomp%xtypes_xyc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer."
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%xcnts_xy(rk+1)=1
        subsize=min(decomp%xen(3),decomp%z2en(i))-max(decomp%xst(3),decomp%z2st(i))+1
        offset =max(decomp%xst(3),decomp%z2st(i))-decomp%xst(3)

#ifdef MPI3
        index_dest=index_dest+1
        xranks(index_dest)=rk
        xweights(index_dest)=decomp%x1dist(k)*decomp%xsz(2)*subsize
#endif

        call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),decomp%xsz(2),subsize/), &
               (/decomp%x1st(k)-decomp%xst(1),0,offset/), &
               MPI_ORDER_FORTRAN,real_type,decomp%xtypes_xyr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_xyr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),decomp%xsz(2),subsize/), &
               (/decomp%x1st(k)-decomp%xst(1),0,offset/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%xtypes_xyc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_xyc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

    enddo

#ifdef MPI3
    allocate(decomp%xranks_xy(index_dest))
    allocate(decomp%zranks_xy(index_src))

    decomp%xranks_xy=xranks(1:index_dest)+1
    decomp%zranks_xy=zranks(1:index_src)+1

    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      MPI_INFO_NULL,.true.,decomp%xtoyNeighborComm,ierror)
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      MPI_INFO_NULL,.true.,decomp%ytoxNeighborComm,ierror)
#endif

    ! Init local variables (y <=> z transpose)
#ifdef MPI3
    xranks(:) = 0
    zranks(:) = 0
    xweights(:) = 0
    zweights(:) = 0
    index_src=0
    index_dest=0
#endif
      
    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_COORDS")
    k = coord(1)
    do i=0,dims(2)-1

      ! Get rank_x and rank_z
      call MPI_Cart_rank(DECOMP_2D_COMM_CART_X,(/k,i/),rk,ierror)
      if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_RANK")

      !   
      ! Local data
      !   X pencil : [1,nx]x[xst(2),xen(2)]x[xst(3),xen(3)]
      !   Y pencil : [yst(1),yen(1)]x[1,ny]x[yst(3),yen(3)]
      !   Z pencil : [zst(1),zen(1)]x[zst(2),zen(2)]x[1,nz]
      ! Remote data on CPU rk located at (k,i)
      !   X pencil : [1,nx]x[y1st(rk),y1en(rk)]x[z2st(rk),z2en(rk)]
      !   Y pencil : [x1st(rk),x1en(rk)]x[1,ny]x[z2st(rk),z2en(rk)]
      !   Z pencil : [x1st(rk),x1en(rk)]x[y2st(rk),y2en(rk)]x[1,nz]
      !

      !
      ! Transpose Y <=> Z
      ! No checks on Y or Z dimension
      ! Transform from Y into Z pencils, so these always overlap
      !
      ! First define the MPI subarray for Z pencil
      if (decomp%zst(1).le.decomp%x1en(k) .and. &
          decomp%zen(1).ge.decomp%x1st(k)) then

        ! Safety check
        if (decomp%ztypes_yzr(rk+1).ne.MPI_INTEGER .or. &
            decomp%ztypes_yzc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer." 
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%zcnts_yz(rk+1)=1
        subsize=min(decomp%zen(1),decomp%x1en(k))-max(decomp%zst(1),decomp%x1st(k))+1
        offset =max(decomp%zst(1),decomp%x1st(k))-decomp%zst(1)

#ifdef MPI3
        index_src=index_src+1
        zranks(index_src)=rk
        zweights(index_src)=subsize*decomp%zsz(2)*decomp%z2dist(i)
#endif

        call MPI_Type_create_subarray(3,decomp%zsz, &
               (/subsize,decomp%zsz(2),decomp%z2dist(i)/), &
               (/offset,0,decomp%z2st(i)-decomp%zst(3)/), &
               MPI_ORDER_FORTRAN,real_type,decomp%ztypes_yzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_yzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%zsz, &
               (/subsize,decomp%zsz(2),decomp%z2dist(i)/), &
               (/offset,0,decomp%z2st(i)-decomp%zst(3)/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%ztypes_yzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%ztypes_yzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

      ! Then define the MPI subarray for Y pencil
      if (decomp%yst(1).le.decomp%x1en(k) .and. &
          decomp%yen(1).ge.decomp%x1st(k)) then

        ! Safety check
        if (decomp%xtypes_yzr(rk+1).ne.MPI_INTEGER .or. &
            decomp%xtypes_yzc(rk+1).ne.MPI_INTEGER) then
          print *, "Rank ", nrank, " : error in prepare_buffer."
          call decomp_2d_abort(13, "prepare_buffer: collision detected.")
        endif

        decomp%xcnts_yz(rk+1)=1
        subsize=min(decomp%yen(1),decomp%x1en(k))-max(decomp%yst(1),decomp%x1st(k))+1
        offset =max(decomp%yst(1),decomp%x1st(k))-decomp%yst(1)

#ifdef MPI3
        index_dest=index_dest+1
        xranks(index_dest)=rk
        xweights(index_dest)=subsize*decomp%y2dist(i)*decomp%xsz(3)
#endif

        call MPI_Type_create_subarray(3,decomp%ysz, &
               (/subsize,decomp%y2dist(i),decomp%xsz(3)/), &
               (/offset,decomp%y2st(i)-decomp%yst(2),0/), &
               MPI_ORDER_FORTRAN,real_type,decomp%xtypes_yzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_yzr(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

        call MPI_Type_create_subarray(3,decomp%ysz, &
               (/subsize,decomp%y2dist(i),decomp%xsz(3)/), &
               (/offset,decomp%y2st(i)-decomp%yst(2),0/), &
               MPI_ORDER_FORTRAN,complex_type,decomp%xtypes_yzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_CREATE_SUBARRAY")
        call MPI_Type_commit(decomp%xtypes_yzc(rk+1),ierror)
        if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_TYPE_COMMIT")

      endif

    enddo

#ifdef MPI3
    allocate(decomp%xranks_yz(index_dest))
    allocate(decomp%zranks_yz(index_src))

    decomp%xranks_yz=xranks(1:index_dest)+1
    decomp%zranks_yz=zranks(1:index_src)+1

    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      MPI_INFO_NULL,.true.,decomp%ytozNeighborComm,ierror)
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      MPI_INFO_NULL,.true.,decomp%ztoyNeighborComm,ierror)
#endif

    return
  end subroutine prepare_buffer

#ifdef OCC
  ! For non-blocking communication code, progress the comminication stack
  subroutine transpose_test(handle)

    implicit none

    integer :: handle, ierror

    call NBC_TEST(handle,ierror)

    return
  end subroutine transpose_test
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transposition routines 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "transpose_x_to_y.inc"
#include "transpose_y_to_z.inc"
#include "transpose_x_to_z.inc"
#include "transpose_z_to_y.inc"
#include "transpose_y_to_x.inc"
#include "transpose_z_to_x.inc"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    double precision :: t1, t2, best_time
    integer :: nfact, i, row, col, ierror, errorcode

    real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3

    TYPE(DECOMP_INFO) :: decomp

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    best_time = huge(t1)
    best_p_row = -1
    best_p_col = -1

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    do i=1, nfact

       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx_global,ny_global)>=row .and. &
            min(ny_global,nz_global)>=col) then

          ! 2D Cartesian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .false.
          periodic(2) = .false.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_CREATE")
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false., DECOMP_2D_COMM_CART_Z, ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_CREATE")
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_COORDS")

          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_SUB")
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_CART_SUB")

          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          allocate(u3(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_z(u2,u3,decomp)
          call transpose_z_to_y(u3,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2,u3)
          call decomp_info_finalize(decomp)

          call MPI_COMM_FREE(DECOMP_2D_COMM_COL,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
          call MPI_COMM_FREE(DECOMP_2D_COMM_ROW,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
          call MPI_COMM_FREE(DECOMP_2D_COMM_CART_X,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")
          call MPI_COMM_FREE(DECOMP_2D_COMM_CART_Z,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_COMM_FREE")

          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
               MPI_COMM_WORLD,ierror)
          if (ierror.ne.0) call decomp_2d_abort(ierror, "MPI_ALLREDUCE")
          t1 = t1 / dble(nproc)

          if (nrank==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          end if

          if (best_time > t1) then
             best_time = t1
             best_p_row = row
             best_p_col = col
          end if

       end if

    end do ! loop through processer grid

    deallocate(factors)

    if (best_p_row/=-1) then
       if (nrank==0) then
          write(*,*) 'the best processor grid is probably ', &
               best_p_row, ' by ', best_p_col
       end if
    else
       errorcode = 9
       call decomp_2d_abort(errorcode, &
            'The processor-grid auto-tuning code failed. ' // &
            'The number of processes requested is probably too large.')
    end if

    return
  end subroutine best_2d_grid

#include "factor.inc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.inc"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_abort(errorcode, msg)

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    integer :: ierror

    if (nrank==0) then
       write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
    end if
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)

    return
  end subroutine decomp_2d_abort


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.inc"


end module decomp_2d

