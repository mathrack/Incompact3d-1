!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
module visu

  implicit none

  ! True to activate the XDMF output
  logical, save :: use_xdmf = .true.
  ! True to use the new enumeration
  logical, save :: filenamedigits = .true.
  ! output2D is defined in the input.i3d file
  !        0 for 3D output (default)
  !        1 for 2D output with X average
  !        2 for 2D output with Y average
  !        3 for 2D output with Z average
  integer, save :: output2D
  integer :: ioxdmf
  character(len=9) :: ifilenameformat = '(I9.9)'

  private
  public :: output2D, visu_init, write_snapshot

contains

  !
  ! Initialize the visu module
  !
  subroutine visu_init()

#ifdef MPI3
    USE MPI_f08
#else
    USE MPI
#endif
    use param, only : ilmn, iscalar, ilast, ifirst, ioutput, istret
    use variables, only : numscalar, prec, nvisu
    use decomp_2d, only : nrank, mytype, xszV, yszV, zszV

    implicit none

    ! Local variables
    logical :: dir_exists
    integer :: noutput, nsnapout
    real(mytype) :: memout

    ! Create folder if needed
    if (nrank==0) then
      inquire(file="data", exist=dir_exists)
      if (.not.dir_exists) then
        call system("mkdir data 2> /dev/null")
      end if
    end if

    ! HDD usage of visu module
    if (nrank==0) then
      noutput = (ilast - ifirst +1)/ioutput

      nsnapout = 4
      if (ilmn)         nsnapout = nsnapout + 1
      if (iscalar.ne.0) nsnapout = nsnapout + numscalar

      memout = prec * nsnapout * noutput
      if (output2D.eq.0) then
        memout = memout * xszV(1) * yszV(2) * zszV(3)
      else if (output2D.eq.1) then
        memout = memout *           yszV(2) * zszV(3)
      else if (output2D.eq.2) then
        memout = memout * xszV(1)           * zszV(3)
      else if (output2D.eq.3) then
        memout = memout * xszV(1) * yszV(2)
      endif
      print *,'==========================================================='
      print *,'Visu module requires ',real(memout*1e-9,4),'GB'
      print *,'==========================================================='
    end if

    ! Safety check
    if (output2D.lt.0 .or. output2D.gt.3 &
        .or. (output2d.eq.2.and.istret.ne.0)) then
      if (nrank.eq.0) print *, "Visu module: incorrect value for output2D."
      call MPI_ABORT(MPI_COMM_WORLD, 0, noutput)
      stop
    endif

  end subroutine visu_init

  !
  ! Write a snapshot
  !
  subroutine write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)

    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : nrank
    use decomp_2d, only : fine_to_coarsev
    use decomp_2d_io, only : decomp_2d_write_one

    use param, only : ivisu, ioutput, nrhotime, ilmn, iscalar, iibm, istret
    use param, only : xlx, yly, zlz

    use variables, only : derx, dery, derz
    use variables, only : ffx, ffxp, fsx, fsxp, fwx, fwxp
    use variables, only : ffy, ffyp, fsy, fsyp, fwy, fwyp, ppy
    use variables, only : ffz, ffzp, fsz, fszp, fwz, fwzp
    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : nx, ny, nz, beta, numscalar

    use var, only : one
    use var, only : uvisu
    use var, only : pp1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1, nxmsize
    use var, only : pp2, ta2, tb2, tc2, td2, te2, tf2, ppi2, di2, dip2, ph2, nymsize
    use var, only : ppi3, ta3, tb3, tc3, td3, te3, tf3, di3, dip3, ph3, nzmsize
    use var, only : npress
    use var, only : dt,t

    use tools, only : rescale_pressure

    implicit none

    !! inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    integer, intent(in) :: itime

    ! Local variables
    integer :: is
    character(len=32) :: fmt2, fmt3, fmt4, num
    real :: tstart, tend

    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tstart)
      print *,'Writing snapshots =>',itime/ioutput
    end if

    ! Snapshot number
    if (filenamedigits) then
      ! New enumeration system, it works integrated with xcompact3d_toolbox
      write(num, ifilenameformat) itime
    else
      ! Classic enumeration system
      write(num, ifilenameformat) itime/ioutput
    endif

    ! Write XDMF header
    if (use_xdmf) call write_xdmf_header(".", "snapshot", trim(num))

    ! Write velocity
    call write_field(ux1, ".", "ux", trim(num))
    call write_field(uy1, ".", "uy", trim(num))
    call write_field(uz1, ".", "uz", trim(num))

    ! Interpolate pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
    (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)

    ! Rescale pressure
    call rescale_pressure(ta1)

    ! Write pressure
    call write_field(ta1, ".", "pp", trim(num), .true.)

    ! LMN - write density
    if (ilmn) call write_field(rho1(:,:,:,1), ".", "rho", trim(num), .true.)

    ! Write scalars
    if (iscalar.ne.0) then
      do is = 1, numscalar
        call write_field(phi1(:,:,:,is), ".", "phi"//char(48+is), trim(num), .true.)
      enddo
    endif

    !
    ! At this location, we could add extra variables in the XDMF file
    !

    ! Write XDMF footer
    if (use_xdmf) call write_xdmf_footer()

    ! Add metadata if XDMF is not used
    if (.not. use_xdmf) then
      if (nrank==0) then

        write(fmt2,'("(A,I16)")')
        write(fmt3,'("(A,F16.4)")')
        write(fmt4,'("(A,F16.12)")')

        open(newunit=is,file="./data/snap"//trim(num)//".ini",action='write',status='replace')
        write(is,'(A)')'[domain]'
        write(is,fmt2) 'nx=      ',nx
        write(is,fmt2) 'ny=      ',ny
        write(is,fmt2) 'nz=      ',nz
        write(is,fmt2) 'istret=  ',istret
        write(is,fmt4) 'beta=    ',beta
        write(is,fmt3) 'Lx=      ',xlx
        write(is,fmt3) 'Ly=      ',yly
        write(is,fmt3) 'Lz=      ',zlz
        write(is,'(A)')'[time]'
        write(is,fmt2) 'itime=   ',itime
        write(is,fmt3) 'dt=      ',dt
        write(is,fmt3) 't =      ',t
        close(is)

      endif
    endif

    ! Update log file
    if (nrank.eq.0) then
      call cpu_time(tend)
      write(*,'(" Time for writing snapshots (s): ",F12.8)') tend-tstart
    endif

  end subroutine write_snapshot

  !
  ! Output binary data associated with the pressure
  !
  subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,ta3,di3,nxmsize,nymsize,nzmsize,uvisu,pre1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use var, only : zero
    use tools, only : mean_plane_z

    implicit none

    integer :: nxmsize,nymsize,nzmsize

    real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
    !Z PENCILS NXM NYM NZM-->NXM NYM NZ
    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
    !Y PENCILS NXM NYM NZ -->NXM NY NZ
    real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
    real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
    !X PENCILS NXM NY NZ  -->NX NY NZ
    real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1,pre1

    character(len=30) filename

    !WORK Z-PENCILS
    call interzpv(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
    call interypv(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
    call interxpv(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    pre1=tb1

    if (save_pre.eq.1) then
      uvisu = zero
      call fine_to_coarseV(1,pre1,uvisu)
      write(filename,"('./data/pre',I4.4)") itime/ioutput
      call decomp_2d_write_one(1,uvisu,filename,2)
    endif

    if (save_prem.eq.1) then
      write(filename,"('./data/prem',I4.4)") itime/ioutput
      call decomp_2d_write_plane(1,pre1,3,-1,filename)
    endif

    return

  end subroutine VISU_PRE

  !
  ! Write the header of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_header(pathname, filename, num)

    use variables, only : nvisu, yp
    use param, only : dx,dy,dz,istret
    use decomp_2d, only : mytype, nrank, xszV, yszV, zszV, ystV

    implicit none

    ! Arguments
    character(len=*), intent(in) :: pathname, filename, num

    ! Local variables
    integer :: i,k
    real(mytype) :: xp(xszV(1)), zp(zszV(3))

    if (nrank.eq.0) then
      OPEN(newunit=ioxdmf,file="./data/"//pathname//"/"//filename//'-'//num//'.xdmf')

      write(ioxdmf,'(A22)')'<?xml version="1.0" ?>'
      write(ioxdmf,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf,*)'<Domain>'
      if (istret.ne.0) then
        do i=1,xszV(1)
          xp(i) = real(i-1,mytype)*dx*nvisu
        enddo
        do k=1,zszV(3)
          zp(k) = real(k-1,mytype)*dz*nvisu
        enddo
        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),xszV(1),'">'
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),1,'">'
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),1,xszV(1),'">'
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        Dimensions="',1,yszV(2),xszV(1),'">'
        endif
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="VXVYVZ">'
        if (output2D.ne.1) then
          write(ioxdmf,*)'        <DataItem Dimensions="',xszV(1),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.2) then
          write(ioxdmf,*)'        <DataItem Dimensions="',yszV(2),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(ystV(1)::nvisu)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',yp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        if (output2D.ne.3) then
          write(ioxdmf,*)'        <DataItem Dimensions="',zszV(3),'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(:)
        else
          write(ioxdmf,*)'        <DataItem Dimensions="1" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',zp(1)
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      else
        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),xszV(1),'">'
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),yszV(2),1,'">'
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        Dimensions="',zszV(3),1,xszV(1),'">'
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        Dimensions="',1,yszV(2),xszV(1),'">'
        endif
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        write(ioxdmf,*)'        <!-- Origin -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        write(ioxdmf,*)'        0.0 0.0 0.0'
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <!-- DxDyDz -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        if (output2D.eq.0) then
          write(ioxdmf,*)'        ',nvisu*dz,nvisu*dy,nvisu*dx
        else if (output2D.eq.1) then
          write(ioxdmf,*)'        ',dz,dy,1.
        else if (output2D.eq.2) then
          write(ioxdmf,*)'        ',dz,1.,dx
        else if (output2D.eq.3) then
          write(ioxdmf,*)'        ',1.,dy,dx
        endif
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      endif
      write(ioxdmf,*)'    <Grid Name="'//num//'" GridType="Uniform">'
      write(ioxdmf,*)'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(ioxdmf,*)'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    endif
  end subroutine write_xdmf_header

  !
  ! Write the footer of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_footer()

    use decomp_2d, only : nrank

    implicit none

    if (nrank.eq.0) then
      write(ioxdmf,'(/)')
      write(ioxdmf,*)'    </Grid>'
      write(ioxdmf,*)'</Domain>'
      write(ioxdmf,'(A7)')'</Xdmf>'
      close(ioxdmf)
    endif

  end subroutine write_xdmf_footer

  !
  ! Write the given field for visualization
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_field(f1, pathname, filename, num, skip_ibm)

    use var, only : di1, ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
    use decomp_2d, only : mytype, xsize, xszV, yszV, zszV
    use decomp_2d, only : nrank, fine_to_coarseV
    use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_write_plane

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: f1
    character(len=*), intent(in) :: pathname, filename, num 
    logical, optional, intent(in) :: skip_ibm

    if (use_xdmf) then
      if (nrank.eq.0) then
        write(ioxdmf,*)'        <Attribute Name="'//filename//'" Center="Node">'
        write(ioxdmf,*)'           <DataItem Format="Binary"'
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
        if (output2D.eq.0) then
          write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
        else
          write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
        endif
#else
        write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
#endif
#else
        write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
#endif
        if (output2D.eq.0) then
          write(ioxdmf,*)'            Dimensions="',zszV(3),yszV(2),xszV(1),'">'
        else if (output2D.eq.1) then
          write(ioxdmf,*)'            Dimensions="',zszV(3),yszV(2),1,'">'
        else if (output2D.eq.2) then
          write(ioxdmf,*)'            Dimensions="',zszV(3),1,xszV(1),'">'
        else if (output2D.eq.3) then
          write(ioxdmf,*)'            Dimensions="',1,yszV(2),xszV(1),'">'
        endif
        write(ioxdmf,*)'              ./'//pathname//"/"//filename//'-'//num//'.bin'
        write(ioxdmf,*)'           </DataItem>'
        write(ioxdmf,*)'        </Attribute>'
      endif
    endif

    if (iibm==2 .and. .not.present(skip_ibm)) then
      di1(:,:,:) = (one - ep1(:,:,:)) * f1(:,:,:)
    else
      di1(:,:,:) = f1(:,:,:)
    endif

    if (output2D.eq.0) then
      uvisu = zero
      call fine_to_coarseV(1,di1,uvisu)
      call decomp_2d_write_one(1,uvisu,"./data/"//pathname//'/'//filename//'-'//num//'.bin',2)
    else
      call decomp_2d_write_plane(1,di1,output2D,-1,"./data/"//pathname//'/'//filename//'-'//num//'.bin')
    endif

  end subroutine write_field

end module visu
