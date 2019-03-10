subroutine filter(af)

  USE param
  USE derivX
  USE derivY
  USE derivZ
  USE variables
  USE var
!=================================================
! Discrete low-pass filter according to 
!=================================================	
  implicit none
  real(mytype),intent(in) :: af 

#ifdef DEBG
  if (nrank .eq. 0) print *,'# filter calculation start'
#endif

  ! Filter functions
  if (nclx1.eq.0.and.nclxn.eq.0) filx => filx_00
  if (nclx1.eq.1.and.nclxn.eq.1) filx => filx_11
  if (nclx1.eq.1.and.nclxn.eq.2) filx => filx_12
  if (nclx1.eq.2.and.nclxn.eq.1) filx => filx_21
  if (nclx1.eq.2.and.nclxn.eq.2) filx => filx_22
  !
  if (ncly1.eq.0.and.nclyn.eq.0) fily => fily_00
  if (ncly1.eq.1.and.nclyn.eq.1) fily => fily_11
  if (ncly1.eq.1.and.nclyn.eq.2) fily => fily_12
  if (ncly1.eq.2.and.nclyn.eq.1) fily => fily_21
  if (ncly1.eq.2.and.nclyn.eq.2) fily => fily_22
  !
  if (nclz1.eq.0.and.nclzn.eq.0) filz => filz_00
  if (nclz1.eq.1.and.nclzn.eq.1) filz => filz_11
  if (nclz1.eq.1.and.nclzn.eq.2) filz => filz_12
  if (nclz1.eq.2.and.nclzn.eq.1) filz => filz_21
  if (nclz1.eq.2.and.nclzn.eq.2) filz => filz_22

return 

end subroutine filter


subroutine set_filter_coefficients(af,alfa1,a1,b1,c1,d1,alfa2,a2,b2,c2,d2,alfa3,a3,b3,c3,d3,e3,f3,&
                                   alfan,an,bn,cn,dn,alfam,am,bm,cm,dm,alfap,ap,bp,cp,dp,ep,fp,&
                                   alfai,ai,bi,ci,di,ff,fs,fw,ffp,fsp,fwp,n,ncl1,ncln)

  use decomp_2d, only : mytype, nrank
  use param

  implicit none

  real(mytype),intent(in) :: af
  integer,intent(in) :: n,ncl1,ncln
  real(mytype),dimension(n),intent(out) :: ff,fs,fw,ffp,fsp,fwp
  real(mytype),intent(out) :: alfa1,a1,b1,c1,d1,alfa2,a2,b2,c2,d2,alfa3,a3,b3,c3,d3,e3,f3,&
                                   alfan,an,bn,cn,dn,alfam,am,bm,cm,dm,alfap,ap,bp,cp,dp,ep,fp,&
                                   alfai,ai,bi,ci,di
  integer :: i
  real(mytype),dimension(n) :: fb,fc

    ! Set the coefficient for the discrete filter following 
    ! the tridiagonal filtering of Motheau and Abraham, JCP 2016 
    ! Filter should be -0.5<filax<0.5
    
    ! General Case (entire points)
    ! alpha*fhat(i-1)+fhat(i)+alpha*fhat(i+1)=af(i)+b/2*[f(i+1)+f(i-1)] + ...
    
    ! Coefficients are calculated according to the report of Gaitonde & Visbal, 1998,
    ! "High-order schemes for Navier-Stokes equations: Algorithm and implementation into FDL3DI"

  
    alfai=af                          ! alpha_f
    !Interior points
    ai=(eleven + ten*af)/sixteen     ! a
    bi=half*(af +thirtyfour*af)/32.  ! b/2 
    ci=half*(-three + six*af)/16.    ! c/2
    di=half*(one - two*af)/thirtytwo ! d/2
    ! Explicit third/fifth-order filters near the boundaries!
    !Boundary point 1 (no-filtering)
    alfa1=zero
    a1=one                           ! a1=7./8.+af/8.! a1/2
    b1=zero                          ! b1=3./8.+5.*af/8.  
    c1=zero                          ! c1=-3./8.+3./8.*af 
    d1=zero                          ! d1=1./8.-1./8.*af 
    !Boundary point 2 (Third order)
    alfa2=af
    a2=one/eight+three/four*af          ! a2
    b2=five/eight+three/four*af          !  b2
    c2=three/eight+af/four               ! c2
    d2=-one/eight+af/four                ! d2
    !Boundary point 3 (Fifth order)
    alfa3=af
    a3= -one/thirtytwo+af/sixteen         ! a3
    b3= five/thirtytwo+eleven/sixteen*af  ! b3
    c3= eleven/sixteen+five*af/eight      ! c3
    d3= five/sixteen+three*af/eight       ! d3
    e3=-five/thirtytwo+five*af/sixteen    ! e3
    f3= one/thirtytwo-af/sixteen          ! f3
    !Boundary point n (no-filtering)
    alfan=zero
    an=one                                !an = 7./8.+af/8.! a1/2
    bn=zero                               !bn = 3./8.+5.*af/8.
    cn=zero                               !cn =-3./8.+3./8.*af    
    dn=zero                               !dn = 1./8.-1./8.*af    
    !Boundary point 2 (Third order)
    alfam=af
    am=one/eight+three/four*af          ! a2
    bm=five/eight+three/four*af          !  b2
    cm=three/eight+af/four               ! c2
    dm=-one/eight+af/four                ! d2
    !Boundary point 3 (Fifth order)
    alfap=af
    ap=-one/thirtytwo+af/sixteen         ! a3
    bp= five/thirtytwo+eleven/sixteen*af  ! b3
    cp= eleven/sixteen+five*af/eight      ! c3
    dp= five/sixteen+three*af/eight       ! d3
    ep=-five/thirtytwo+five*af/sixteen    ! e3
    fp= one/thirtytwo-af/sixteen          ! f3

    ff=zero;fs=zero;fw=zero;ffp=zero;fsp=zero;fwp=zero
    fb=zero;fc=zero
  
  if     (ncl1.eq.0) then !Periodic
     ff(1)   =alfai
     ff(2)   =alfai
     fc(1)   =two
     fc(2)   =one
     fb(1)   =alfai
     fb(2)   =alfai
  elseif (ncl1.eq.1) then !Free-slip
     ff(1)   =alfai+alfai
     ff(2)   =alfai
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfai 
     fb(2)   =alfai
  elseif (ncl1.eq.2) then !Dirichlet
     ff(1)   =alfa1
     ff(2)   =alfa2
     fc(1)   =one
     fc(2)   =one
     fb(1)   =alfa2 
     fb(2)   =alfai
  endif
  if (ncln.eq.0) then !Periodic
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one+alfai*alfai
     fb(n-2)=alfai
     fb(n-1)=alfai
     fb(n  )=zero
  elseif (ncln.eq.1) then !Free-slip
     ff(n-2)=alfai
     ff(n-1)=alfai
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one
     fb(n-2)=alfai
     fb(n-1)=alfai+alfai
     fb(n  )=zero
  elseif (ncln.eq.2) then !Dirichlet
     ff(n-2)=alfai
     ff(n-1)=alfam
     ff(n)  =zero
     fc(n-2)=one
     fc(n-1)=one
     fc(n  )=one
     fb(n-2)=alfam
     fb(n-1)=alfan
     fb(n  )=zero
  endif
  do i=3,n-3
     ff(i)=alfai
     fc(i)=one
     fb(i)=alfai
  enddo
  
  do i=1,n
     ffp(i)=ff(i)
  enddo

  call prepare (fb,fc,ff ,fs ,fw ,n)

  if (ncl1.eq.1) then
     ffp(1)=zero
  endif
  if (ncln.eq.1) then
     fb(n-1)=zero
  endif

  call prepare (fb,fc,ffp,fsp,fwp,n)

  return

end subroutine set_filter_coefficients

subroutine filx_00(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 

USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx
  
if(iibm.eq.2) call lagpolx(ux)

   do k=1,nz 
   do j=1,ny 
      tx(1,j,k)=fiaix*ux(1,j,k)+fibix*(ux(2,j,k)+ux(nx,j,k))& 
                               +ficix*(ux(3,j,k)+ux(nx-1,j,k))&
                               +fidix*(ux(4,j,k)+ux(nx-2,j,k)) 
      rx(1,j,k)=-1.
      tx(2,j,k)=fiaix*ux(2,j,k)+fibix*(ux(3,j,k)+ux(1,j,k))&
                               +ficix*(ux(4,j,k)+ux(nx,j,k))& 
                               +fidix*(ux(5,j,k)+ux(nx-1,j,k)) 
      rx(2,j,k)=0. 
      tx(3,j,k)=fiaix*ux(3,j,k)+fibix*(ux(4,j,k)+ux(2,j,k))&
          
                               +ficix*(ux(5,j,k)+ux(1,j,k))& 
                               +fidix*(ux(6,j,k)+ux(nx,j,k)) 
      rx(3,j,k)=0. 
      do i=4,nx-3
         tx(i,j,k)=fiaix*ux(i,j,k)+fibix*(ux(i+1,j,k)+ux(i-1,j,k))& 
                                  +ficix*(ux(i+2,j,k)+ux(i-2,j,k))&
                                  +fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
         rx(i,j,k)=0. 
      enddo
      tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+fibix*(ux(nx-3,j,k)+ux(nx-1,j,k))&
                                     +ficix*(ux(nx-4,j,k)+ux(nx,j,k))& 
                                     +fidix*(ux(nx-5,j,k)+ux(1,j,k)) 
      rx(nx-2,j,k)=0. 
      tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+fibix*(ux(nx-2,j,k)+ux(nx,j,k))&
                                     +ficix*(ux(nx-3,j,k)+ux(1,j,k))& 
                                     +fidix*(ux(nx-4,j,k)+ux(2,j,k)) 
      rx(nx-1,j,k)=0. 
      tx(nx,j,k)=fiaix*ux(nx,j,k)+fibix*(ux(nx-1,j,k)+ux(1,j,k))&
                                 +ficix*(ux(nx-2,j,k)+ux(2,j,k))& 
                                 +fidix*(ux(nx-3,j,k)+ux(3,j,k)) 
      rx(nx,j,k)=fialx           
      do i=2, nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fifsx(i) 
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fifsx(i) 
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fifwx(nx) 
      rx(nx,j,k)=rx(nx,j,k)*fifwx(nx) 
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-fiffx(i)*tx(i+1,j,k))*fifwx(i) 
         rx(i,j,k)=(rx(i,j,k)-fiffx(i)*rx(i+1,j,k))*fifwx(i) 
      enddo
        fisx(j,k)=(tx(1,j,k)-fialx*tx(nx,j,k))&
           /(1.+rx(1,j,k)-fialx*rx(nx,j,k)) 
      do i=1,nx 
         tx(i,j,k)=tx(i,j,k)-fisx(j,k)*rx(i,j,k) 
      enddo
   enddo
   enddo

return  

end subroutine filx_00

subroutine filx_11(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 

USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

   if(iibm.eq.2) call lagpolx(ux)

   print *, 'Not ready yet'
   stop

   return

end subroutine filx_11

subroutine filx_12(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 

USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

   if(iibm.eq.2) call lagpolx(ux)

   print *, 'Not ready yet'
   stop

   return

end subroutine filx_12

subroutine filx_21(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 

USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

   if(iibm.eq.2) call lagpolx(ux)

   print *, 'Not ready yet'
   stop

   return

end subroutine filx_21


subroutine filx_22(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 
  
USE param  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

   if(iibm.eq.2) call lagpolx(ux)

   do k=1,nz
   do j=1,ny 
      tx(1,j,k)=ux(1,j,k)
      tx(2,j,k)=fia2x*ux(1,j,k)+fib2x*ux(2,j,k)+fic2x*ux(3,j,k)+&
                fid2x*ux(4,j,k)
      tx(3,j,k)=fia3x*ux(1,j,k)+fib3x*ux(2,j,k)+fic3x*ux(3,j,k)+&
                fid3x*ux(4,j,k)+fie3x*ux(5,j,k)+fif3x*ux(6,j,k)
      do i=4,nx-3
         tx(i,j,k)=fiaix*ux(i,j,k)+fibix*(ux(i+1,j,k)+ux(i-1,j,k))& 
                                  +ficix*(ux(i+2,j,k)+ux(i-2,j,k))&
                                  +fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
      enddo
      tx(nx,j,k)=ux(nx,j,k)
      tx(nx-1,j,k)=fiamx*ux(nx,j,k)+fibmx*ux(nx-1,j,k)+ficmx*ux(nx-2,j,k)+&
                            fidmx*ux(nx-3,j,k)
      tx(nx-2,j,k)=fiapx*ux(nx,j,k)+fibpx*ux(nx-1,j,k)+ficpx*ux(nx-2,j,k)+&
                fidpx*ux(nx-3,j,k)+fiepx*ux(nx-4,j,k)+fifpx*ux(nx-5,j,k)
      do i=2,nx 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fifsx(i) 
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fifwx(nx) 
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-fiffx(i)*tx(i+1,j,k))*fifwx(i) 
      enddo
   enddo 
   enddo 

return  
end subroutine filx_22

subroutine fily_00(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy
    
if(iibm.eq.2) call lagpoly(uy)

    do k=1,nz 
    do i=1,nx 
       ty(i,1,k)=fiajy*uy(i,1,k)+fibjy*(uy(i,2,k)+uy(i,nx,k))& 
                                +ficjy*(uy(i,3,k)+uy(i,nx-1,k))&
                                +fidjy*(uy(i,4,k)+uy(i,nx-2,k)) 
       ry(i,1,k)=-1.
       ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))&
                                +ficjy*(uy(i,4,k)+uy(i,nx,k))& 
                                +fidjy*(uy(i,5,k)+uy(i,nx-1,k)) 
       ry(i,2,k)=0. 
       ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))&
                                +ficjy*(uy(i,5,k)+uy(i,1,k))& 
                                +fidjy*(uy(i,6,k)+uy(i,nx,k)) 
       ry(i,3,k)=0. 
       do j=4,ny-3
          ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                   +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                   +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
          ry(i,j,k)=0. 
       enddo
       ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-3,k)+uy(i,ny-1,k))&
                                      +ficjy*(uy(i,ny-4,k)+uy(i,ny,k))& 
                                      +fidjy*(uy(i,ny-5,k)+uy(i,1,k)) 
       ry(i,ny-2,k)=0. 
       ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny-2,k)+uy(i,ny,k))&
                                      +ficjy*(uy(i,ny-3,k)+uy(i,1,k))& 
                                      +fidjy*(uy(i,ny-4,k)+uy(i,2,k)) 
       ry(i,ny-1,k)=0. 
       ty(i,ny,k)=fiajy*uy(i,ny,k)+fibjy*(uy(i,ny-1,k)+uy(i,1,k))&
                                  +ficjy*(uy(i,ny-2,k)+uy(i,2,k))& 
                                  +fidjy*(uy(i,ny-3,k)+uy(i,3,k)) 
       ry(i,ny,k)=fialy           
       do j=2, ny
          ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
          ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*fifsy(j) 
       enddo
       ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
       ry(i,ny,k)=ry(i,ny,k)*fifwy(ny) 
       do j=ny-1,1,-1
          ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
          ry(i,j,k)=(ry(i,j,k)-fiffy(j)*ry(i,j+1,k))*fifwy(j) 
       enddo
         fisy(i,k)=(ty(i,1,k)-fialy*ty(i,ny,k))&
            /(1.+ry(i,1,k)-fialy*ry(i,ny,k)) 
       do j=1,ny 
          ty(i,j,k)=ty(i,j,k)-fisy(i,k)*ry(i,j,k) 
       enddo
    enddo
    enddo

    if (istret.ne.0) then   
        do k=1,nz 
           do j=1,ny 
              do i=1,nx 
                 ty(i,j,k)=ty(i,j,k)*ppy(j) 
              enddo
           enddo
        enddo
    endif
    
return

end subroutine fily_00

!********************************************************************
!
subroutine fily_11(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy
    
if(iibm.eq.2) call lagpoly(uy)

    if (npaire==1) then 
    do k=1,nz 
    do i=1,nx 
         ty(i,1,k)=fiajy*uy(i,1,k)+fibjy*(uy(i,2,k)+uy(i,2,k))&
                                  +ficjy*(uy(i,3,k)+uy(i,3,k))&
                                  +fidjy*(uy(i,4,k)+uy(i,4,k))
         ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))& 
                                  +ficjy*(uy(i,4,k)+uy(i,2,k))&
                                  +fidjy*(uy(i,5,k)+uy(i,3,k)) 
         ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))& 
                                  +ficjy*(uy(i,5,k)+uy(i,1,k))&
                                  +fidjy*(uy(i,6,k)+uy(i,2,k)) 
    enddo
    enddo 
    do k=1,nz 
    do j=4,ny-3 
    do i=1,nx 
       ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
    enddo
    enddo 
    enddo 
    do k=1,nz 
    do i=1,nx 
       ty(i,ny,k)=fiajy*uy(i,ny,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-1,k))&
                                  +ficjy*(uy(i,ny-2,k)+uy(i,ny-2,k))&
                                  +fidjy*(uy(i,ny-3,k)+uy(i,ny-3,k))
       ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)  +uy(i,ny-2,k))& 
                                      +ficjy*(uy(i,ny-1,k)+uy(i,ny-3,k))&
                                      +fidjy*(uy(i,ny-2,k)+uy(i,ny-4,k)) 
       ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k))& 
                                      +ficjy*(uy(i,ny,k)+uy(i,ny-4,k))&
                                      +fidjy*(uy(i,ny-1,k)+uy(i,ny-5,k)) 
    enddo
    enddo 
    do k=1,nz
    do j=2,ny  
    do i=1,nx 
       ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
    enddo
    enddo
    enddo 
    do k=1,nz 
    do i=1,nx 
       ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
    enddo 
    enddo 
    do k=1,nz
    do j=ny-1,1,-1  
    do i=1,nx 
    ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
    enddo 
    enddo 
    enddo 
   endif
   if (npaire==0) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,1,k)=fiajy*uy(i,1,k)
         ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))& 
                                  +ficjy*(uy(i,4,k)-uy(i,2,k))&
                                  +fidjy*(uy(i,5,k)-uy(i,3,k)) 
         ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))& 
                                  +ficjy*(uy(i,5,k)+uy(i,1,k))&
                                  +fidjy*(uy(i,6,k)-uy(i,2,k)) 
      enddo
      enddo 
      do k=1,nz 
      do j=4,ny-3 
      do i=1,nx 
         ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                  +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                  +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=fiajy*uy(i,ny,k)
         ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)+uy(i,ny-2,k))& 
                                        +ficjy*(-uy(i,ny-1,k)+uy(i,ny-3,k))&
                                        +fidjy*(-uy(i,ny-2,k)+uy(i,ny-4,k)) 
         ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k))& 
                                        +ficjy*(uy(i,ny,k)+uy(i,ny-4,k))&
                                        +fidjy*(-uy(i,ny-1,k)+uy(i,ny-5,k)) 
      enddo
      enddo 
      do k=1,nz
      do j=2,ny  
      do i=1,nx 
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
      enddo
      enddo
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
      enddo 
      enddo 
      do k=1,nz
      do j=ny-1,1,-1  
      do i=1,nx 
      ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
      enddo 
      enddo 
      enddo 
   endif

    if (istret.ne.0) then   
        do k=1,nz 
           do j=1,ny 
              do i=1,nx 
                 ty(i,j,k)=ty(i,j,k)*ppy(j) 
              enddo
           enddo
        enddo
    endif
    
return

end subroutine fily_11


subroutine fily_12(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy
    
   if(iibm.eq.2) call lagpoly(uy)
   print *, 'Not ready yet'
   stop

end subroutine fily_12


subroutine fily_21(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy
    
   if(iibm.eq.2) call lagpoly(uy)
   print *, 'Not ready yet'
   stop

end subroutine fily_21

subroutine fily_22(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
  
USE param  
USE parfiY 

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy
    
   if(iibm.eq.2) call lagpoly(uy)
   print *, 'Not ready yet'
   stop

end subroutine fily_22

subroutine filz_00(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

  if(iibm.eq.2) call lagpolz(uz)

   do j=1,ny 
   do i=1,nx 
      tz(i,j,1)=fiakz*uz(i,j,1)+fibkz*(uz(i,j,2)+uz(i,j,nz))& 
                               +fickz*(uz(i,j,3)+uz(i,j,nz-1))&
                               +fidkz*(uz(i,j,4)+uz(i,j,nz-2)) 
      rz(i,j,1)=-1.
      tz(i,j,2)=fiakz*uz(i,j,2)+fibkz*(uz(i,j,3)+uz(i,j,1))&
                               +fickz*(uz(i,j,4)+uz(i,j,nz))& 
                               +fidkz*(uz(i,j,5)+uz(i,j,nz-1)) 
      rz(i,j,2)=0. 
      tz(i,j,3)=fiakz*uz(i,j,3)+fibkz*(uz(i,j,4)+uz(i,j,2))&
                               +fickz*(uz(i,j,5)+uz(i,j,1))& 
                               +fidkz*(uz(i,j,6)+uz(i,j,nz)) 
      rz(i,j,3)=0.
   enddo
   enddo 
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
         tz(i,j,k)=fiakz*uz(i,j,k)+fibkz*(uz(i,j,k+1)+uz(i,j,k-1))& 
                                  +fickz*(uz(i,j,k+2)+uz(i,j,k-2))&
                                  +fidkz*(uz(i,j,k+3)+uz(i,j,k-3)) 
         rz(i,j,k)=0. 
   enddo
   enddo
   enddo 
   do j=1,ny 
   do i=1,nx 
      tz(i,j,nz-2)=fiakz*uz(i,j,nz-2)+fibkz*(uz(i,j,nz-3)+uz(i,j,nz-1))&
                                     +fickz*(uz(i,j,nz-4)+uz(i,j,nz))& 
                                     +fidkz*(uz(i,j,nz-5)+uz(i,j,1)) 
      rz(i,j,nz-2)=0. 
      tz(i,j,nz-1)=fiakz*uz(i,j,nz-1)+fibkz*(uz(i,j,nz-2)+uz(i,j,nz))&
                                     +fickz*(uz(i,j,nz-3)+uz(i,j,1))& 
                                     +fidkz*(uz(i,j,nz-4)+uz(i,j,2)) 
      rz(i,j,nz-1)=0. 
      tz(i,j,nz)=fiakz*uz(i,j,nz)+fibkz*(uz(i,j,nz-1)+uz(i,j,1))&
                                 +fickz*(uz(i,j,nz-2)+uz(i,j,2))& 
                                 +fidkz*(uz(i,j,nz-3)+uz(i,j,3)) 
      rz(i,j,nz)=fialz           
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fifsz(k) 
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fifsz(k) 
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx 
      tz(i,j,nz)=tz(i,j,nz)*fifwz(nz) 
      rz(i,j,nz)=rz(i,j,nz)*fifwz(nz) 
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx    
         tz(i,j,k)=(tz(i,j,k)-fiffz(k)*tz(i,j,k+1))*fifwz(k) 
         rz(i,j,k)=(rz(i,j,k)-fiffz(k)*rz(i,j,k+1))*fifwz(k) 
   enddo
   enddo
   enddo   
   do j=1,ny
   do i=1,nx      
    fisz(i,j)=(tz(i,j,1)-fialz*tz(i,j,nz))&
           /(1.+rz(i,j,1)-fialz*rz(i,j,nz)) 
   enddo
   enddo    
   do k=1,nz 
   do j=1,ny
   do i=1,nx      
         tz(i,j,k)=tz(i,j,k)-fisz(i,j)*rz(i,j,k) 
   enddo
   enddo
   enddo

return  
end subroutine filz_00

subroutine filz_11(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

  if(iibm.eq.2) call lagpolz(uz)
   
  print *, 'Not ready yet'
  stop

end subroutine filz_11

subroutine filz_12(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

  if(iibm.eq.2) call lagpolz(uz)
   
  print *, 'Not ready yet'
  stop

end subroutine filz_12

subroutine filz_21(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

  if(iibm.eq.2) call lagpolz(uz)
   
  print *, 'Not ready yet'
  stop

end subroutine filz_21


subroutine filz_22(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

  if(iibm.eq.2) call lagpolz(uz)
   
  print *, 'Not ready yet'
  stop

end subroutine filz_22

