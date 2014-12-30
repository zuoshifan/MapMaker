!! BINNING.F90
!
! BINNING: A SET OF ROUTINES FOR BINNING TOD INTO MAPS AND BASELINES.
!
!

! BASELINE BINNING ROUTINES:

subroutine bin_to_baselines(m,pix,bl,cn,hitmap,limit,nsamp,pixels,nb,x)
  implicit none

  integer nsamp,pixels,nb,bl,limit
!f2py intent(in) nsamp, pixels, nb, bl
  real*8 m(pixels),hitmap(pixels),cn(nb)
!f2py intent(in) m,hitmap,cn
  integer pix(nsamp)
!f2py intent(in) pix
  real*8 x(nb)
!f2py intent(out) x


  integer i,xi,ip
  
  do i=1,nsamp
     xi = (i-1)/bl + 1
     ip = pix(i) + 1

     if (hitmap(ip) > limit) then 
        x(xi) = x(xi) + m(ip)/cn(xi)
     else
        x(xi) = 0.0
     end if
  enddo

end subroutine bin_to_baselines

! MAP BINNING ROUTINES:

subroutine bin_pix(x,pix,cn,bl,sw,hw,nsamp,pixels,nb)
  !x     = input TOD (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)

  implicit none

  integer nsamp
  integer pixels
  integer nb
  integer bl
!f2py intent(in) nsamp, pixels,bl,nb

  real*8 cn(nb)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw


  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( ((i-1)/bl) +1 )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( ((i-1)/bl) +1 )
  end do

end subroutine bin_pix

subroutine bin_pix_hits(x,pix,cn,bl,sw,hw,hits,nsamp,pixels,nb)
  !x     = input TOD (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)
  !hits = hits map (output)

  implicit none

  integer nsamp
  integer pixels
  integer nb
  integer bl
!f2py intent(in) nsamp, pixels,bl,nb

  real*8 cn(nb)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
  real*8 hits(pixels)
!f2py intent(in,out) sw,hw,hits


  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( ((i-1)/bl) +1 )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( ((i-1)/bl) +1 )
     hits(pix(i)+1) = hits(pix(i)+1) + 1
  end do

end subroutine bin_pix_hits

subroutine bin_pix_ext(x,pix,cn,sw,hw,bl,nsamp,pixels,nb)
  !x     = input TOD (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)

  implicit none

  integer nsamp
  integer pixels
  integer nb
  integer bl
!f2py intent(in) nsamp, pixels,bl,nb

  real*8 cn(nb)
  real*8 x(nb)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw

  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x( ((i-1)/bl) +1 )/cn( ((i-1)/bl) +1 )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( ((i-1)/bl) +1 )
  end do

end subroutine bin_pix_ext

subroutine bin_pix_with_ext(x,a,pix,cn,sw,hw,bl,nsamp,pixels,nb)
  !x     = input TOD (input)
  !x     = input baseline values (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)

  implicit none

  integer nsamp
  integer pixels
  integer nb
  integer bl
!f2py intent(in) nsamp, pixels,bl,nb

  real*8 cn(nb)
  real*8 a(nb)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn,a, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw

  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + (x(i) - a( ((i-1)/bl) +1 ))/cn( ((i-1)/bl) +1 )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( ((i-1)/bl) +1 )
  end do

end subroutine bin_pix_with_ext
