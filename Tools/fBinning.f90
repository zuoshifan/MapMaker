!! BINNING.F90
!
! BINNING: A SET OF ROUTINES FOR BINNING TOD INTO MAPS AND BASELINES.
!
!

! BASELINE BINNING ROUTINES:

subroutine bin_to_baselines(q,u,phi,pix,bl,cn,hitmap,limit,nsamp,pixels,nb,x)
  implicit none

  integer nsamp,pixels,bl,limit,nb
!f2py intent(in) nsamp, pixels, bl,nb
  real*8 q(pixels),u(pixels),hitmap(pixels),cn(nsamp),phi(nsamp)
!f2py intent(in) q,u,hitmap,cn,phi
  integer pix(nsamp)
!f2py intent(in) pix
  real*8 x(nb)
!f2py intent(out) x


  integer i,xi,ip
  
  do i=1,nsamp
     xi = (i-1)/bl + 1
     ip = pix(i) + 1

     x(xi) = x(xi) + ( q(ip)*sin(phi(i)) + u(ip)*cos(phi(i)) )/cn(i)/real(bl)
  enddo

end subroutine bin_to_baselines

! MAP BINNING ROUTINES:

subroutine bin_pix(x,pix,cn,bl,sw,hw,nsamp,pixels)
  !x     = input TOD (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)

  implicit none

  integer nsamp
  integer pixels
  integer bl
!f2py intent(in) nsamp, pixels,bl

  real*8 cn(nsamp)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw


  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( i )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
  end do

end subroutine bin_pix

subroutine bin_pix_hits(x,pix,cn,bl,sw,hw,hits,nsamp,pixels)
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
  integer bl
!f2py intent(in) nsamp, pixels,bl

  real*8 cn(nsamp)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
  real*8 hits(pixels)
!f2py intent(in,out) sw,hw,hits


  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( i )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
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

  real*8 cn(nsamp)
  real*8 x(nb)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw

  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + x( ((i-1)/bl) +1 )/cn( i )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
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

  real*8 cn(nsamp)
  real*8 a(nb)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn,a, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw

  integer i
 
  do i=1, nsamp
     sw(pix(i)+1)   = sw(pix(i)+1)   + (x(i) - a( ((i-1)/bl) +1 ))/cn( i )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
  end do

end subroutine bin_pix_with_ext


subroutine DownSample(a,newlen,alen,out,err)
  implicit none

  integer alen, newlen
!f2py intent(in) alen, newlen
  
  real*8 a(alen)
!f2py intent(in) a

  real*8 out(newlen),err(newlen)
!f2py intent(out) out, err


  integer nsteps,i, k , count

  real*8 suma,suma2,fnsteps,laststep

  nsteps = alen/newlen
  fnsteps = real(nsteps)

  suma = 0
  suma2= 0
  k = 1
  laststep = fnsteps + real(alen - newlen*nsteps)
  count = 0

  do i=1,alen

     if (i .eq. alen-laststep+1) then
        nsteps  = int(laststep)
        fnsteps = laststep
     end if

     !Calc mean and stderr of samples in step
     suma  = suma  + a(i)
     suma2 = suma2 + a(i)*a(i)
     
     if (mod(i-count,nsteps) .eq. 0) then             
        out(k) = suma /fnsteps
        err(k) = suma2/fnsteps - suma*suma/fnsteps/fnsteps
        !If last err only has one element
        ! Use error of previous box.
        if (err(k) .eq. 0 .and. k .gt. 1) then
           err(k) = err(k-1)
        end if
        
        !reset suma and suma2
        suma = 0
        suma2= 0
        
        !increase output index k
        k = k + 1
        count = i
     end if

        
  end do

end subroutine DownSample
