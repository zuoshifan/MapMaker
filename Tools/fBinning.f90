SUBROUTINE sort(a,n)
  INTEGER n
  REAL*8 a(n)
  REAL*8 temp
  INTEGER :: i, j
  LOGICAL :: swapped = .TRUE.

  DO j = SIZE(a)-1, 1, -1
    swapped = .FALSE.
    DO i = 1, jx1
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE sort

!! BINNING.F90
!
! BINNING: A SET OF ROUTINES FOR BINNING TOD INTO MAPS AND BASELINES.
!
!

! BASELINE BINNING ROUTINES:

subroutine bin_to_baselines(m,pix,bl,cn,nsamp,pixels,nb,x)
  implicit none

  integer nsamp,pixels,nb,bl
!f2py intent(in) nsamp, pixels, nb, bl
  real*8 m(pixels),cn(nsamp)
!f2py intent(in) m,cn
  integer pix(nsamp)
!f2py intent(in) pix
  real*8 x(nb)
!f2py intent(out) x


  integer i,xi,ip
  
  do i=1,nsamp
     xi = (i-1)/bl + 1
     ip = pix(i) + 1

     if (cn(i) < 1e20) then 
        x(xi) = x(xi) + m(ip)/cn(i)
     endif

  enddo

end subroutine bin_to_baselines

subroutine bin_to_baselines_pmask(m,pix,bl,pmask,cn,nsamp,pixels,nb,x)
  implicit none

  integer nsamp,pixels,nb,bl
!f2py intent(in) nsamp, pixels, nb, bl
  real*8 m(pixels),cn(nsamp)
  integer pmask(nsamp)
!f2py intent(in) m,cn,pmask
  integer pix(nsamp)
!f2py intent(in) pix
  real*8 x(nb)
!f2py intent(out) x


  integer i,xi,ip
  
  do i=1,nsamp
     xi = (i-1)/bl + 1
     ip = pix(i) + 1

     if (cn(i) < 1e20) then 
        if (pmask(i) .eq. 1) then
           x(xi) = x(xi) + m(ip)/cn(i)
        endif
     endif

  enddo

end subroutine bin_to_baselines_pmask

subroutine bin_ft(vals,cn,bl,nb,nsamp,x)
  implicit none

  integer nb,nsamp,bl
!f2py intent(in) nb,nsamp,bl
  real*8 vals(nsamp),cn(nsamp)
!f2py intent(in) vals,cn
  real*8 x(nb)
!f2py intent(out) x
  
  integer i,xi

  do i=1,nsamp
     xi = (i-1)/bl + 1

     if (cn(i) < 1e20) then 
        x(xi) = x(xi) + vals(i)/cn(i)
     endif

  enddo

end subroutine bin_ft

subroutine bin_ft_pmask(vals,cn,bl,nb,pmask,nsamp,x)
  implicit none

  integer nb,nsamp,bl
!f2py intent(in) nb,nsamp,bl
  real*8 vals(nsamp),cn(nsamp)
  integer pmask(nsamp)
!f2py intent(in) vals,cn,pmask
  real*8 x(nb)
!f2py intent(out) x
  
  integer i,xi

  do i=1,nsamp
     xi = (i-1)/bl + 1

     if (cn(i) < 1e20) then 
        if (pmask(i) .eq. 1) then 
           x(xi) = x(xi) + vals(i)/cn(i)
        endif
     endif

  enddo

end subroutine bin_ft_pmask


subroutine bin_ft_ext(a,cn,bl,nb,nsamp,x)
  implicit none

  integer nb,nsamp,bl
!f2py intent(in) nb,nsamp,bl
  real*8 a(nb),cn(nsamp)
!f2py intent(in) a,cn
  real*8 x(nb)
!f2py intent(out) x
  
  integer i,xi

  do i=1,nsamp
     xi = (i-1)/bl + 1

     if (cn(i) < 1e20) then 
        x(xi) = x(xi) + a(xi)/cn(i)
     endif

  enddo

end subroutine bin_ft_ext

subroutine bin_ft_ext_pmask(a,cn,bl,pmask,nb,nsamp,x)
  implicit none

  integer nb,nsamp,bl
!f2py intent(in) nb,nsamp,bl
  real*8 a(nb),cn(nsamp)
  integer pmask(nsamp)
!f2py intent(in) a,cn,pmask
  real*8 x(nb)
!f2py intent(out) x
  
  integer i,xi
  real*8 atemp

  print *, 'HELLO'

  atemp = 0.
  do i=1,nsamp

     if (cn(i) < 1e20) then 

        if (pmask(i) .eq. 1) then 
           xi = ((i-1)/bl) +1 
           atemp = a( xi )
           x(xi) = x(xi) + atemp/cn(i)
        endif
     endif

  enddo

end subroutine bin_ft_ext_pmask


! MAP BINNING ROUTINES:

subroutine bin_pix(x,pix,cn,sw,hw,nsamp,pixels)
  !x     = input TOD (input)
  !pix   = array of pixel numbers (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (in/output)
  !hw   = weights map (in/output)

  implicit none

  integer nsamp 
  integer pixels
!f2py intent(in) nsamp, pixels

  real*8 cn(nsamp)
  real*8 x(nsamp)
  integer pix(nsamp)
!f2py  intent(in)  cn, x ,pix

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw


  integer i
 
  do i=1, nsamp

     if (cn(i) < 1e20) then 
        sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( i )
        hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
     endif
  end do

end subroutine bin_pix

subroutine bin_pix_hits(x,pix,cn,sw,hw,hits,nsamp,pixels)
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
!f2py intent(in) nsamp, pixels

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
     if (cn(i) < 1e20) then 
        sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn(i)
        hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn(i)
        hits(pix(i)+1) = hits(pix(i)+1) + 1
     endif

  end do

end subroutine bin_pix_hits

subroutine bin_pix_hits_pmask(x,pix,cn,sw,hw,hits,pmask,nsamp,pixels)
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
!f2py intent(in) nsamp, pixels

  real*8 cn(nsamp)
  real*8 x(nsamp)
  integer pix(nsamp)
  integer pmask(nsamp)
!f2py  intent(in)  cn, x ,pix,pmask

  real*8 sw(pixels)
  real*8 hw(pixels)
  real*8 hits(pixels)
!f2py intent(in,out) sw,hw,hits


  integer i
 
  do i=1, nsamp
     if (cn(i) < 1e20) then 
        if (pmask(i) .eq. 1) then
           sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn(i)
           hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn(i)
           hits(pix(i)+1) = hits(pix(i)+1) + 1
        endif
     endif

  end do

end subroutine bin_pix_hits_pmask


subroutine median(a,med,nsamp)
  !USE m_mrgrnk
  implicit none

  integer  nsamp
!f2py intent(in) nsamp

  real*8 ,intent(in) :: a(nsamp)
!f2py intent(in) a

  real*8 ,intent(out) :: med
!f2py intent(out) median

  real*8 ac(nsamp)
  integer indices(nsamp)

  ac = a
  call sort(ac,size(ac))

  if (nsamp > 0) then
     if ( mod(nsamp,2) == 0) then
        med = (ac(nsamp/2 + 1) + ac(nsamp/2))/2.0
     else
        med = ac(nsamp/2 + 1)
     end if
  else
     med = 1e-24
  end if
end subroutine median


subroutine bin_pix_hits_median(x,pix,cn,sw,hw,hits,nsamp,pixels)
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
!f2py intent(in) nsamp, pixels

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
     if (cn(i) < 1e20) then 
        sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn(i)
        hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn(i)
        hits(pix(i)+1) = hits(pix(i)+1) + 1
     endif

  end do

end subroutine bin_pix_hits_median


subroutine bin_pix_ext_pmask(x,pix,cn,sw,hw,pmask,bl,nsamp,pixels,nb)
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
  integer pmask(nsamp)
!f2py  intent(in)  cn, x ,pix,pmask

  real*8 sw(pixels)
  real*8 hw(pixels)
!f2py intent(in,out) sw,hw

  integer i,xi
  real*8 xtemp
 

  !print *, 'Hello'
  print *,'HELLO'

  xtemp = 0.
  do i=1, nsamp
     if (cn(i) < 1e20) then 

        if (pmask(i) .eq. 1) then 
           xi = ((i-1)/bl) + 1 
           xtemp = x( xi )
           sw(pix(i)+1)   = sw(pix(i)+1)   + xtemp/cn( i )
           hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
        endif

     endif
  end do

end subroutine bin_pix_ext_pmask


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
     if (cn(i) < 1e20) then 
        sw(pix(i)+1)   = sw(pix(i)+1)   + x( ((i-1)/bl) +1 )/cn( i )
        hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( i )
     endif
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
