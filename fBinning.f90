!! BINNING.F90
!
! BINNING: A SET OF ROUTINES FOR BINNING TOD INTO MAPS AND BASELINES.
!
!

! DOWN SAMPLING BINNING ROUTINES:

subroutine AvgCalSig(caltod,ncalsigs,cal,nsamp,nchans)
  implicit none

  integer ncalsigs,nsamp,nchans

!f2py intent(in) ncalsigs,nsamp,nchans

  real*8 caltod(nsamp,nchans)

!f2py intent(in) caltod
 
  real*8 cal(nchans)

!f2py intent(in,out) cal

  integer upper_lo , lower_lo ,steps, callen

  integer i,j, c
  
  real*8 calpower, baseline, topline

  callen   = 50
  lower_lo = 30
  upper_lo = 10
  steps = 12
  
  do c=1,nchans
     calpower = 0

     do i=0,ncalsigs-1
        baseline = 0
        topline  = 0

        do j=1,steps
           baseline = baseline + caltod(i*callen + lower_lo + j,c)/real(steps,8)
           topline  = topline  + caltod(i*callen + upper_lo + j,c)/real(steps,8)
        enddo
           
        calpower = calpower + (topline - baseline)/real(ncalsigs,8)
     
     enddo

     cal(c) = calpower

  enddo

end subroutine AvgCalSig

subroutine downsample(tod,az,el,jd,ma,outtod,outaz,outel,outjd,outma,bl,file_nb,nb,nsamp,nchan,nhorn)

  !,outtod,outaz,outel,outjd,outma
  implicit none

  !integer, intent(in) :: nsamp
  !integer, intent(in) :: nchan
  !integer, intent(in) :: nhorn
  !integer, intent(in) :: nb

  integer nsamp,nchan,nhorn,nb

!f2py intent(in) nsamp,nchan,nhorn,nb

  real*8 tod(nsamp,nchan)
  real*8 az(nsamp)
  real*8 el(nsamp)
  real*8 jd(nsamp)
  real*8 ma(nsamp,nhorn)
  
  real*8 bl
  integer file_nb


!f2py intent(in) bl, file_nb, tod, az, el ,jd, ma

  real*8 outtod(nb,nchan)
  real*8 outaz(nb)
  real*8 outel(nb)
  real*8 outjd(nb)
  real*8 outma(nb,nhorn)

!f2py intent(in,out) outtod, outaz, outel, outjd, outma

  integer i,j,c

  !Loop through the baselines
  do i=0,file_nb-1
     
     !For each baseline loop through each sample
     do j=1, bl

        !Down sample tod:
        do c=1,nchan
           outtod(i+1,c) = outtod(i+1,c) + tod(i*bl + j,c)/bl
        enddo

        !Down sample mod angle
        do c=1,nhorn
           outma(i+1,c) = outma(i+1,c) + ma(i*bl + j,c)/bl
        enddo


        outaz(i+1) = outaz(i+1) + az(i*bl + j)/bl
        outel(i+1) = outel(i+1) + el(i*bl + j)/bl
        outjd(i+1) = outjd(i+1) + jd(i*bl + j)/bl

     enddo

  enddo

end subroutine downsample


! BASELINE BINNING ROUTINES:

subroutine bin_to_baselines(map,pix,bl,cn,hitmap,limit,nsamp,pixels,nb,x)
  implicit none

  integer nsamp,pixels,nb,bl,limit
!f2py intent(in) :: nsamp, pixels, nb, bl
  real*8 map(pixels),hitmap(pixels),cn(nb)
!f2py intent(in) ::  map,hitmap,cn
  integer pix(nsamp)
!f2py intent(in) ::  pix
  real*8 x(nb)
!f2py intent(out) x


  integer i,xi,ip
  
  do i=1,nsamp
     xi = (i-1)/bl + 1
     ip = pix(i) + 1

     if (hitmap(ip) > limit) then 
        x(xi) = x(xi) + map(ip)/cn(xi)
     else
        x(xi) = 0.0
     end if
  enddo

  !do i=1,nb
  !   x(i) = x(i)/cn(i)
  !enddo

end subroutine bin_to_baselines

! MAP BINNING ROUTINES:

subroutine bin_pix(x,upix,nupix,ips,cn,sw,hw,bl,nsamp,pixels,upix_len,nb)
  !x     = input TOD (input)
  !upix  = array of unique pixel numbers (input)
  !nupix = array of hits per unique pixel (input)
  !ips   = array of TOD pixel values in ascending order (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (output)
  !hw   = weights map (output)

  implicit none

  integer, intent(in) :: nsamp
  integer, intent(in) :: pixels
  integer, intent(in) :: upix_len
  integer, intent(in) :: bl
  integer, intent(in) :: nb
!f2py intent(in) nsamp, pixels, upix_len,bl,nb
  real*8, intent(in) :: cn(nb)
  real*8, intent(in) :: x(nsamp)
!f2py real*8 x,cn
  integer, intent(in) :: ips(nsamp)
  integer, intent(in) :: upix(upix_len)
  integer, intent(in) :: nupix(upix_len)
!f2py integer upix,nupix,ips
  real*8, intent(inout) :: sw(pixels)
  real*8, intent(inout) :: hw(pixels)
!f2py real*8 sw,hw


  integer i,j,last
  
  last = 0
  do i=1, upix_len
     
     do j=1,nupix(i)
        sw(upix(i)+1) = sw(upix(i)+1) + x(ips(last + j) + 1)/cn( (ips(last+j))/bl +1)
        hw(upix(i)+1) = hw(upix(i)+1) + 1.0/cn( (ips(last+j))/bl +1 )
     end do

     last = last + nupix(i)

  end do

end subroutine bin_pix


subroutine bin_pix_hits_jackknife(x,upix,nupix,ips,cn,sw1,hw1,sw2,hw2,nmap,hits,bl,nsamp,pixels,upix_len,nb)
  !x     = input TOD (input)
  !upix  = array of unique pixel numbers (input)
  !nupix = array of hits per unique pixel (input)
  !ips   = array of TOD pixel values in ascending order (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (output)
  !hw   = weights map (output)
  !nmap = noise map (output)
  !hits = hits map (output)

  implicit none

  integer, intent(in) :: nsamp
  integer, intent(in) :: pixels
  integer, intent(in) :: upix_len
  integer, intent(in) :: bl
  integer, intent(in) :: nb
!f2py integer nsamp, pixels, upix_len,bl,nb
  real*8, intent(in) :: cn(nb)
  real*8, intent(in) :: x(nsamp)
!f2py real*8 x,cn
  integer, intent(in) :: ips(nsamp)
  integer, intent(in) :: upix(upix_len)
  integer, intent(in) :: nupix(upix_len)
!f2py integer upix,nupix,ips
  real*8, intent(inout) :: sw1(pixels)
  real*8, intent(inout) :: hw1(pixels)

  real*8, intent(inout) :: sw2(pixels)
  real*8, intent(inout) :: hw2(pixels)

  real*8, intent(inout) :: nmap(pixels)
  real*8, intent(inout) :: hits(pixels)
!f2py real*8 sw1,hw1,sw2,hw2,nmap,hits


  integer i,j,last
  real*8 xmap(pixels),mean,x2map(pixels)

  last = 0
  do i=1, upix_len
     

     xmap(upix(i)+1) = 0
     x2map(upix(i)+1) = 0

     do j=1,nupix(i)
        if (j < nupix(i)/2) then
           sw1(upix(i)+1) = sw1(upix(i)+1) + x(ips(last + j) + 1)/cn( (ips(last+j))/bl +1)
           hw1(upix(i)+1) = hw1(upix(i)+1) + 1.0/cn( (ips(last+j))/bl +1 )
        else
           sw2(upix(i)+1) = sw2(upix(i)+1) + x(ips(last + j) + 1)/cn( (ips(last+j))/bl +1)
           hw2(upix(i)+1) = hw2(upix(i)+1) + 1.0/cn( (ips(last+j))/bl +1 )
        end if

        hits(upix(i)+1) = hits(upix(i)+1) + 1
        xmap(upix(i)+1) = xmap(upix(i)+1) + x(ips(last + j) + 1)
        x2map(upix(i)+1) = x2map(upix(i)+1) + x(ips(last + j) + 1)*x(ips(last + j) + 1)

     end do

     mean = xmap(upix(i)+1)/hits(upix(i)+1)
     nmap(upix(i)+1) = (x2map(upix(i)+1)/hits(upix(i)+1) - mean*mean)

     last = last + nupix(i)

  end do

end subroutine bin_pix_hits_jackknife


subroutine bin_pix_hits(x,pix,cn,sw,hw,hits,bl,nsamp,pixels,nb)
  !x     = input TOD (input)
  !upix  = array of unique pixel numbers (input)
  !nupix = array of hits per unique pixel (input)
  !ips   = array of TOD pixel values in ascending order (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (output)
  !hw   = weights map (output)
  !nmap = noise map (output)
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

!f2py  intent(in)  x, cn , pix

  real*8 sw(pixels)
  real*8 hw(pixels)
  real*8 hits(pixels)
!f2py intent(in,out) sw,hw,hits


  integer i
 
  do i=1, nsamp
     
     sw(pix(i)+1)   = sw(pix(i)+1)   + x(i)/cn( ((i-1)/bl) +1 )
     hw(pix(i)+1)   = hw(pix(i)+1)   + 1.0 /cn( ((i-1)/bl) +1 )
     hits(pix(i)+1) = hits(pix(i)+1) + 1
     !print*,x(i)/cn( i/bl +1),1.0 /cn( i/bl +1 ),i/bl +1,nb

  end do

end subroutine bin_pix_hits


subroutine bin_pix_pol(x,a,p,cn,sw,hws,hwc,cc,ss,sc,hits,bl,nsamp,pixels,nb)
  !x     = input TOD (input)
  !upix  = array of unique pixel numbers (input)
  !nupix = array of hits per unique pixel (input)
  !ips   = array of TOD pixel values in ascending order (input)
  !cn    = array of weights for each TOD value (input)
  !bl    = baseline length (input)

  !sw   = signal*weights map (output)
  !hw   = weights map (output)
  !nmap = noise map (output)
  !hits = hits map (output)

  implicit none

  integer, intent(in) :: nsamp
  integer, intent(in) :: pixels
  integer, intent(in) :: nb

  integer, intent(in) :: bl
!f2py integer nsamp, pixels,bl,nb
  real*8, intent(in) :: cn(nb)
  real*8, intent(in) :: x(nsamp)
  real*8, intent(in) :: a(nsamp)

!f2py real*8 x,cn,a

  integer, intent(in) :: p(nsamp)
!f2py integer p

  real*8, intent(inout) :: sw(pixels)
  real*8, intent(inout) :: hwc(pixels)
  real*8, intent(inout) :: hws(pixels)

  real*8, intent(inout) :: cc(pixels)
  real*8, intent(inout) :: ss(pixels)
  real*8, intent(inout) :: sc(pixels)

  real*8, intent(inout) :: hits(pixels)
!f2py real*8 sw,hwc,hws,cc,ss,sc,hits

  integer i
  
  do i=1, nsamp

     sw(p(i)+1)  = sw(p(i)+1)  + 1.0/cn( i/bl +1 )
     hwc(p(i)+1) = hwc(p(i)+1) + x(i)/cn( i/bl +1)*cos(a(i))
     hws(p(i)+1) = hws(p(i)+1) + x(i)/cn( i/bl +1)*sin(a(i))

     cc(p(i) +1) = cc(p(i)+1) + cos(a(i))*cos(a(i))
     ss(p(i) +1) = ss(p(i)+1) + sin(a(i))*sin(a(i))
     sc(p(i) +1) = sc(p(i)+1) + sin(a(i))*cos(a(i))

     hits(p(i)+1) = hits(p(i)+1) + 1
  end do

end subroutine bin_pix_pol
