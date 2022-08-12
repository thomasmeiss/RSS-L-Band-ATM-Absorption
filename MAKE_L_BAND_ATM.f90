! Thomas Meissner
! Remote Sensing Systems
! RSS L-band ATMOSPHERIC Absorption 
! 08/11/2022
      
! MAKE_L_BAND_ATM
! Create ancillary fileds for AO, AV, AL, TRAN, TBUP, TBDW
! at SMAP frequency and relevant Earth Incidence (EIA) range
! ingest NCEP GFS/GDAS 0.25 deg surface and atmospheric profiles 
! for temperature, height, relative humidity, cloud water mixing ratio.
! write out as raw binary without header
! no space or time interpolation.
!      
! RSS-Compile: Menu -> Tools -> compile64

! internal version history
! Change ilat and ilon to go from 1-721 and 1-1440, respectively - AManaster 11/2021
! v2: change back from Wentz-Meissner 2016 to Liebe O2 absorption


include 'fdabscoeff_v2.f90'
include 'dielectric.f90'
include 'fdcldabs.f90'
include 'goff_gratch_vap.f90'
include 'column.f90'
include 'get_atm.f90'
include 'findncep_025deg.f90'
 


    program MAKE_L_BAND_ATM 
	implicit none

	real(4), parameter              :: freq=1.413  ! center frequency
	integer(4), parameter           :: icld=1      ! 1 = include coukld liquid water absorption

    character(len=200), parameter   :: outpath  = 'sample_data\'  ! output directory. specified by user.
    character(len=200)              :: outfile, str 
    
	real(4)                         :: xtbup(39:41), xtbdw(39:41), xtran(39:41)
	real(4)                         :: xao, xav, xal
	real(4)                         :: tbup(1:1440,1:721,39:41), tbdw(1:1440,1:721,39:41), tran(1:1440,1:721,39:41) 
	real(4)                         :: ao(1:1440,1:721), av(1:1440,1:721), al(1:1440,1:721)
	
	integer(4)                      :: lyear, imon, iday, ihour
	integer(4)                      :: ilat, ilon


	write(*,*) ' L-BAND ATM for SMAP'
	write(*,*) ' create ATM ABSORPTION, TRAN, TBUP, TBDW form 0.25 DEG NCEP ATMOSPHERIC PROFILES and SURFACE.'
	write(*,*) ' output on regular 0.25 DEG NCEP GRID.'
	write(*,*) 
	
	
	write(*,*)' enter year'
	read(*,*) lyear
	
	write(*,*)' enter month'
	read(*,*) imon

	write(*,*)' enter day of month'
	read(*,*) iday
	
	write(*,*)'enter hour of day of analysis time: 00, 06, 12, 18'
	read(*,*) ihour
	
	write(*,*)
	write(*,*) ' start processing'	

	write(*,*) ' time stamp ',lyear,imon,iday,ihour

	do ilat=1,721
	write(*,*) ' lat ',ilat
    
    do ilon=1,1440
   
    call  get_atm(icld,lyear,imon,iday,ihour,ilat,ilon,freq, xtbup,xtbdw,xtran,xao,xav,xal) 

    ao(ilon,ilat)=xao
    av(ilon,ilat)=xav
    al(ilon,ilat)=xal
 
    tran(ilon,ilat,:)=xtran
    tbup(ilon,ilat,:)=xtbup
    tbdw(ilon,ilat,:)=xtbdw
	       
	enddo !ilon
	enddo !ilat


	write(str,9001) lyear,imon,iday,ihour
    9001 format('atmos_ncep_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
    outfile = trim(outpath)//trim(str)
    
    open(unit=4,form='binary',file=outfile,action='write')
	write(4) ao
    write(4) av
    write(4) al
	write(4) tran
	write(4) tbup
	write(4) tbdw
	close(4)		
	
	stop ' norm end pgm MAKE_L_BAND_ATM'
	end program MAKE_L_BAND_ATM   