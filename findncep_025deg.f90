!   11/2021 (AManaster)
!   Changed to read in 0.25 degree NCEP files. Changed 'ilon' and 'ilat'
!   to indicies between 1-1440 and 1-721, respectively instead of integers
!   representative of the actual lon and lat values (cannot for 0.25 degrees).


!     03/04/2013 
!	changed path from Nasserver5 to ops1p
!
!     01/05/2006
!     was: if t(ipr) < 100 k : ierr=-10, bail out
!     changed: if t(ipr) < 100 k, set rhocwat=0.0
!
!     01/28/2003
!     rhol = liquid cloud water density 
!
!     ordered ncep profiles, surtep, wind, pwat ,colvap
!     
!     input:
!     nmax: maximum number of levels
!     lyear:   long year (e.g. 1999)
!     imon:    month
!     iday:    day of the month 
!     ihour: uct hour
!     ilat: latitude  [between -90 and 90]
!     ilon: longitude	[between 0 and 360]


!    output:
!     colvap: ncep water vapor integrated [mm]
!     colwat: ncep columnar liquid cloud water integrated [mm]
!     pwat:   ncep value for precipitable water [mm]
!     cwat:   ncep value for columnar cloud water	(liquid + ice) [mm]
!     p: air pressure profile	   [mb] (0:nmax)   ordered
!     t: air temperature profile [k]  (0:nmax)
!     z: elevation profile [m]
!     pv: water vapor pressure profile [mb]  (0:nmax)
!     rhov: water vapor density [g /m**3]
!     rhol: liquid water density [g /m**3] 
!     ibegin: index of surface level
 

    subroutine findncep_025deg(lyear,imon,iday,ihour,ilat,ilon,  colvap,colwat,pwat,cwat,p,t,pv,rhov,rhol,z,ibegin)
  	implicit none

    integer(4), parameter                       :: nmax =26
    real(4), parameter                          :: r_e = 6371000, rd=287.05, epsilon = 0.622 

	integer(4), intent(in)                      :: lyear,imon,iday,ihour,ilat,ilon	
	integer(4), intent(out)                     :: ibegin
	
	real(4), intent(out)                        :: colvap,colwat,pwat,cwat
	real(4),    dimension(0:nmax), intent(out)  :: p,t,pv,rhov,rhol,z 
	
	integer(4)                                  :: ipr
	real(4),    dimension(0:nmax)               :: rhocwat,rhoi,clwmr,hgt,rh
	
	real(4)                                     :: p_sfc
	real(4)                                     :: p0,rhoair


    p(0:nmax) = (/ 0., &
	     1000.,975.,950.,925.,900.,850.,800.,750.,700.,650.,600., &
 	     550.,500.,450.,400.,350.,300.,250.,200.,150.,100., &
         70., 50., 30., 20. ,10.  /)


    call fnd_ncep_025deg(lyear,imon,iday,ihour,ilat,ilon,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)

	p0 = 0.01*p_sfc    ! hpa

    ibegin=-1
    do ipr=1,nmax
	    if(p(ipr).le.p0) then  !le is used rather than lt to be consistent with old sorting method
	        ibegin=ipr-1
	        exit
	    endif
	enddo
	if(ibegin.eq.-1) stop 'ibegin.eq.-1 in findncep, pgm stopped' 

	p(    ibegin) =p0
	t(    ibegin) =t(0)
	hgt(  ibegin) =hgt(0)
	rh(   ibegin) =rh(0)
	clwmr(ibegin) =clwmr(0)
  

    z = hgt * r_e / (r_e - hgt)
	if(z(ibegin) .ge. z(ibegin+1)) z(ibegin) = z(ibegin+1) - 0.1 

!	transform rh -> water vapor pressure and density
    call goff_gratch_vap(nmax,t,rh,p,  pv,rhov)

!   convert ncep clwmr to rhocwat
    clwmr(ibegin) = clwmr(ibegin+1) ! clwmr at surface  

	do ipr = 0,nmax
	    if (t(ipr) >= 100) then
	        rhoair = p(ipr)/(rd*t(ipr)) ! density of dry air
	        rhoair = rhoair*(1.0 - (1.0-epsilon)*pv(ipr)/p(ipr)) 
	        ! density of moist air (without cloud), unusual ncep definition for mixing ratio

            rhocwat(ipr) = 1.0e5*clwmr(ipr)*rhoair	! [g/m**3]
	        call cldwatice(t(ipr),rhocwat(ipr),  rhol(ipr),rhoi(ipr))  
	    else ! unphysical temperature
            rhocwat(ipr)=0.0
	        rhol(ipr) = 0.0
            rhoi(ipr) = 0.0
	    endif

	enddo

    where (rhov < 0) rhov=0.0
    where (rhol < 0) rhol=0.0

    call column(nmax-ibegin,z(ibegin:nmax),rhov(ibegin:nmax),2, colvap)
    call column(nmax-ibegin,z(ibegin:nmax),rhol(ibegin:nmax),1, colwat)

	colvap = colvap *1.e-3   ! g/m**2 -> mm	= kg/m**2
	colwat = colwat *1.e-3   ! g/m**2 -> mm	= kg/m**2
	
	return
	end subroutine findncep_025deg




	subroutine cldwatice(t,rhoc,  rhol,rhoi)
!     t: temperature [k]
!     rhoc:   cloud water density [arbitrary unit]
! 
!     rhol : liquid part [same unit as rhoc]
!     rhoi : ice    part [same unit as rhoc]
    implicit none
	real, parameter         :: twat = 273.15, tice = 253.15
	real, intent(in)        :: t,rhoc
	real(4), intent(out)    :: rhol,rhoi

!     water or ice 
	if (t >= twat) then
	         rhol = rhoc ! all liquid water
             rhoi = 0.0 
	else if (t <= tice) then 
	         rhol = 0.0       ! all ice
             rhoi = rhoc
	else
	         rhol = rhoc*(t-tice)/(twat-tice)
             rhoi = rhoc - rhol  
			 ! linear temperature interpolation in between
	endif
	
	return
	end subroutine cldwatice

!    01/23/2006
!    kyle changed subroutine yread
!    read fileheader between each profile level 

!
!    ncep variables
!    1/4 deg resolution
!    4 times daily: 00z 06z 12z 18z
! 
!    thomas meissner : apr 1999 
!
!    needs to be linked with time_routines.f 
!
!
!         varaible    level (subfolder)
!           sst               sfc    (0) 
!           t                 sfc    (0) 2 m above ground
!           rh                sfc    (0)
!           hgt               sfc    (0)
!           p                 sfc    (0)
!           pwat              col
!           cwat              col
!           t                 1:nmax  (prf)
!           rh                1:nmax   (prf)
!           hgt               1:nmax  (prf)
!           clwmr             1:nmax  (prf) 
!
!    named as year_month_day_xxz.dat :
!    binary file:
!    year month day hour (integer(4))
!    real*4 array 
!    grid: lat from +90 to -90 by -0.25 deg
!          lon from 0 to 360 by 0.25 deg
!    total size 4*4 + 360*181*4 = 260656  
!
!    input: lyear  (year, long form, e.g. 2005)                   integer
!           imon   (month, 1-12)                                  integer
!           iday (day of month, 1-31)                             integer
!           ihour   (hour of day)                                 integer
!           xlat   (latitude)                                     real(4)
!           xlon   (longitude)                                    real(4)
!    
!    output: var (interpolated variable) real(4) 
!                                   
!

      subroutine fnd_ncep_025deg(lyear,imon,iday,ihour,ilat,ilon,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
	  implicit none

      integer, parameter        ::  nmax =26

	  integer(4), intent(in)    ::  lyear,imon,iday,ihour,ilat,ilon

      real(4), intent(out)      ::  t(0:nmax),hgt(0:nmax),rh(0:nmax),clwmr(0:nmax),p_sfc,pwat,cwat

      integer(4), save          ::  lyearsv=-999,imonsv=-999,idaysv=-999,ihoursv=-999
	  
	  integer(4)                ::  j,k

      real(4), save             ::  t_map(1440,721,0:nmax), hgt_map(1440,721,0:nmax), rh_map(1440,721,0:nmax),clwmr_map(1440,721,0:nmax)
      real(4), save             ::  p_sfc_map(1440,721),pwat_map(1440,721),cwat_map(1440,721)

      if(lyear.ne.lyearsv .or. imon.ne.imonsv .or. iday.ne.idaysv .or. ihour.ne.ihoursv ) then
          lyearsv=lyear
          imonsv=imon
	      idaysv=iday
          ihoursv=ihour
          call read_ncep_025deg(lyear,imon,iday,ihour,  t_map,hgt_map,rh_map,clwmr_map,p_sfc_map,pwat_map,cwat_map)
	  endif

      j= 722 - ilat 
      ! flips indicies around correctly since NCEP goes from +90 to -90 e.g., ilat = 1 (-90 on RSS grid) is NCEP idx 721 (-90 on NCEP grid)
	  k= ilon

      t=        t_map(k,j,:)                                                              
      hgt=    hgt_map(k,j,:)
      rh =     rh_map(k,j,:)                                                              
      clwmr=clwmr_map(k,j,:)  
	                                                            
      p_sfc=p_sfc_map(k,j)                                                            
      pwat = pwat_map(k,j)                                                            
      cwat = cwat_map(k,j)                                                            

      return
      end subroutine fnd_ncep_025deg



      subroutine read_ncep_025deg(lyear,imon,iday,ihour,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
      implicit none

      integer(4), parameter     ::   nmax = 26, nrh=21

      integer(4), intent(in)    ::   lyear,imon,iday,ihour
      real(4), intent(out)      ::   t(1440,721,0:nmax),hgt(1440,721,0:nmax),rh(1440,721,0:nmax),clwmr(1440,721,0:nmax)
	  real(4), intent(out)      ::   p_sfc(1440,721),pwat(1440,721),cwat(1440,721)

	  character(len=200)        ::   filename
      logical(4)                ::   lexist

      integer(4)                ::   kyear,kmon,kday,khour     
      integer(4)                ::   ilevel,icase,isleep

        ! location and name of NCEP files
        ! specified by user
        9001 format('sample_data\NCEP_PROF_TMP_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9002 format('sample_data\NCEP_SURF_TMP_2m_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9003 format('sample_data\NCEP_PROF_RH_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9004 format('sample_data\NCEP_SURF_RH_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9005 format('sample_data\NCEP_PROF_HGT_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9006 format('sample_data\NCEP_SURF_HGT_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9007 format('sample_data\NCEP_SURF_PRES_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9008 format('sample_data\NCEP_SURF_PWAT_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9009 format('sample_data\NCEP_PROF_CLWMR_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
        9010 format('sample_data\NCEP_SURF_CWAT_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')

        do icase=1,10
        if(icase.eq. 1) write(filename,9001) lyear,imon,iday,ihour
        if(icase.eq. 2) write(filename,9002) lyear,imon,iday,ihour
        if(icase.eq. 3) write(filename,9003) lyear,imon,iday,ihour
        if(icase.eq. 4) write(filename,9004) lyear,imon,iday,ihour
        if(icase.eq. 5) write(filename,9005) lyear,imon,iday,ihour
        if(icase.eq. 6) write(filename,9006) lyear,imon,iday,ihour
        if(icase.eq. 7) write(filename,9007) lyear,imon,iday,ihour
        if(icase.eq. 8) write(filename,9008) lyear,imon,iday,ihour
        if(icase.eq. 9) write(filename,9009) lyear,imon,iday,ihour
        if(icase.eq.10) write(filename,9010) lyear,imon,iday,ihour

        isleep=1

        10 continue
 	    inquire(file=filename,exist=lexist)

	    if(.not.lexist) then
	    if(isleep.eq.1) write(*,1001) trim(filename)
        1001 format(' waiting for ',a100)
        isleep=0
        call sleepqq(60000) !wait one minute 
	    goto 10
	    endif

        call sleepqq(10000) !wait 10 seconds for possible file download to complete 

        open(unit=3,file=filename,status='old',form='binary',action='read')

        read(3) kyear, kmon, kday, khour
        if(kyear.ne.lyear .or. kmon.ne.imon .or. kday.ne.iday .or. khour.ne.ihour) then
            write(*,*) filename
            write(*,*) kyear,kmon,kday,khour
            write(*,*) lyear,imon,iday,ihour
            write(*,*) ' header error in read_ncep_025deg, pgm stopped'
            stop
        endif


        if(icase.eq.1) then
        do ilevel=1,nmax  
        read(3) t(:,:,ilevel)
        enddo
	    endif

        if(icase.eq.2) then
        read(3) t(:,:,0)
	    endif

        if(icase.eq.3) then
        do ilevel=1,nrh   
        read(3) rh(:,:,ilevel)
        enddo
	    rh(:,:,nrh+1:nmax) = 0.0
	    endif

        if(icase.eq.4) then
        read(3) rh(:,:,0)
	    endif

        if(icase.eq.5) then
        do ilevel=1,nmax 
        read(3) hgt(:,:,ilevel)
        enddo
	    endif

        if(icase.eq.6) then
        read(3) hgt(:,:,0)
	    endif

        if(icase.eq.7) then
        read(3) p_sfc
	    endif

        if(icase.eq.8) then
        read(3) pwat
	    endif

        if(icase.eq.9) then
	    do ilevel=1,nrh 
	    read(3) clwmr(:,:,ilevel)
        enddo
	    clwmr(:,:,nrh+1:nmax) = 0.0
        clwmr(:,:,0) = 0.0
	    endif

        if(icase.eq.10) then
        read(3) cwat
	    endif

        close(3)

      enddo !icase

    return
	end subroutine read_ncep_025deg 
