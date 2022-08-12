!     adapted for SMAP
!     EIA between 39 - 41 deg


!     ncep 0.25 deg profiles 
!    
!     rss vapor absorption 
!     liebe o2 absorption
!     meissner - wentz dielectric model for cloud water
!   
!     icld = 0: no clouds
!     icld = 1: ncep cloud density, only liquid part 
!

	subroutine get_atm(icld,lyear,imonth,iday,ihour,ilat,ilon,freq,   tbup,tbdw,tran,ao,av,al)
	implicit none
      
	integer(4), parameter   ::  nmax = 26
    integer(4), intent(in)  ::  icld,lyear,imonth,iday,ihour,ilat,ilon
	real(4), intent(in)     ::  freq

	real(4), intent(out)    ::  tbup(39:41),tbdw(39:41),tran(39:41)
	real(4), intent(out)    ::  ao,av,al

 	integer(4)              ::  ipr,ibegin,itht
	real(4)                 ::  tht,vap,cld,pwat,cwat
	real(4)                 ::  t(0:nmax),p(0:nmax),pv(0:nmax),z(0:nmax),rhov(0:nmax),rhol(0:nmax),rhol0(0:nmax)
	real(4)                 ::  abh2o(0:nmax),abo2(0:nmax),abcld(0:nmax),tabs(0:nmax) 

    !     ncep parameters
    call findncep_025deg(lyear,imonth,iday,ihour,ilat,ilon, vap,cld,pwat,cwat,p,t,pv,rhov,rhol0,z,ibegin)

	if (icld == 0) then ! no cloud
	    rhol= 0.0
	else
        rhol= rhol0
    endif

    !    atmospheric absorption
    !    o2 from rosenkranz 
    !    h2o vapor from amsr atbd
    !    cloud with meissner dielectric constant 


	do ipr  = ibegin,nmax

        call fdabscoeff(freq,p(ipr),t(ipr),pv(ipr),  abh2o(ipr),abo2(ipr)) ! neper/km

	    if (rhol(ipr) > 1.0e-7) then
	        call fdcldabs(freq,t(ipr),rhol(ipr),   abcld(ipr)) ! nepers/km
	    else
	        abcld(ipr) = 0.0
	    endif

	    abh2o(ipr)=abh2o(ipr)/1000.0	! neper/m
	    abo2(ipr) = abo2(ipr)/1000.0    ! neper/m
	    abcld(ipr)=abcld(ipr)/1000.0    ! neper/m

	enddo ! ipr

    !     vertical integrals

	call column(nmax-ibegin,z(ibegin:nmax),abh2o(ibegin:nmax),2, av)
	call column(nmax-ibegin,z(ibegin:nmax), abo2(ibegin:nmax),2, ao)
	call column(nmax-ibegin,z(ibegin:nmax),abcld(ibegin:nmax),1, al)

    !     total absorption 
	tabs(ibegin:nmax) = abh2o(ibegin:nmax) + abo2(ibegin:nmax) + abcld(ibegin:nmax)


    do itht=39,41
	    tht=float(itht)
	    call atm_tran(nmax-ibegin,tht,t(ibegin:nmax),z(ibegin:nmax),tabs(ibegin:nmax), tran(itht),tbdw(itht),tbup(itht))
	enddo

	return
	end subroutine get_atm


    ! 	compute atmospheric downwelling and upwelling brightness temperatures
    !	and upward transmittance at each pressure level (altitude) 

    !	input:
    !     nlev           number of atmosphere levels
    !     tht            earth incidence angle [in deg]
    !     tabs(0:nlev)   atmosphric absorptrion coefficients [nepers/m]
    !     t(0:nlev)      temperature profile[in k]
    !	  z(0:nlev)      elevation (m) 


    !     output:
    !     tran        	total atmospheric transmission
    !     tbdw			downwelling brightness temperature t_bd [in k]
    !     tbup			upwelling   brightness temperature t_bu [in k]  

    
    subroutine atm_tran(nlev,tht,t,z,tabs,  tran,tbdw,tbup)
    implicit none

	real(4), parameter          :: re=6378.135, delta=0.00035

	integer(4), intent(in)      :: nlev
	real(4), intent(in)         :: t(0:nlev),z(0:nlev),tabs(0:nlev)
	real(4), intent(out)        :: tran, tbdw, tbup
	
	real(4)                     :: opacty(nlev),tavg(nlev),ems(nlev)
	real(4)                     :: sumop, sumdw, sumup, tbavg, dsdh, tht

    integer(4)                  :: i

	dsdh = (1.0+delta)/sqrt(cosd(tht)**2 + delta*(2+delta))  

    do i=1,nlev
	  opacty(i)=-dsdh*0.5*(tabs(i-1)+tabs(i))*(z(i)-z(i-1)) 
      tavg(i)  =0.5*(t(i-1)+t(i))
      ems(i)   =1.-exp(opacty(i))
	enddo

    sumop=0 
	sumdw=0
    do i=1,nlev
      sumdw=sumdw+(tavg(i)-t(1))*ems(i)*exp(sumop) 
      sumop=sumop+opacty(i)
	enddo

    sumop=0 
	sumup=0.
    do i=nlev,1,-1
      sumup=sumup+(tavg(i)-t(1))*ems(i)*exp(sumop)
      sumop=sumop+opacty(i)
	enddo

    tran=exp(sumop)
    tbavg=(1.-tran)*t(1)
    tbdw=tbavg+sumdw
    tbup=tbavg+sumup

	return
	end subroutine atm_tran