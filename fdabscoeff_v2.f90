!     TM 09/19/2016 
!     changed O2 absorption back to Liebe 1992
!  
!     same as 'O:\skytemp3\atms_abs_routines.f' dated July 17 2015 except
!     1.  routine fdabscoeff is added
!     2.  fdaray is removed
!     3.  ga=5.6e-3*(pdry+1.1*pwet)*1.097687*tht**0.8   !1.102393=(300./262.6)**0.7, 262.6 is avg teff from zonal_adjustment.m


!	input:
!     oxygen absorption from rosenkranz
!     freq  frequency [in ghz]
!     p      pressure [in h pa]
!     t      temperature [in k]
!     pv     water vapor pressure  [in hpa]
!
!     output:	
!     av          water vapor absorption coefficients [neper/km]
!     ao          oxygen absortption coefficient		[neper/km]

    subroutine fdabscoeff(freq,p,t,pv, av,ao)
	implicit none

    real(4), parameter :: xnaper=0.2302585094 !convert db/km to naper/km

    real(4) freq,p,t,pv
	real(4)	av,ao
	real(4) gamoxy,gamh2o

	call fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)  !gamoxy is db/km
    call abh2o_rk_modified(     p,t,pv,freq, gamh2o)  !gamh2o is db/km

	ao=xnaper*gamoxy
	av=xnaper*gamh2o 
 
    return
    end subroutine fdabscoeff


!     june 11 2015 changed july 17 2015.  oxyopc adjustment above 37 ghz changed slightly. see 'O:\gmi\abs_cal\memo20.txt'

!     this module contains the oxygen absoprtion routine, the vapor absorption routine and fdaray
!     there were used for the april-june 2009 skytemp update.  see 'memo8.txt'


!     ================================================================================================================
!     ========================== modified version of Liebe 1992 oxygen model =========================================
!     ================================================================================================================

!     This is from Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford, 1992
!     coded: June 2009 1992 by f.wentz
!           inputs: t, temperature (k)
!                   p, total pressure (mb)
!                   pv, water vapor pressure (mb)
!                   freq, frequency (ghz)
!           output: gamoxy, oxygen absorption coefficient (db/km)


!     It is the same as fdabsoxy_1989 except for the a5 and a6 coefs for finding delta have different values.
!     Also this 1992 versions says delta is proprotional to total pressure p rather than pdry. 
!     Also note in this routine, the 1.e-3 scaling is done in the startup block.

!     compared to abo2_rk (the code Rosenkranz sent us), this routine gives very similar results if you set the apterm to 0.
!     for my freqs 6-85 ghz, the largest dif was 0.0003 at the coldest vapor at 85.5 GHz.
!     the apterm adds 0.003 at 85 ghz
!     Apart from the apterm, you can essentially say fdabsoxy_1992 and abo2_rk are the same for my purposes.

!     this routine has been modified in the following ways (see 'memo8.txt')
!     1.  non-resonance continuum temperature coef changed from 0.8 to 1.5
!     2.  a p*p continuum was added
!     these modifications were done June 22 2009 

    subroutine fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)
	implicit none

    integer(4), parameter :: nlines=44
    integer(4) istart,i
    real(4) p,t,pv,freq,gamoxy
	real(4) tht,pwet,pdry,ga,gasq,delta,rnuneg,rnupos,ff,zterm,apterm,sftot,xterm
    real(4) h(6,nlines),f0(nlines),a1(nlines),a2(nlines),a3(nlines),a4(nlines),a5(nlines),a6(nlines)

    real(8) sum

    data istart/1/
	data a4/38*0., 6*0.6/

!          freq          a1      a2       a3       a5          a6
      data h/ &
       50.474238,    0.94e-6,  9.694,  8.60e-3,  0.210,  0.685, &
       50.987749,    2.46e-6,  8.694,  8.70e-3,  0.190,  0.680, &
       51.503350,    6.08e-6,  7.744,  8.90e-3,  0.171,  0.673, &
       52.021410,   14.14e-6,  6.844,  9.20e-3,  0.144,  0.664, &
       52.542394,   31.02e-6,  6.004,  9.40e-3,  0.118,  0.653, &
       53.066907,   64.10e-6,  5.224,  9.70e-3,  0.114,  0.621, &
       53.595749,  124.70e-6,  4.484, 10.00e-3,  0.200,  0.508, &
       54.130000,  228.00e-6,  3.814, 10.20e-3,  0.291,  0.375, &
       54.671159,  391.80e-6,  3.194, 10.50e-3,  0.325,  0.265, &
       55.221367,  631.60e-6,  2.624, 10.79e-3,  0.224,  0.295, &
       55.783802,  953.50e-6,  2.119, 11.10e-3, -0.144,  0.613, &
       56.264775,  548.90e-6,  0.015, 16.46e-3,  0.339, -0.098, &
       56.363389, 1344.00e-6,  1.660, 11.44e-3, -0.258,  0.655, &
       56.968206, 1763.00e-6,  1.260, 11.81e-3, -0.362,  0.645, &
       57.612484, 2141.00e-6,  0.915, 12.21e-3, -0.533,  0.606, &
       58.323877, 2386.00e-6,  0.626, 12.66e-3, -0.178,  0.044, &
       58.446590, 1457.00e-6,  0.084, 14.49e-3,  0.650, -0.127, &
       59.164207, 2404.00e-6,  0.391, 13.19e-3, -0.628,  0.231, &
       59.590983, 2112.00e-6,  0.212, 13.60e-3,  0.665, -0.078, &
       60.306061, 2124.00e-6,  0.212, 13.82e-3, -0.613,  0.070, &
       60.434776, 2461.00e-6,  0.391, 12.97e-3,  0.606, -0.282, &
       61.150560, 2504.00e-6,  0.626, 12.48e-3,  0.090, -0.058, &
       61.800154, 2298.00e-6,  0.915, 12.07e-3,  0.496, -0.662, &
       62.411215, 1933.00e-6,  1.260, 11.71e-3,  0.313, -0.676, &
       62.486260, 1517.00e-6,  0.083, 14.68e-3, -0.433,  0.084, &
       62.997977, 1503.00e-6,  1.665, 11.39e-3,  0.208, -0.668, &
       63.568518, 1087.00e-6,  2.115, 11.08e-3,  0.094, -0.614, &
       64.127767,  733.50e-6,  2.620, 10.78e-3, -0.270, -0.289, &
       64.678903,  463.50e-6,  3.195, 10.50e-3, -0.366, -0.259, &
       65.224071,  274.80e-6,  3.815, 10.20e-3, -0.326, -0.368, &
       65.764772,  153.00e-6,  4.485, 10.00e-3, -0.232, -0.500, &
       66.302091,   80.09e-6,  5.225,  9.70e-3, -0.146, -0.609, &
       66.836830,   39.46e-6,  6.005,  9.40e-3, -0.147, -0.639, &
       67.369598,   18.32e-6,  6.845,  9.20e-3, -0.174, -0.647, &
       67.900867,    8.01e-6,  7.745,  8.90e-3, -0.198, -0.655, &
       68.431005,    3.30e-6,  8.695,  8.70e-3, -0.210, -0.660, &
       68.960311,    1.28e-6,  9.695,  8.60e-3, -0.220, -0.665, &
      118.750343,  945.00e-6,  0.009, 16.30e-3, -0.031,  0.008, &
      368.498350,   67.90e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
      424.763124,  638.00e-6,  0.044, 19.16e-3,  0.0,    0.0,   &
      487.249370,  235.00e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
      715.393150,   99.60e-6,  0.145, 18.10e-3,  0.0,    0.0,   &
      773.839675,  671.00e-6,  0.130, 18.10e-3,  0.0,    0.0,   &
      834.145330,  180.00e-6,  0.147, 18.10e-3,  0.0,    0.0/  

      if(istart.eq.1) then
      istart=0
      f0(:)=h(1,:)
      a1(:)=h(2,:)/h(1,:)
      a2(:)=h(3,:)
      a3(:)=h(4,:)
      a5(:)=0.001*h(5,:)
      a6(:)=0.001*h(6,:)
      endif

      tht = 300/t
      pwet=0.1*pv
      pdry=0.1*p-pwet
      xterm=1-tht

      sum = 0.
      do i=1,nlines
      ga = a3(i)*(pdry*tht**(0.8-a4(i)) + 1.1*tht*pwet)
      gasq=ga*ga
	  delta=(a5(i) + a6(i)*tht)*p*tht**0.8
      rnuneg = f0(i)-freq
      rnupos = f0(i)+freq
      ff = (ga-rnuneg*delta)/(gasq+rnuneg**2) +  (ga-rnupos*delta)/(gasq+rnupos**2)
      sum = sum + ff*a1(i)*exp(a2(i)*xterm)
	  enddo
      if(sum.lt.0) sum=0

!     add nonresonant contribution

!     ga=5.6e-3*(pdry+1.1*pwet)*tht**0.8  
!x    ga=5.6e-3*(pdry+1.1*pwet)*tht**1.5  !modification 1
!     changed O2 absorption back to Liebe
      ga=5.6e-3*(pdry+1.1*pwet)*1.097687*tht**0.8   
      !1.102393=(300./262.6)**0.7, 262.6 is avg teff from zonal_adjustment.m
      
      zterm=ga*(1.+(freq/ga)**2)
      apterm=1.4e-10*(1-1.2e-5*freq**1.5)*pdry*tht**1.5
      if(apterm.lt.0) apterm=0
      sftot=pdry*freq*tht**2 * (tht*sum + 6.14e-4/zterm + apterm)

      gamoxy=0.1820*freq*sftot
!x    if(freq.gt.37) gamoxy=gamoxy + 0.1820*43.e-10 *pdry**2*tht**3*(freq-37.)**1.7  !prior to 7/17/2015
      if(freq.gt.37) gamoxy=gamoxy + 0.1820*26.e-10 *pdry**2*tht**3*(freq-37.)**1.8  !implemented 7/17/2015.

      return
      end subroutine fdabsoxy_1992_modified


!     ================================================================================================================
!     ========================== modified version of Rosenkranz water vapor model ====================================
!     ================================================================================================================


!   purpose- compute absorption coef in atmosphere due to water vapor
!
!   calling sequence parameters-
!   specifications
!   name    units    i/o  descripton            valid range
!
!      t       kelvin    i   temperature
!      p       millibar  i   pressure              .1 to 1000
!      f       ghz       i   frequency             0 to 800
!      gamh2o  db/km     o   absorption coefficient
!
!   references-
!   p.w. rosenkranz, radio science v.33, pp.919-928 (1998).
!
!   line intensities selection threshold=
!   half of continuum absorption at 1000 mb.
!   widths measured at 22,183,380 ghz, others calculated.
!   a.bauer et al.asa workshop (sept. 1989) (380ghz).
!
!   revision history-
!   date- oct.6, 1988  p.w.rosenkranz - eqs as publ. in 1993.
!          oct.4, 1995  pwr- use clough's definition of local line
!                   contribution,  hitran intensities, add 7 lines.
!          oct. 24, 95  pwr -add 1 line.
!          july 7, 97   pwr -separate coeff. for self-broadening,
!                       revised continuum.
!          dec. 11, 98  pwr - added comments

!   the routine is a modified version of abh2o_rk_reformat. 
!   this routine has been modified in the following three ways (see 'memo8.txt')
!   1.  b1(1)=1.01*b1(1)  :22 ghz line strength increase slightly
!   2.  22 ghz line shape below 22 ghz has been modified
!   3.  foreign and self broadening continuum has been adjusted 
!   these modification were done June 22 2009 
 




      subroutine abh2o_rk_modified(p,t,pv,freq,  gamh2o)
      implicit none
  
      integer(4), parameter :: nlines=15

      integer(4) istart,i
      real(4) t,p,freq
      real(4) b1(nlines),b2(nlines),b3(nlines),f0(nlines),b4(nlines),b5(nlines),b6(nlines)
      real(4) pv,s,base,gamh2o
	  real(4) tht,pwet,pdry,ga,gasq,sftot,xterm,rnuneg,rnupos
	  real(8) sum

	  real(4) chi,chisq,freqsq,f0sq,u

      data istart/1/

!     line frequencies:
      data f0/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
      443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, &
      620.7008, 752.0332, 916.1712/
!     line intensities at 300k:
      data b1/ .1310e-13, .2273e-11, .8036e-13, .2694e-11, .2438e-10, &
              .2179e-11, .4624e-12, .2562e-10, .8369e-12, .3263e-11, .6659e-12, &
              .1531e-08, .1707e-10, .1011e-08, .4227e-10/
!     t coeff. of intensities:
      data b2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
!     air-broadened width parameters at 300k:
      data b3/.0281, .0281, .023, .0278, .0287, .021, .0186, .0263, .0215, .0236, .026, .0321, .0244, .0306, .0267/
!     self-broadened width parameters at 300k:
      data b5/.1349, .1491, .108, .135, .1541, .090, .0788, .1275, .0983, .1095, .1313, .1320, .1140, .1253, .1275/
!     t-exponent of air-broadening:
      data b4 /.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, .71, .68, .70/
!     t-exponent of self-broadening:
      data b6/.61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, 1.0, .68, .84, .78/

      if(istart.eq.1) then
	  istart=0
      b1=1.8281089E+14*b1/f0**2
	  b5=b5/b3  !convert b5 to Liebe notation
      b1(1)=1.01*b1(1)  !modification 1
	  endif
 
      if(pv.le.0.) then
      gamh2o=0
	  return
      endif
 

      pwet=0.1*pv
      pdry=0.1*p-pwet
      tht = 300./t
      xterm=1-tht
	  freqsq=freq*freq

      sum = 0.
      do i=1,nlines
	  f0sq=f0(i)*f0(i)
      ga=b3(i)*(pdry*tht**b4(i) + b5(i)*pwet*tht**b6(i))
      gasq = ga*ga
      s = b1(i)*exp(b2(i)*xterm)
      rnuneg = f0(i)-freq
      rnupos = f0(i)+freq
      base = ga/(562500. + gasq)  !use clough's definition of local line contribution

      if(i.ne.1) then
      if(abs(rnuneg).lt.750) sum = sum + s*(ga/(gasq + rnuneg**2) - base)
      if(abs(rnupos).lt.750) sum = sum + s*(ga/(gasq + rnupos**2) - base)

      else
	  chi=0

	  if(freq.lt.19) then
	  u=abs(freq-19.)/16.5
	  if(u.lt.0) u=0
	  if(u.gt.1) u=1
	  chi=ga*u*u*(3-2*u)  !modification 2
	  endif
	
	  chisq=chi*chi
      sum=sum +     s*2*((ga-chi)*freqsq + (ga+chi)*(f0sq+gasq-chisq))/((freqsq-f0sq-gasq+chisq)**2 + 4*freqsq*gasq)
	  endif

 	  enddo
      if(sum.lt.0) sum=0
      
!x    sftot=pwet*freq*tht**3.5*(sum +     1.2957246e-6*pdry/tht**0.5 +                    4.2952193e-5*pwet*tht**4)
      sftot=pwet*freq*tht**3.5*(sum + 1.1*1.2957246e-6*pdry/tht**0.5 + 0.425*(freq**0.10)*4.2952193e-5*pwet*tht**4) !modification 3

      gamh2o=0.1820*freq*sftot
	  
	  return
      end subroutine abh2o_rk_modified

