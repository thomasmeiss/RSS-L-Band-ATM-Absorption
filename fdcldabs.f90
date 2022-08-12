!     liquid cloud water absorption
!     rayleigh 
!     freq:      frequency [ghz]
!     t:         temperature [k]
!     rhol:      liquid cloud water density [g/m**3]
!     output:
!     al:        cloud water absorption coefficient [neper/km]

    subroutine fdcldabs(freq,t,rhol,   al)
    implicit none 

	real(4) :: freq,t,rhol,rhol0,al,wavlen,c=29.979,pi=3.14159
	complex(4) :: permit

	rhol0 = 1.0e-6*rhol ![g/cm**3]
	call meissner(freq,t,0.0,   permit)	
	wavlen = c/freq
	al = (6.0*pi*rhol0/wavlen)*aimag((1.0-permit)/(2.0+permit))  ! nepers/cm
	al = 1.0e5*al ! nepers/km
	
	return
	end subroutine fdcldabs

