************************************************************************

	subroutine reconsfirst(ampshum,s1gev)

        common /comgeom/igeomflag
        include 'event_format.inc'
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'
	include 'comgam.for'
	include 'comggen.for'

	data ifl/0/
	if(ifl.eq.0) then
		call iniwl(ampshum,s1gev)
	ifl=1
	endif
        igeomflag=1     ! Rectangular geometry

        if( amount_of_crystals_on_Z.NE.104 .OR.
     +      amount_of_crystals_on_PHI.NE.88 ) THEN
	   print *, 'Nz x Nphi ',amount_of_crystals_on_Z,
     +         ' x ',amount_of_crystals_on_PHI
         stop 'Reconstruction: Cradle size must be Nz x Nphi = 104 x 88'
	 endif

        nmmmax=crystals_matrix_amount_PHOS
        if(crystals_matrix_amount_PHOS.ne.1) then
          write (*,*) ' WRONG GEOMETRY: WRONG NUMBER OF CRADLES '
          stop
        endif

        call clrcompart ! nbpart=0

        i=1     ! Cradle number
	nggl=crystals_amount_with_amplitudes(i)
	call vzero(wl,nt)

	do n=1,nggl
	  eggl(n)=crystals_amplitudes_Iad(1,n,i)*crystal_amplitudes_unit
     	  mggl(n)=crystals_amplitudes_Iad(2,n,i)+1

	  m=mggl(n)
	  ewlm=eggl(n)
	  wl(m)=ewlm
	enddo


	eporog =0.005
	call extrgw(nt,wl,ipgl,mgl,egl,eporog)
	nglw(1)=ipgl
	call vzero(wl,nt)

	call mygl

	return
	end

************************************************************************
