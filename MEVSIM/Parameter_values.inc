CCC   Set array dimension sizes here:

      integer npid,nmax_integ,n_mult_max_steps,nflowterms
      parameter (npid = 30)         ! max # of particle ID types
      parameter (nmax_integ = 100)  ! max # integration steps in parameter
CCC                                 ! variance calculation.
      parameter (n_mult_max_steps = 1000)
CCC                                 ! max # integration steps in multiplicity
CCC                                 ! variance calculation (this must be an
CCC                                 ! even integer).
      parameter (nflowterms = 6)    ! max # of terms in the anisotropic
CCC                                 ! flow model for azimuthal (phi angle)
CCC                                 ! dependence.

