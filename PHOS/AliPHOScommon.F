*CMZ :  2.03/01 18/09/98  11.43.30  by  Federico Carminati
*-- Author :
      SUBROUTINE PHOS_DATA
*KEEP,SCXXCOM.
	parameter (NGp=1000,nsps=10,nvertmax=1000)
        COMMON /RCGAMMA/KG,MW(ngp),ID(ngp),JD(ngp),E(ngp),E4(ngp),
     ,  XW(ngp),YW(ngp),ES(nsps,ngp),ET(nsps,ngp),ISsd(ngp),
     ,  IGDEV(ngp),ZGDEV(ngp),sigexy(3,ngp),Emimx(2,nsps,ngp),
     ,  kgfix,igfix(ngp),cgfix(3,ngp),sgfix(3,ngp),hiw(ngp),
     ,  wsw(nsps,ngp),h1w(ngp),h0w(ngp),raxay(5,ngp),
     ,  sigmaes0(nsps,ngp),dispeces(nsps,ngp),
     ,  igamvert(ngp)


      integer*4 crystals_amount_max,crystals_in_matrix_amount_max,
     +          crystals_matrix_amount_max
      parameter (crystals_matrix_amount_max=4)
      parameter (crystals_in_matrix_amount_max=40000)
      parameter (crystals_amount_max =crystals_matrix_amount_max*
     +                                crystals_in_matrix_amount_max)

* All units are in GeV, cm, radian
      real       crystal_amplitudes_unit, radius_unit,
     +           crystal_size_unit, crystal_length_unit,
     +           matrix_coordinate_Z_unit, matrix_coordinate_PHI_unit
      integer    crystal_amplitudes_in_units_min
      parameter (crystal_amplitudes_in_units_min        = 1)
      parameter (crystal_amplitudes_unit                = 0.001 ) ! 1.0  MeV
      parameter (radius_unit                            = 0.1   ) ! 0.1  cm
      parameter (crystal_size_unit                      = 0.01  ) ! 0.01 cm
      parameter (crystal_length_unit                    = 0.01  ) ! 0.01 cm
      parameter (matrix_coordinate_Z_unit               = 0.1   ) ! 0.1  cm
      parameter (matrix_coordinate_PHI_unit             = 1e-4  ) ! 1e-4 radian

      integer*2 crystals_matrix_amount_PHOS, crystal_matrix_type,
     +          amount_of_crystals_on_Z, amount_of_crystals_on_PHI,
     +          crystals_amount_with_amplitudes, crystals_amplitudes_Iad
      integer*4 event_number

      real      radius, crystal_size, crystal_length,
     +          matrix_coordinate_Z, matrix_coordinate_PHI

      real      crystals_amplitudes, crystals_energy_total
      integer   event_file_unit_number

      common /common_for_event_storing/
     + ! Event-independent information
     +        crystals_matrix_amount_PHOS,
     +        crystal_matrix_type,
     +        amount_of_crystals_on_Z,
     +        amount_of_crystals_on_PHI,
     +        radius,
     +        crystal_size,
     +        crystal_length,
     +        matrix_coordinate_Z     (crystals_matrix_amount_max),
     +        matrix_coordinate_PHI   (crystals_matrix_amount_max),
     +
     + ! Event-dependent information
     +        event_number,
     +        crystals_amount_with_amplitudes
     +                                (crystals_matrix_amount_max),
     +        crystals_amplitudes_Iad (2,crystals_in_matrix_amount_max,
     +                                 crystals_matrix_amount_max),
     +
     + ! These information don't store in data file
     +        crystals_amplitudes     (crystals_amount_max),
     +        crystals_energy_total,
     +        event_file_unit_number



        INTEGER MAXCRAD
        PARAMETER (MAXCRAD=100)
        INTEGER PHOSsize,PHOS_Ndiv_magic
        REAL    PHOSflags,PHOScell,PHOSradius,PHOSCPV,
     +          PHOScradlesA,PHOSTXW,PHOSAIR,PHOSFTI,
     +          PHOSextra, PHOSangle
        COMMON /PHOS_PARS/ PHOSflags(9),
     +                           PHOScell(9),PHOSradius,PHOSCPV(9),
     +                           PHOSsize(3), PHOScradlesA,
     +                           PHOSTXW(3),PHOSAIR(3),PHOSFTI(4),
     +                           PHOSextra(9), PHOSangle(MAXCRAD),
     +                           PHOS_Ndiv_magic

*KEND.
      END
