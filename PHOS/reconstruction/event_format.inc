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
