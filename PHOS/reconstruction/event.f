*****************************************************************************

      subroutine print_event_information
        implicit none
        include 'event_format.inc'

        write (*,*)'***********************************************'
        write (*,*)'************** event_number = ',
     +             event_number,' *******'
        write (*,*)'crystals_matrix_amount_PHOS =',
     +          crystals_matrix_amount_PHOS
        write (*,*)'crystal_matrix_type         =',
     +              crystal_matrix_type
        write (*,*)'amount_of_crystals_on_Z     =',
     +              amount_of_crystals_on_Z
        write (*,*)'amount_of_crystals_on_PHI   =',
     +              amount_of_crystals_on_PHI
        write (*,*)'radius                      =',radius
        write (*,*)'crystal_size                =',crystal_size
        write (*,*)'crystal_length              =',crystal_length

        call print_event_data

      end

*****************************************************************************

      subroutine print_event_data
        implicit none
        include 'event_format.inc'

        integer i,j,k

        do i=1,crystals_matrix_amount_PHOS
          write (*,*)'---- matrix number ',i,
     +           '  Z = ', matrix_coordinate_Z(i),
     +           '  PHI = ', matrix_coordinate_PHI(i)
          write (*,*) 'crystals_amount_with_amplitudes=',
     +             crystals_amount_with_amplitudes(i)

          do j=1,crystals_amount_with_amplitudes(i)
            write (*,*)j,'  amplitude/address=',
     +             crystals_amplitudes_Iad(1,j,i),
     +             crystals_amplitudes_Iad(2,j,i)
          enddo
        enddo

      end

*******************************************************************

      subroutine save_event_to_file
        implicit none
        include 'event_format.inc'

        integer n,i,j
        integer*2 integer_radius,
     +            integer_crystal_size, integer_crystal_length,
     +            integer_coordinate_Z, integer_coordinate_PHI,
     +            event_number__low_bytes, event_number__high_bytes

c        write (*,*)'Event saving...'

        integer_radius          = radius         /radius_unit
        integer_crystal_size    = crystal_size   /crystal_size_unit
        integer_crystal_length  = crystal_length /crystal_length_unit

        event_number__high_bytes = event_number/32767
        event_number__low_bytes  = event_number -
     +                             event_number__high_bytes*32767

        write(event_file_unit_number,err=1),
     +        event_number__high_bytes, event_number__low_bytes,
     +        crystals_matrix_amount_PHOS,
     +        crystal_matrix_type,
     +        amount_of_crystals_on_Z, amount_of_crystals_on_PHI,
     +        integer_radius,
     +        integer_crystal_size, integer_crystal_length

        do i=1,crystals_matrix_amount_PHOS
          integer_coordinate_Z   = matrix_coordinate_Z(i) /
     +                             matrix_coordinate_Z_unit
          integer_coordinate_PHI = matrix_coordinate_PHI(i) /
     +                             matrix_coordinate_PHI_unit
          write(event_file_unit_number,err=1),
     +          integer_coordinate_Z, integer_coordinate_PHI
        enddo

        n = 0
        do i=1,crystals_matrix_amount_PHOS
          crystals_amount_with_amplitudes(i) = 0
          do j=1,amount_of_crystals_on_Z*amount_of_crystals_on_PHI
            n = n + 1
c            write (*,*)n,crystals_amplitudes(n)
            if( crystals_amplitudes(n)/crystal_amplitudes_unit
     +          .GE. crystal_amplitudes_in_units_min ) then
              crystals_amount_with_amplitudes(i) = 
     +          crystals_amount_with_amplitudes(i) + 1
              crystals_amplitudes_Iad
     +                    (1,crystals_amount_with_amplitudes(i),i)=
     +          crystals_amplitudes(n)/crystal_amplitudes_unit
              crystals_amplitudes_Iad
     +                    (2,crystals_amount_with_amplitudes(i),i)=j-1
            endif
          enddo
        enddo

        do i=1,crystals_matrix_amount_PHOS
          write(event_file_unit_number,err=1),
     +          crystals_amount_with_amplitudes(i)
          write(event_file_unit_number,err=1),
     +         ((crystals_amplitudes_Iad(n,j,i),n=1,2),
     +                        j=1,crystals_amount_with_amplitudes(i))
        enddo

c        write (*,*)'Saving - DONE'
        return

1       stop 'Error in writing event data file'

      end

*****************************************************************************

      subroutine read_event_from_file
! at end of file or if there is error - event_number=-1
        implicit none
        include 'event_format.inc'

        integer n,i,j
        integer*2 integer_radius,
     +            integer_crystal_size, integer_crystal_length,
     +            integer_coordinate_Z, integer_coordinate_PHI,
     +            event_number__low_bytes, event_number__high_bytes

c        write (*,*)'Event reading...'

        read(event_file_unit_number,err=1,end=2),
     +        event_number__high_bytes, event_number__low_bytes,
     +        crystals_matrix_amount_PHOS,
     +        crystal_matrix_type,
     +        amount_of_crystals_on_Z, amount_of_crystals_on_PHI,
     +        integer_radius,
     +        integer_crystal_size, integer_crystal_length

        event_number = event_number__high_bytes
        event_number = event_number*32767 + event_number__low_bytes

        radius          = integer_radius          * radius_unit
        crystal_size    = integer_crystal_size    * crystal_size_unit
        crystal_length  = integer_crystal_length  * crystal_length_unit

        do i=1,crystals_matrix_amount_PHOS
          read(event_file_unit_number,err=1),
     +          integer_coordinate_Z, integer_coordinate_PHI
          matrix_coordinate_Z(i)   = integer_coordinate_Z   *
     +                               matrix_coordinate_Z_unit
          matrix_coordinate_PHI(i) = integer_coordinate_PHI *
     +                               matrix_coordinate_PHI_unit
        enddo

        do i=1,crystals_matrix_amount_PHOS
          read(event_file_unit_number,err=1),
     +          crystals_amount_with_amplitudes(i)
          read(event_file_unit_number,err=1),
     +         ((crystals_amplitudes_Iad(n,j,i),n=1,2),j=1,
     +                        crystals_amount_with_amplitudes(i))
        enddo

c        write (*,*)'Reading - DONE'
        return


2       event_number=-1
        return
1       stop 'Error in reading from data file.'

      end

*****************************************************************************

      subroutine open_event_file(unit_number,file_name,open_status)
* open_status  = 1   - new
* open_status <> 1   - old
        implicit none
        include 'event_format.inc'
        integer unit_number, open_status, flag
        character*(*) file_name

c        write (*,*)file_name

        event_number=0
        event_file_unit_number = -1
        if( open_status.EQ.1 ) then
          open(unit=unit_number,form='unformatted',name=file_name,
     +         status='new',err=1,iostat=flag)
        else
          open(unit=unit_number,form='unformatted',name=file_name,
     +         status='old',err=1,iostat=flag)
        endif
        event_file_unit_number = unit_number
        return
1       write (*,*)'Error in openning event data file. flag=',flag
        stop
      end

*****************************************************************************

      subroutine close_event_file
        implicit none
        include 'event_format.inc'
        if( event_file_unit_number.NE.-1 ) then
          close(event_file_unit_number)
        endif
      end

*****************************************************************************
