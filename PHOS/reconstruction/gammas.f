*****************************************************************************

      subroutine open_gamma_impulses_file
     +                     (unit_number,file_name,open_status)
* open_status  = 1   - new
* open_status <> 1   - old
        implicit none
        integer unit_number, open_status, flag
        character*(*) file_name

        if( open_status.EQ.1 ) then
          open(unit=unit_number,form='unformatted',name=file_name,
     +         status='new',err=1,iostat=flag)
        else
          open(unit=unit_number,form='unformatted',name=file_name,
     +         status='old',err=1,iostat=flag)
        endif
        return

1       write (*,*) 'Error in openning gammas file. flag=',flag
        stop
      end

*****************************************************************************

      subroutine close_gamma_impulses_file(unit_number)
        implicit none
        integer unit_number
        close(unit_number)
      end

*****************************************************************************

      subroutine save_gamma_impulse(unit_number,event_number,
     +                              gammas_amount,gamma_impulse)
        implicit none
        integer unit_number
        integer event_number
        integer*2 event_number__high_bytes, event_number__low_bytes

        integer gammas_amount, i,j
        real gamma_impulse(3,*)

        event_number__high_bytes = event_number/32767
        event_number__low_bytes  = event_number -
     +                             event_number__high_bytes*32767

        write(unit_number,err=1),
     +        event_number__high_bytes, event_number__low_bytes,
     +        gammas_amount
        write(unit_number,err=1),
     +        ((gamma_impulse(i,j),i=1,3),j=1,gammas_amount)

        return
1       stop 'Error in writing to gammas file.'
      end

*****************************************************************************

      subroutine read_gamma_impulse(unit_number,event_number,
     +                              gammas_amount,gamma_impulse)
        implicit none
        integer gammas_amount, i,j, unit_number
        integer event_number
        integer*2 event_number__high_bytes, event_number__low_bytes
        real gamma_impulse(3,*)

        read(unit_number,err=1,end=2)
     +        event_number__high_bytes, event_number__low_bytes,
     +        gammas_amount
        read(unit_number,err=1)
     +        ((gamma_impulse(i,j),i=1,3),j=1,gammas_amount)

        event_number = event_number__high_bytes
        event_number = event_number*32767 + event_number__low_bytes

        return

2       event_number = -1
        return

1       stop 'Error in reading from gammas file.'
      end

*****************************************************************************
