c-------------------------------------------------------------------------------
c  Jul 02 1998 P.Murat: routines to be called from C++ side
c-------------------------------------------------------------------------------
      subroutine tpythia6_open_fortran_file(lun, name)
      implicit none
      integer lun
      character*(*) name

      open (lun, file=name)
      return
      end

      subroutine tpythia6_close_fortran_file(lun)
      implicit none
      integer lun
      close (lun)
      return
      end


