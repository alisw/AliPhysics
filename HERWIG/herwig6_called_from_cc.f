
      subroutine herwig6_open_fortran_file(lun, name)
      implicit none
      integer lun
      character*(*) name

      open (lun, file=name)
      return
      end

      subroutine herwig6_close_fortran_file(lun)
      implicit none
      integer lun
      close (lun)
      return
      end


