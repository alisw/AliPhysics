*======================================================================*
* Routines to open and close a given file with a given FORTRAN unit    *
* from the C++ world                                                   *
*======================================================================*
*
*=== FLUKA_OPENINP ====================================================*
*
      SUBROUTINE FLUKA_OPENINP(IOUNIT,FILNAM)
*
*----------------------------------------------------------------------*
* Opens a file with a given unit number
*
*
* IOUNIT: Input unit to be assiged to the file
* FILNAM: Name of the file
*
*----------------------------------------------------------------------*
*

      IMPLICIT NONE
      INTEGER IOUNIT
      CHARACTER*(*) FILNAM

      PRINT *, '==> FLUKA_OPENINP(',IOUNIT,',',FILNAM,')'

      OPEN (UNIT=IOUNIT, FILE=FILNAM, STATUS="OLD")

      PRINT *, '<== FLUKA_OPENINP(',IOUNIT,',',FILNAM,')'

      RETURN
 9999 END

      SUBROUTINE FLUKA_OPENOUT(IOUNIT,FILNAM)
*
*----------------------------------------------------------------------*
* Opens a file with a given unit number
*
*
* IOUNIT: Input unit to be assiged to the file
* FILNAM: Name of the file
*
*----------------------------------------------------------------------*
*

      IMPLICIT NONE
      INTEGER IOUNIT
      CHARACTER*(*) FILNAM

      PRINT *, '==> FLUKA_OPENOUT(',IOUNIT,',',FILNAM,')'

      OPEN (UNIT=IOUNIT, FILE=FILNAM, STATUS="UNKNOWN")

      PRINT *, '<== FLUKA_OPENOUT(',IOUNIT,',',FILNAM,')'

      RETURN
 9999 END


*
*=== FLUKA_CLOSEINP ====================================================*
*
      SUBROUTINE FLUKA_CLOSEINP(IOUNIT)
*
*----------------------------------------------------------------------*
* Closes the given unit number
*
*
* IOUNIT: Input unit to be assiged to the file
*
*----------------------------------------------------------------------*
*
      IMPLICIT NONE
      INTEGER IOUNIT

      PRINT *, '==> FLUKA_CLOSEINP(',IOUNIT,')'

      CLOSE (UNIT=IOUNIT)

      PRINT *, '<== FLUKA_CLOSEINP(',IOUNIT,')'
      RETURN
 9999 END
