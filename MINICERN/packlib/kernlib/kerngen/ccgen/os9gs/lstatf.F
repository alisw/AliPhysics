*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
* Kernlib
*
*
*	T=XCC.
      INTEGER FUNCTION LSTATF (NAME, INFO)
C
C CERN PROGLIB# Z265    LSTATF          .VERSION KERNOS9  1.01  940906
C ORIG. 14/03/91, RDM
C MOD.  06/09/94, MM
C Fortran interface routine to lstat:
C     The lstat system call is not available on OS-9 system. OS-9 does not
C     support links. Therefore the call to LSTATI is replaced by a simple
C     call to STATI instead.  
C
C        NAME   the name of the file of which information is required
C     INFO(12)  buffer returning the information about "file"
C
C        A succesfull call returns 0. A failed call returns -1.
C
      CHARACTER*(*)  NAME
      INTEGER        INFO(12), STATI

      LGN   = LNBLNK (NAME)
#if defined(CERNLIB_QMVAX)
      STATF = STATI (%REF(NAME), INFO, LGN)
#endif
#if !defined(CERNLIB_QMVAX)
      STATF = STATI (NAME, INFO, LGN)
#endif
      RETURN
      END
