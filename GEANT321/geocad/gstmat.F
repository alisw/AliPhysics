*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GSTMAT(FNAME)
*
* *** Opens a material list file for CADD output
*
#include "geant321/gcsetf.inc"
*
* *** name of the SET-file:
      CHARACTER*24    FNAME
* *** name of the MAT-file:
      CHARACTER*80    FNAMEM
*
* *** Append filename FNAME to FNAMEM with extension .mat:
      IF(INDEX(FNAME,'.').EQ.0) THEN
         IT=LNBLNK(FNAME)
      ELSE
         IT=INDEX(FNAME,'.')-1
      ENDIF
#if !defined(CERNLIB_IBM)
      FNAMEM=FNAME(1:IT)//'.mat'
#endif
#if defined(CERNLIB_IBM)
      FNAMEM='/'//FNAME(1:MIN(IT,8))//' MAT A1'
      CALL CLTOU(FNAMEM)
#endif
*
* *** Open material file:
      OPEN  (UNIT=NUNIT2, FILE=FNAMEM, STATUS='UNKNOWN', ERR=10 )
*
      WRITE (NUNIT2,*)'GEANT-SET MATERIAL LISTING FILE'
      WRITE (NUNIT2,*)'--------------------------------'
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)'Materials in the geometry described in'
      WRITE (NUNIT2,*)'the .SET file:  ',FNAME
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,10000)
      WRITE (NUNIT2,*)
      RETURN
*
   10 WRITE(*,*)'Error opening the .mat file for materials'
      WRITE(*,*)'Check use of logical units.'
*
10000 FORMAT('Volume name    ','Tracking media ',9X
     + ,'Material  ',14X,'Density  ')
*
      END
