*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE MATISO(IZ,IA,NI,IDISO,FSINGL,NUNIT)
C
C Search array MATIDS for the isotopes which have to be taken
C into account for the element described by IZ and IA
C
#include "geant321/mmicap.inc"
#include "geant321/minput.inc"
C
      DIMENSION IDISO(20,2)
      LOGICAL FSINGL
C
      IF(IZ.GT.0.AND.IZ.LE.100.and.MATIDS(IZ,1,1).GT.0) THEN
        ID   = IZ*1000+IA
        IF = 0
        IC = 0
        IDIFF = 1000000
C
C check first if selected isotope available
        DO 10 I=2,MATIDS(IZ,1,1)+1
          IF( MATIDS(IZ,I,1).EQ.ID  .AND.
     +       (MATIDS(IZ,I,2).EQ.100 .OR. FSINGL)) IF = I
          IF( IABS(MATIDS(IZ,I,1)-ID).LT.IDIFF) THEN
             IDIFF = IABS(MATIDS(IZ,I,1)-ID)
             IC = I
          ENDIF
          IDISO(I-1,1) = MATIDS(IZ,I,1)
          IDISO(I-1,2) = MATIDS(IZ,I,2)
   10   CONTINUE
        NI = 1
C the unit number on which the x-section is stored
        NUNIT = MATIDS(IZ,1,2)
        IF(.NOT. FSINGL) THEN
          IF(IF .EQ. 0) THEN
C no matching isotope found. Look for closest one
             IF(MATIDS(IZ,2,2).NE.100) NI = MATIDS(IZ,1,1)
          ELSE
C matching isotope found
            IDISO(1,1) = MATIDS(IZ,IF,1)
            IDISO(1,2) = 100
          ENDIF
        ELSE
          IDISO(1,1) = MATIDS(IZ,IC,1)
          IDISO(1,2) = 100
        ENDIF
      ELSE
         WRITE(IOUT,'('' MATISO: Error in neutron x-section '',
     +                ''file detected - Z = '',I4)') IZ
         WRITE(6,'('' MICAP : Error in x-section file '',
     +             '' detected -> STOP '')')
         STOP
      ENDIF
      RETURN
      END
