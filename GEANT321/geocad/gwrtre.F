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
      SUBROUTINE GWRTRE (VLNAME, NVOL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to write out the GEANT tree in the .mat file
C.    *       starting from the given volume                           *
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GTXSET
C.    *
C.    *       Author: Jouko Vuoskoski                                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdlin.inc"
#include "geant321/gcsetf.inc"
C
      CHARACTER*4 VLNAME,VOLUNA
      CHARACTER*80   HELPST
C
C    this has to be changed:
      IF(JCADNT.EQ.0) THEN
         CALL MZBOOK(IXSTOR,JCADNT,JCADNT,1,'CADI',1,1,0,2,-1)
         CALL MZBOOK(IXSTOR,JBUF1,
     +               JCADNT,-1,'CAD1',0,0,NVOLUM,2,-1)
      ENDIF
      DO 10 JV=1,NVOLUM
         IQ(JBUF1+JV)=0
   10 CONTINUE
C
C     Load IVO numbers of this particular part of the tree
C     First the volume where to be started
C
      DO 20 IVO=1, NVOL
         JVVOLU=IQ(JVOLUM+IVO)
         CALL UHTOC(JVVOLU,4,VOLUNA,4)
         IF(VOLUNA.EQ.VLNAME) THEN
            IQ(JBUF1+1)=IVO
         ENDIF
   20 CONTINUE
C
C      Then all the others
C
      JV=1
      DO 70 NH1=1, 15
         DO 60 NH2=1, JV
            JVO=LQ(JVOLUM-IQ(JBUF1+NH2))
            NIN=Q(JVO+3)
            IF(NIN.GT.0) THEN
               DO 40 IIN=1, NIN
                  JIN=LQ(JVO-IIN)
                  JVFLAG=0
                  DO 30 NH3=1, JV
                     IF(Q(JIN+2).EQ.IQ(JBUF1+NH3)) JVFLAG=1
   30             CONTINUE
                  IF(JVFLAG.EQ.0) THEN
                     JV=JV+1
                     IQ(JBUF1+JV)=Q(JIN+2)
                  ENDIF
   40          CONTINUE
            ELSEIF(NIN.LT.0) THEN
               JDIV=LQ(JVO-1)
               JVFLAG=0
               DO 50 NH3=1, JV
                  IF(Q(JDIV+2).EQ.IQ(JBUF1+NH3)) JVFLAG=1
   50          CONTINUE
               IF(JVFLAG.EQ.0) THEN
                  JV=JV+1
                  IQ(JBUF1+JV)=Q(JDIV+2)
               ENDIF
            ENDIF
   60    CONTINUE
   70 CONTINUE
C
C
C     write out the volumes into the .mat file
C
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)'GEANT TREE'
      WRITE (NUNIT2,*)'----------'
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)'The GEANT tree starting from the given volume'
      WRITE (NUNIT2,*)
C
      DO 90 NH1=1, JV
         JVO=LQ(JVOLUM-IQ(JBUF1+NH1))
         NIN=Q(JVO+3)
         IF(NIN.GT.0) THEN
            I1=10
            I2=16
            WRITE (HELPST,'(A80)')' '
            WRITE (HELPST(1:4),10000)IQ(JVOLUM+IQ(JBUF1+NH1))
            WRITE (HELPST(6:10),10100)NIN
            DO 80 IIN=1, NIN
               JIN=LQ(JVO-IIN)
               IVO=Q(JIN+2)
               WRITE (HELPST(I1:I2),10200)IQ(JVOLUM+IVO)
               I1=I1+6
               I2=I2+6
               IF (I2.GE.73) THEN
                  WRITE(NUNIT2,'(A80)')HELPST
                  I1=10
                  I2=16
                  WRITE (HELPST,'(A80)')' '
               ENDIF
   80       CONTINUE
            WRITE(NUNIT2,'(A80)')HELPST
C
         ELSEIF(NIN.LT.0) THEN
            JDIV=LQ(JVO-1)
            IVO=Q(JDIV+2)
            NUMDIV=Q(JDIV+3)
            WRITE(NUNIT2,10300)IQ(JVOLUM+IQ(JBUF1+NH1)),-NUMDIV,
     +      Q(JDIV+ 4),Q(JDIV+5),IQ(JVOLUM+IVO)
 
         ENDIF
   90 CONTINUE
C
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)
      WRITE (NUNIT2,*)'  ------ end of file -------'
C
10000 FORMAT(A4)
10100 FORMAT(I4)
10200 FORMAT(2X,A4)
10300 FORMAT(A4,1X,I4,2X,E15.8,2X,E15.8,2X,A4)
C
      END
