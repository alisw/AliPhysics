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
      SUBROUTINE GTXSET(FNAME,ANAME,NBINS,LUNIT,LUNIT2,
     +INST,SITE,DEPT,RESP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Writes out a description of the detector in SET file format.  *
C.    *  This is the main routine in the GEANT-SET interface           *
C.    *                                                                *
C.    *                    CADINT 1.20                                 *
C.    *                    ***********                                 *
C.    *    In this version VISIBLE VOLUMES are written out in a        *
C.    *    flat assembly.                                              *
C.    *    For questions, contact the authors.                         *
C.    *                                                                *
C.    *       Authors : NJ Hoimyr, J Vuoskoski     *********           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcsetf.inc"
 
#include "geant321/gcdlin.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
*
      CHARACTER*24    INST,SITE,DEPT,RESP
      CHARACTER*24    FNAME
      CHARACTER*4     ANAME,NAME
      INTEGER         LUNIT,LUNIT2,NBINS
*     ------------------------------------------------------------
*
      WRITE(CHMAIL,'('' Starting to write SET file...'')')
      CALL GMAIL(1,0)
*
*     local zebra links
      CALL MZLINT(IXSTOR,'/GCDLIN/',NCWORD,JCDFRS,JCDLST)
*
*     number of volumes:
      NVOL=0
*     counter for SET blocks:
      N1=2
*     Assign logical unit for .set file:
      NUNIT1=LUNIT
*     Assign logical unit for .mat file:
      NUNIT2=LUNIT2
*
* *** I n i t i a l i z e  S E T f i l e
*
*     Starting blocks:
      CALL GSTSET(FNAME,ANAME,INST,SITE,DEPT,RESP)
      NPLACE=1
*
*     Starts material listing file:
      CALL GSTMAT(FNAME)
*
*     c a l q.  n u m. o f  v o l u m e s
   10 IF (LQ(JVOLUM-NVOL-1).EQ.0.OR.IQ(JVOLUM-2).EQ.NVOL) GOTO 20
      NVOL=NVOL+1
      GOTO 10
   20 CONTINUE
*
*     Is NAME an existing volume ?
      NAME = ANAME
      CALL GLOOK (NAME,IQ(JVOLUM+1),NVOLUM,IVO)
*
* *** Uses information given by the drawing package **
*
      CALL GDEXCA (NAME,NBINS)
*
* *** e n d  o f  S E T  f i l e
*
* *** write the tree
*
      CALL GWRTRE (NAME, NVOL)
      CALL GENSET
*
      WRITE(CHMAIL,'('' '')')
      CALL GMAIL(0,0)
      WRITE(CHMAIL,'('' The selected geometry is now converted '')')
      CALL GMAIL(0,0)
      WRITE(CHMAIL,'('' into the SET file format and the file '')')
      CALL GMAIL(0,0)
      WRITE(CHMAIL,'('' is in your current working directory '')')
      CALL GMAIL(0,0)
*
      IF(JCADNT.NE.0) THEN
         CALL MZDROP(IXSTOR, JCADNT, ' ')
      ENDIF
      NCWORD(1) = 0
*
      END
