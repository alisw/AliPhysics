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
      SUBROUTINE GSTSET(FNAME,ANAME,INST,SITE,DEPT,RESP)
*
*
*      This subroutine starts a SET file.
*
*      It writes the beginning of the file, with a relevant header
*      A SET data assembly is intialized after the header.
*
#include "geant321/gcsetf.inc"
#include "geant321/gcunit.inc"
 
*
* *** name of the SET-file:
      CHARACTER*24    FNAME
      CHARACTER*80    FNAMES
* *** name of the data assembly:
      CHARACTER*4    ANAME
      CHARACTER*8     DATE,TIME
      CHARACTER*24    INST,SITE,DEPT,RESP
*
* *** File name and other general input to be read in by the user **
      IF(INDEX(FNAME,'.').EQ.0) THEN
         IT=LNBLNK(FNAME)
      ELSE
         IT=INDEX(FNAME,'.')-1
      ENDIF
#if !defined(CERNLIB_IBM)
      FNAMES=FNAME(1:IT)//'.set'
#endif
#if defined(CERNLIB_IBM)
      FNAMES='/'//FNAME(1:MIN(IT,8))//' SET A1'
      CALL CLTOU(FNAMES)
#endif
*
      OPEN (UNIT=NUNIT1, FILE=FNAMES, STATUS='UNKNOWN', ERR=10 )
*
* *** library call for date and time:
      CALL DATIMH (DATE,TIME)
*
      NPLACE = 1
      BLKSTR = ' '
      JVSTR  = ' '
      WRITE (BLKSTR,10000) INST,SITE,DEPT,RESP
      CALL GJWRIT
*   the last line:
      WRITE (NUNIT1,'(A80)') JVSTR
      JVSTR  = ' '
*
      WRITE (NUNIT1,10100) DATE(7:8),DATE(4:5),DATE(1:2),
     +                     TIME(1:2),TIME(4:5),TIME(7:8)
*
      NPLACE = 1
      BLKSTR = ' '
      WRITE (BLKSTR,10200) ANAME
      CALL GJWRIT
*   the last line:
      WRITE(NUNIT1,'(A80)')JVSTR
      JVSTR = ' '
*
      WRITE (NUNIT1,'(A)') ',1,0.01,0,-300.,300.,-300.,300.,-300.,300.'
      WRITE (NUNIT1,'(A)') '#9920,:12,0.005,:14,0,:15,0'
      WRITE (NUNIT1,'(A)') '@301,2#301,1.,0.,0.,0.,1.
     + ,0.,0.,0.,1.,0.,0.,0.'
*
      RETURN
*
   10 WRITE(CHMAIL,'('' Error opening the .SET file. '')')
      CALL GMAIL(0,0)
      WRITE(CHMAIL,'(''Check your directory and unit number in use'')')
      CALL GMAIL(0,0)
      WRITE(CHMAIL,'('' CADINT cannot proceed. Program terminated '')')
      CALL GMAIL(0,0)
      STOP
*
10000 FORMAT('@9900#9900,''Z68-300.89.06'',''2.0'',''
     +CADINT'',''1.20'',''',A,''',''',A,''',''',A,
     +''',''',A,'''')
10100 FORMAT('#9901,19',A2,',',A2,',',A2,',',A2,',',A2,',',A2)
10200 FORMAT('@9901,1,''',A4,'''#9905,3,1.')
*
      END
