*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_QMUIX)
      SUBROUTINE RZSTRIP(CHFILE,NST,NRS,LREC,ISTAT)
*
**********************************************************************
*
*        Routine called by RZOPEN and RZIODO for striped files
*
*       CHFILE  Name of the header description file (Ascii)
*       NST*    Number of striped files
*       NRS*    Number of records in each file
*       LREC*   Record length in words of the files
*       ISTAT*  return sattus
*
*           Description of the /RZCSTR/ structure
*       ISLAST  File number currently described in the structure
*       ISTRIP  Current striped file number for file ISLAST
*       NSTRIP  Number of striped files for ISLAST
*       NRSTRIP Number of records in each file of file ISLAST
*
*            Example of a descriptor file (starts in column 1)
*            ============================
*
*          The first line after the ###FTRIPE line contains
*       nstrip, Total_number_of_records, Record_Length_in_words, File_size
*
*       ###FSTRIPE-V001###
*       5 3872 1024 15859712
*       775 3174400 /user/brun/hbook/hart.hbook_1
*       775 3174400 /user/brun/hbook/hart.hbook_2
*       775 3174400 /user/brun/hbook/hart.hbook_3
*       775 3174400 /user/brun/hbook/hart.hbook_4
*       772 3162112 /user/brun/hbook/hart.hbook_5
*
*
*  Author  : R.Brun CERN/CN
*  Written : 10.06.94
*
************************************************************************
*
#include "zebra/rzcstr.inc"
      character*(*) chfile
      character*128 file
      character*10 strip
      integer istat,nst,nrs
      parameter (lunhead=121)
*
*-----------------------------------------------------------------------
*
      istat=1
      open(unit=lunhead,file=chfile,status='old')
      read(lunhead,'(a)',err=20,end=20)strip
      if(strip.ne.'###FSTRIPE')go to 20
      read(lunhead,*,err=20,end=20)nst,nrtot,lrec,nwords
      if(nst.le.0)go to 20
      nrs=(nrtot+nst-1)/nst
      if(nrs.le.1)go to 20
      if(nst.ge.maxstrip)go to 20
      do 10 i=1,nst
         read(lunhead,'(a)',err=20,end=20)file
         isl=index(file,'/')
         if(isl.eq.0)go to 20
         rzsfile(i)=file(isl:)
  10  continue
      rzsfile(maxstrip)=chfile
      istat=0
  20  close(lunhead)
*
       end
#endif
