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
      SUBROUTINE RZSTRIR(LUNRZ,IREC)
*
**********************************************************************
*
*        Routine called by RZIODO for striped files
*       This routine checks if the input value for IREC is on the
*       current striped file. If not, it closes the current file
*       and open the file corresponding to the input record IREC.
*       If LUNRZ does not correspond to ISLAST, the routine opens
*       the descriptor text file corresponding to LUNPTR and loads
*       the commons /RZCSTR./ with the corresponding file names
*       and attributes
*
*       IREC*   Record number in current striped file
*
*  Author  : R.Brun CERN/CN
*  Written : 10.06.94
*
************************************************************************
*
#include "zebra/rzclun.inc"
#include "zebra/rzcstr.inc"
      integer istat,nst,nrs,irec
*
*-----------------------------------------------------------------------
*
      lunc=lunrz
      nrs=nrstrip(lunc)
      irlow =(istrip(lunc)-1)*nrs+1
      irnext=irlow+nrs
      if(irec.ge.irlow.and.irec.lt.irnext)then
         irec=irec-irlow+1
      else
         if(imodec.eq.0)then
            close(lunc)
         else
            call cfclos(lunc,0)
         endif
         is=(irec+nrs-1)/nrs
         istrip(lunc)=is
         irec=irec-(is-1)*nrs
         if(islast.ne.lunc)then
            call rzstrip(rznames(lunc),nst,nrs,lrec,istat)
         endif
         if(imodec.ne.0)then
            call cfopen(lunc,0,izrecl,'r',0,rzsfile(is),istat)
         else
            OPEN(UNIT=LUNC,FILE=rzsfile(is),FORM='UNFORMATTED',
#if defined(CERNLIB_RZBYTES)
     +         RECL=4*lrec,ACCESS='DIRECT',STATUS='old',IOSTAT=ISTAT)
#endif
#if !defined(CERNLIB_RZBYTES)
     +         RECL=lrec,ACCESS='DIRECT',STATUS='old',IOSTAT=ISTAT)
#endif
         endif
         islast=lunc
      endif
*
      end
#endif
