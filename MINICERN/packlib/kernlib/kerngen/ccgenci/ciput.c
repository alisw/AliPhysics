/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:54  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:35:19  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:46  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:18  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIPUT
  CERN PROGLIB# Z311    CIPUT           .VERSION KERNFOR  4.37  930715
  ORIG. 12/10/91, JZ
      CALL CIPUT (LUNDES, MBUF, NBPUT, ISTAT)
      write to the file :
       LUNDES  file descriptor
       MBUF    vector to be written
       NBPUT   number of bytes to be written
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_reaw.h"
#ifndef WIN32
#  include <errno.h>
#else
#  include <stdlib.h>
#endif
#include "kerngen/cf_xaft.h"
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QMVAX)
#include <descrip.h>
#endif
#if defined(CERNLIB_QX_SC)
void type_of_call ciput_(lundes, mbuf, nbput, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciput(lundes, mbuf, nbput, stat)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifdef CERNLIB_MSSTDCALL
     void type_of_call CIPUT(lundes, mbuf, lmbuf, nbput, stat)
     int lmbuf;
#  else
     void type_of_call CIPUT(lundes, mbuf, nbput, stat)
#  endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd mbuf;
#endif
#if defined(CERNLIB_QMVAX)
      struct dsc$descriptor_s  *mbuf;
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMVAX))
      char *mbuf;
#endif
      int  *lundes, *nbput, *stat;
{
      char *ubuf;
      int  fildes;
      int  nbdn, nbdo;

      *stat = 0;
      if (*nbput <= 0)            return;
#if defined(CERNLIB_QMCRY)
      ubuf = _fcdtocp(mbuf);
#endif
#if defined(CERNLIB_QMVAX)
      ubuf = mbuf->dsc$a_pointer;
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMVAX))
      ubuf = mbuf;
#endif

/*        write the file     */

      fildes = *lundes;
      nbdo   = *nbput;
      nbdn   = write (fildes, ubuf, nbdo);
      if (nbdn < 0)               goto trouble;
      return;

#if defined(CERNLIB_PROJSHIFT)
trouble:  *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
trouble:  *stat = errno;
#endif
          perror (" error in CIPUT");
          return;
}
/*> END <----------------------------------------------------------*/
