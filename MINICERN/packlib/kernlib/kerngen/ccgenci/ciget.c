/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:49  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:35:17  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:44  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:16  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:37  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIGET
  CERN PROGLIB# Z311    CIGET           .VERSION KERNFOR  4.37  930715
  ORIG. 12/10/91, JZ
      CALL CIGET (LUNDES, MBUF, NBDO, NBDONE, ISTAT)
      read from the file :
       LUNDES  file descriptor
      *MBUF    vector to be read into
       NBDO    number of bytes to be read
      *NBDONE  number of bytes actually read
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
void type_of_call ciget_(lundes, mbuf, nbdo, nbdone, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciget(lundes, mbuf, nbdo, nbdone, stat)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifdef CERNLIB_MSSTDCALL
     void type_of_call CIGET(lundes, mbuf, lmbuf, nbdo, nbdone, stat)
     int lmbuf;
#  else
     void type_of_call CIGET(lundes, mbuf, nbdo, nbdone, stat)
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
      int  *lundes, *nbdo, *nbdone, *stat;
{
      char *ubuf;
      int  fildes;
      int  nbdn, nbxq;

      *stat = 0;
      if (*nbdo <= 0)            return;
#if defined(CERNLIB_QMCRY)
      ubuf = _fcdtocp(mbuf);
#endif
#if defined(CERNLIB_QMVAX)
      ubuf = mbuf->dsc$a_pointer;
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMVAX))
      ubuf = mbuf;
#endif

/*        read the file      */

      fildes = *lundes;
      nbxq   = *nbdo;
      nbdn   = read (fildes, ubuf, nbxq);
      if (nbdn == 0)               goto heof;
      if (nbdn <  0)               goto herror;
      *nbdone = nbdn;
      return;

/*        Handle exceptions        */

heof:     *stat = -1;
          return;

#if defined(CERNLIB_PROJSHIFT)
herror:   *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
herror:   *stat = errno;
#endif
          perror (" error in CIGET");
          return;
}
/*> END <----------------------------------------------------------*/
