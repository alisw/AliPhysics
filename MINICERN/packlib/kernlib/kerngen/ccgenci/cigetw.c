/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:50  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:35:17  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:44  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:17  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:37  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIGETW
  CERN PROGLIB# Z311    CIGETW          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CIGETW (LUNDES, MBUF, NWDO, NWDONE, ISTAT)
      read from the file :
       LUNDES  file descriptor
      *MBUF    vector to be read into
       NWDO    number of full words to be read
      *NWDONE  number of full words actually read
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_reaw.h"
#ifndef WIN32
#  include <errno.h>
#else
#  include <stdlib.h>
#endif
#include "kerngen/cf_xaft.h"
#include "kerngen/wordsizc.h"
#if defined(CERNLIB_QX_SC)
void type_of_call cigetw_(lundes, mbuf, nwdo, nwdone, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cigetw(lundes, mbuf, nwdo, nwdone, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CIGETW(lundes, mbuf, nwdo, nwdone, stat)
#endif
      int  *mbuf;
      int  *lundes, *nwdo, *nwdone, *stat;
{
      int  fildes;
      int  nbdn, nbxq;

      *stat = 0;
      if (*nwdo <= 0)            return;

/*        read the file      */

      fildes = *lundes;
      nbxq   = *nwdo * NBYTPW;
      nbdn   = read (fildes, mbuf, nbxq);
      if (nbdn == 0)               goto heof;
      if (nbdn <  0)               goto herror;
      *nwdone = nbdn / NBYTPW;
      return;

/*        Handle exceptions        */

heof:     *stat = -1;
          return;

#if defined(CERNLIB_PROJSHIFT)
herror:   *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
herror:   *stat = errno;
#endif
          perror (" error in CIGETW");
          return;
}
/*> END <----------------------------------------------------------*/
