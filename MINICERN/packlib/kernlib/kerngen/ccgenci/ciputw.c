/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:54  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:35:20  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:47  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:19  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIPUTW
  CERN PROGLIB# Z311    CIPUTW          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CIPUTW (LUNDES, MBUF, NWPUT, ISTAT)
      write to the file :
       LUNDES  file descriptor
       MBUF    vector to be written
       NWPUT   number of full words to be written
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
void type_of_call ciputw_(lundes, mbuf, nwput, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciputw(lundes, mbuf, nwput, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CIPUTW(lundes, mbuf, nwput, stat)
#endif
      int  *mbuf;
      int  *lundes, *nwput, *stat;
{
      int  fildes;
      int  nbdn, nbdo;

      *stat = 0;
      if (*nwput <= 0)            return;

/*        write the file     */

      fildes = *lundes;
      nbdo   = *nwput * NBYTPW;
      nbdn   = write (fildes, mbuf, nbdo);
      if (nbdn < 0)               goto trouble;
      return;

#if defined(CERNLIB_PROJSHIFT)
trouble:  *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
trouble:  *stat = errno;
#endif
          perror (" error in CIPUTW");
          return;
}
/*> END <----------------------------------------------------------*/
