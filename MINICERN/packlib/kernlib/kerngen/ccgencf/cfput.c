/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/10/23 16:33:19  mclareni
 * NT mods
 *
 * Revision 1.3  1997/02/04 17:35:12  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:33  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:11  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CFPUT
  CERN PROGLIB# Z310    CFPUT           .VERSION KERNFOR  4.29  910718
  ORIG. 12/01/91, JZ
      CALL CFPUT (LUNDES, MEDIUM, NWREC, MBUF, ISTAT)
      write to the file :
       LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
       NWREC   record size, number of words to be written
       MBUF    vector to be written
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_reaw.h"
#include <errno.h>
#include "kerngen/cf_xaft.h"
#include "kerngen/wordsizc.h"
#if defined(CERNLIB_QX_SC)
void type_of_call cfput_(lundes, medium, nwrec, mbuf, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cfput(lundes, medium, nwrec, mbuf, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CFPUT(lundes, medium, nwrec, mbuf, 
# ifdef CERNLIB_CFPUT_CHARACTER
     lmbuf,
# endif
            stat)
#  ifdef CERNLIB_CFPUT_CHARACTER
     int lmbuf;
#  endif
#endif
      char *mbuf;
      int  *lundes, *medium, *nwrec, *stat;
{
      int  fildes;
      int  nbdn, nbdo;

      *stat = 0;
      if (*nwrec <= 0)            return;

/*        write the file     */

      fildes = *lundes;
      nbdo   = *nwrec * NBYTPW;
      nbdn   = write (fildes, mbuf, nbdo);
      if (nbdn < 0)               goto trouble;
      return;

#if defined(CERNLIB_PROJSHIFT)
trouble:  *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
trouble:  *stat = errno;
#endif
          perror (" error in CFPUT");
          return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_CFPUT
#undef CERNLIB_TCGEN_CFPUT
#endif
