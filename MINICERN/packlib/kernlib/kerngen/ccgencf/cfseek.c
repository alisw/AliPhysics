/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:13  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:12  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE CFSEEK
  CERN PROGLIB# Z310    CFSEEK          .VERSION KERNFOR  4.29  910718
  ORIG. 12/01/91, JZ
      CALL CFSEEK (LUNDES, MEDIUM, NWREC, JCREC, ISTAT)
      reposition the file :
       LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
       NWREC   number of words per record
       JCREC   number of records before current
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"
#include "kerngen/wordsizc.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call cfseek_(lundes, medium, nwrec, jcrec, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cfseek(lundes, medium, nwrec, jcrec, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CFSEEK(lundes, medium, nwrec, jcrec, stat)
#endif
      int  *lundes, *medium, *nwrec, *jcrec, *stat;
{
      int  fildes;
      int  nbdo;
      int  isw;

/*        position the file        */

      fildes = *lundes;
      nbdo   = *jcrec * *nwrec * NBYTPW;
      isw    = lseek (fildes, nbdo, 0);
      if (isw <  0)                goto trouble;
      *stat = 0;
      return;

trouble:  *stat = -1;
          perror (" error in CFSEEK");
          return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_CFSEEK
#undef CERNLIB_TCGEN_CFSEEK
#endif
