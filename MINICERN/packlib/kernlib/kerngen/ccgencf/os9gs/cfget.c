/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:37  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE CFGET
  CERN PROGLIB# Z310    CFGET           .VERSION KERNOS9  1.01  940729
  ORIG. 12/01/91, JZ
  MOD.  29/07/94, MM    (remove label retn: )
      CALL CFGET (LUNDES, MEDIUM, NWREC, NWTAK, MBUF, ISTAT)
      read from the file :
       LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
       NWREC   number of words record size
      *NWTAK*  number of words to be read / actually read
      *MBUF    vector to be read into
      *ISTAT   status, =zero if success
*/
#include <stdio.h>
#include <errno.h>
#define NBYTPW 4       /* Number of bytes per word */

void cfget_(lundes, medium, nwrec, nwtak, mbuf, stat)
      char *mbuf;
      int  *lundes, *medium, *nwrec, *nwtak, *stat;
{
      int  fildes;
      int  nbdn, nbdo;

      *stat = 0;
      if (*nwtak <= 0)            return;

/*        read the file      */

      fildes = *lundes;
      nbdo   = *nwrec * NBYTPW;
      nbdn   = read (fildes, mbuf, nbdo);
      if (nbdn == 0)               goto heof;
      if (nbdn <  0)               goto herror;
      *nwtak = (nbdn - 1) / NBYTPW + 1;
      return;

/*        Handle exceptions        */

heof:     *stat = -1;
          return;

herror:   *stat = errno;
          perror (" error in CFGET");
          return;
}
/*> END <----------------------------------------------------------*/
