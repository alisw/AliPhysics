/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:37  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE CFOPEI
  CERN PROGLIB# Z310    CFOPEI          .VERSION KERNOS9  1.01  940801
  ORIG. 12/01/91, JZ
      CALL CFOPEN (LUNDES, MEDIUM, NWREC, MODE, NBUF, TEXT, ISTAT)
      open a file :
      *LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
       NWREC   record length in number of words
       MODE    string selecting IO mode
               = 'r ', 'w ', 'a ', 'r+ ', ...
       NBUF    number of buffers to be allocated, (not used)
       TEXT    name of the file
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf#open.h"
#include <errno.h>
#include <modes.h>
#include "kerngen/cf#xaft.h"
#include "kerngen/fortchar.h"
#include "kerngen/wordsizc.h"
      int cfopen_perm = 0;

void cfopei_(lundes,medium,nwrec,mode,nbuf,ftext,stat,lgtx)
      char *ftext;
      int  *lundes, *medium, *nwrec, *nbuf, *stat, *lgtx;
      int  *mode;
{
      char *pttext, *fchtak();
      int  flags;
      int  fildes;
      int  perm;

      *lundes = 0;
      *stat   = -1;

      perm = cfopen_perm;
      cfopen_perm = 0;

/*
 *    construct flags :
 *      mode[0] =    0 r    1 w    2 a
 *      mode[1] =    1 +
 */

      if ((*medium == 1) || (*medium == 3))
      {

/*
 *    flags for tape
 */

          if (mode[0] == 0)
          {
              if (mode[1] == 0)
                  flags = FAM_READ;
              else
                  flags = FAM_READ | FAM_WRITE;

          } else if (mode[0] == 1) {

              if (mode[1] == 0)
                  flags = FAM_WRITE;
              else
                  flags = FAM_READ | FAM_WRITE;

          } else if (mode[0] == 2)       return;

      } else {

/*
 *    flags for disk
 */

          if (mode[0] == 0)
          {
              if (mode[1] == 0)
                  flags = FAM_READ;
              else
                  flags = FAM_READ | FAM_WRITE;

          } else if (mode[0] == 1) {

              if (mode[1] == 0)
                  flags = FAM_WRITE;
              else
                  flags = FAM_WRITE | FAM_READ;

          } else if (mode[0] == 2) {

              if (mode[1] == 0)
                  flags = FAM_WRITE | FAM_APPEND;
              else
                  flags = FAM_WRITE | FAM_READ | FAM_APPEND;
          }
      }

/*
 *    open the file
 */

      pttext = fchtak(ftext,*lgtx);
      if (pttext == 0) return;

      if (perm == 0)   perm = FAP_READ | FAP_WRITE | FAP_PREAD;

      if ( (mode[0] == 1) &
           ((*medium == 0) || (*medium == 2))
         ) {
              if ( (fildes = create (pttext, flags, perm)) < 0 )
                  fildes = creat (pttext, flags);

      } else {
          fildes = open (pttext, flags);
          if ((mode[0] == 2) &
              (fildes < 0) ) fildes = create (pttext, flags, perm);
      }

      if (fildes < 0)  goto errm;

      *lundes = fildes;
      *stat   = 0;
      goto done;

errm: *stat = errno;
      perror (" error in CFOPEN");

done: free(pttext);
      return;
}
/*> END <----------------------------------------------------------*/
