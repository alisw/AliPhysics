/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:53  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:35:18  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1997/01/15 16:25:45  cernlib
 * fix from F.Hemmer to return rfio return code
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:17  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMOS9)
#include "os9gs/ciopei.c"
#else
/*>    ROUTINE CIOPEI
  CERN PROGLIB# Z311    CIOPEI          .VERSION KERNFOR  4.39  940228
  ORIG. 12/10/91, JZ
      CALL CIOPEN (LUNDES, MODE, TEXT, ISTAT)
      open a file :
      *LUNDES  file descriptor
       MODE    string selecting IO mode
               = 'r ', 'w ', 'a ', 'r+ ', ...
       TEXT    name of the file
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_open.h"
#ifndef WIN32
#  include <errno.h>
#else
#  include <stdlib.h>
#endif
#include "kerngen/cf_xaft.h"
#include "kerngen/fortchar.h"
      int ciopen_perm = 0;
#if defined(CERNLIB_QX_SC)
void type_of_call ciopei_(lundes,mode,ftext,stat,lgtx)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciopei(lundes,mode,ftext,stat,lgtx)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifdef CERNLIB_MSSTDCALL
   void type_of_call CIOPEI(lundes,mode,ftext,lftext,stat,lgtx)
   int lftext;
#  else
   void type_of_call CIOPEI(lundes,mode,ftext,stat,lgtx)
#  endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif
      int  *lundes, *stat, *lgtx;
      int  *mode;
{
      char *pttext, *fchtak();
      int  flags;
      int  fildes;
      int  perm;

      *lundes = 0;
      *stat   = -1;

      perm = ciopen_perm;
      ciopen_perm = 0;

/*        construct flags :
            mode[0] =    0 r    1 w    2 a
            mode[1] =    1 +
*/
/*        flags for disk     */


      if (mode[0] == 0)
        {if (mode[1] == 0)
          flags = O_RDONLY;
        else
          flags = O_RDWR;}

      else if (mode[0] == 1)
        {if (mode[1] == 0)
          flags = O_WRONLY | O_CREAT | O_TRUNC;
        else
          flags = O_RDWR | O_CREAT | O_TRUNC;}

      else if (mode[0] == 2)
        {if (mode[1] == 0)
          flags = O_WRONLY | O_CREAT | O_APPEND;
        else
          flags = O_RDWR | O_CREAT | O_APPEND;}

/*        open the file      */

      pttext = fchtak(ftext,*lgtx);
      if (pttext == 0)             return;

      if (perm == 0)   perm = 0644;

#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
      fildes = open (pttext, flags | O_BINARY, perm);
#else
      fildes = open (pttext, flags, perm);
#endif
      if (fildes < 0)              goto errm;
      *lundes = fildes;
      *stat   = 0;
      goto done;

#if defined(CERNLIB_PROJSHIFT)
errm: *stat = (serrno ? serrno : (rfio_errno ? rfio_errno : errno));
#else
errm: *stat = errno;
#endif
/*    perror (" error in CIOPEN");  */

done: free(pttext);
      return;
}
/*> END <----------------------------------------------------------*/
#endif
