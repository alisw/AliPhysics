/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:40  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:41  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNFOR  4.36  930602
  ORIG. 24/05/93, JZ
  FORTRAN interface routine to sigvec    */
#include <signal.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
int type_of_call signalf_(signum,funct,flag)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call signalf(signum,funct,flag)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call SIGNALF(signum,funct,flag)
#endif
      int  *signum, *flag;
      void (*funct)();
{
      int  istat, signo;

      struct sigvec newbuf;
      struct sigvec oldbuf;

      signo = *signum;

      if        (*flag < 0)    newbuf.sv_handler = funct;
        else if (*flag == 0)   newbuf.sv_handler = SIG_DFL;
        else if (*flag == 1)   newbuf.sv_handler = SIG_IGN;
        else                   newbuf.sv_handler = (void (*)())*flag;

      newbuf.sv_flags = 0;
      newbuf.sv_mask  = 0;

      istat = sigvec(signo,&newbuf,&oldbuf);
      if (istat == 0)        return (int)oldbuf.sv_handler;
      return -1;
}
/*> END <----------------------------------------------------------*/
