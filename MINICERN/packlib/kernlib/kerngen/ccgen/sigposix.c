/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:42  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:42  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNFOR  4.36  930602
  ORIG. 24/05/93, JZ
  FORTRAN interface routine to sigaction    */
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
      int  sigaction();
      int  istat, signo;

      struct sigaction newbuf;
      struct sigaction oldbuf;
      void (*action)();

#if defined(CERNLIB_QMHPX)
      static  subsq = 0;
#endif

      signo = *signum;

      if        (*flag == -1)  action = funct;
        else if (*flag == 0)   action = SIG_DFL;
        else if (*flag == 1)   action = SIG_IGN;
        else                   action = (void (*)())*flag;

      newbuf.sa_handler = action;
      newbuf.sa_flags   = 0;
      sigemptyset(&newbuf.sa_mask);

#if defined(CERNLIB_QMHPX)
/*    on HP sigaction is not properly initialized
 *    if running with a Fortran main program,
 *    a call to signal will fix this                  */
      if (subsq == 0 )
      {   subsq = 7;
          istat = signal (signo, action);
       }
#endif

      istat = sigaction(signo,&newbuf,&oldbuf);
      if (istat == 0)        return (int)oldbuf.sa_handler;
      return -1;
}
/*> END <----------------------------------------------------------*/
