/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:33  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNIRT  1.03  910314
  ORIG. 12/03/91, JZ
  FORTRAN interface routine to sigaction    */
#include <stdio.h>
#include <signal.h>
#include <errno.h>
#if defined(CERNLIB_IBMRT)&&defined(CERNLIB_QXNO_SC)
int signalf(signum,funct,flag)
#else
int signalf_(signum,funct,flag)
#endif
      long *signum, *flag;
      long  *funct;
{
      int  sigaction();
      int  istat, signo;

      struct mysig {
          int       sa_handler;
          sigset_t  sa_mask;
          int       sa_flags;
         };

      struct mysig newbuf;
      struct mysig oldbuf;

      signo = *signum;

      if        (*flag < 0)    newbuf.sa_handler = funct;
        else if (*flag == 0)   newbuf.sa_handler = SIG_DFL;
        else if (*flag == 1)   newbuf.sa_handler = SIG_IGN;
        else                   newbuf.sa_handler = *flag;

      newbuf.sa_flags   = 0;
      sigemptyset(&newbuf.sa_mask);

      istat = sigaction(signo,&newbuf,&oldbuf);
      if (istat == 0)        return oldbuf.sa_handler;
      return -errno;
}
/*> END <----------------------------------------------------------*/
