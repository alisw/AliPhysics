/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:34  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNSGI  1.04  930120
  ORIG. 12/03/91, JZ
  FORTRAN interface routine to sigaction    */
#include <signal.h>
#include <stdio.h>
#include <errno.h>
int signalf_(signum,funct,flag)
      int  *signum, *flag;
      void (*funct)();
{
      int  sigaction();
      int  istat, signo;

      struct sigaction newbuf;
      struct sigaction oldbuf;

      signo = *signum;

      if        (*flag < 0)    newbuf.sa_handler = funct;
        else if (*flag == 0)   newbuf.sa_handler = SIG_DFL;
        else if (*flag == 1)   newbuf.sa_handler = SIG_IGN;
        else                   newbuf.sa_handler = (void (*)())*flag;

      newbuf.sa_flags   = 0;
      sigemptyset(&newbuf.sa_mask);

      istat = sigaction(signo,&newbuf,&oldbuf);
      if (istat == 0)        return (int)oldbuf.sa_handler;
      return -1;
}
/*> END <----------------------------------------------------------*/
