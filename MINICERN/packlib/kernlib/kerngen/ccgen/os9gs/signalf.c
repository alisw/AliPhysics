/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNOS9  1.01  940722
  ORIG. 12/03/91, JZ
  FORTRAN interface routine to signal

      INTEGER FUNCTION SIGNALF (NUMSIG,PROC,IFLAG)

C-        NUMSIG :  signal number
C-          PROC :  external of the handler, if IFLAG = -1
C-         IFLAG :  < 0  instal PROC
C-                  = 0  default action
C-                  = 1  ignore signal
C-                  > 1  adr of handler as returned earlier
C-        function value = adr of previous handler
*/
#include <signal.h>
#if defined(CERNLIB_QX_SC)
int signalf_(signum,funct,flag)
#endif
#if defined(CERNLIB_QXNO_SC)
int signalf(signum,funct,flag)
#endif
#if defined(CERNLIB_QXCAPT)
int SIGNALF(signum,funct,flag)
#endif
      int  *signum, *flag;
      int  *funct;
{
      int  signo, istat;
      int  handler;
      void *oldhand;

      signo = *signum;

#if defined(CERNLIB_QCCINDAD)
      if (*flag < 0)          handler = *funct;
#endif
#if !defined(CERNLIB_QCCINDAD)
      if (*flag < 0)          handler = (int)funct;
#endif
        else if (*flag == 0)  handler = (int)SIG_DFL;
        else if (*flag == 1)  handler = (int)SIG_IGN;
        else                  handler = *flag;

      oldhand = signal(signo,(void *)handler);
      istat   = (int)oldhand;
#ifndef __GNUC__
      if (oldhand == SIG_ERR)  istat = -1;
#endif
      return istat;
}
/*> END <----------------------------------------------------------*/
