/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:41  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:41  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMSGI)
#include "sgigs/signalf.c"
#elif defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIRTD)
#include "irtgs/signalf.c"
#elif defined(CERNLIB_QMOS9)
#include "os9gs/signalf.c"
#elif defined(CERNLIB_QSIGBSD)
#include "sigbsd.c"
#elif defined(CERNLIB_QSIGPOSIX)
#include "sigposix.c"
#else
/*>    ROUTINE SIGNALF
  CERN PROGLIB#         SIGNALF         .VERSION KERNFOR  4.38  931108
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
int type_of_call signalf_(signum,funct,flag)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call signalf(signum,funct,flag)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call SIGNALF(signum,funct,flag)
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

      oldhand = signal(signo,handler);
      istat   = (int)oldhand;
#ifndef __GNUC__
      if (oldhand == SIG_ERR)  istat = -1;
#endif
      return istat;
}
/*> END <----------------------------------------------------------*/
#endif
