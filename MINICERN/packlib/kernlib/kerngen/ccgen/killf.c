/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:37  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:24  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:32  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:24  mclareni
 * Kernlib
 *
 */
#ifdef WIN32
# include <windows.h>
#endif

#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"


/*>    ROUTINE KILLF (IPID,ISIG)
  CERN PROGLIB# Z265    KILLF           .VERSION KERNFOR  4.26  910313
  ORIG. 22/02/91, JZ
  Fortran interface routine to kill
*/
#if defined(CERNLIB_QX_SC)
int type_of_call killf_(pid, sig)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call killf(pid, sig)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call KILLF(pid, sig)
#endif
      int  *pid, *sig;
{
#ifndef WIN32
      int  kill();
      int  pidu, sigu, istat;

      pidu = *pid;
      sigu = *sig;
      istat = kill(pidu, sigu);
      return istat;
#else
      HANDLE hProcess;
      BOOL TermSucc;

      hProcess= OpenProcess(PROCESS_ALL_ACCESS, TRUE, *pid);
      if (hProcess == NULL)
      TermSucc= TerminateProcess(hProcess, -1);
      return -1;
#endif
}
/*> END <----------------------------------------------------------*/
