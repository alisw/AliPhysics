/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:37  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:38  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMSUN)
#include "sungs/qnexte.c"
#else
/*>    ROUTINE QNEXTE
  CERN PROGLIB# Z041    QNEXTE          .VERSION KERNFOR  4.29  910718
* ORIG. 23/05/91, JZ
*/
#include <setjmp.h>

#if defined(CERNLIB_QX_SC)
void type_of_call qnext_();
#endif

#if defined(CERNLIB_QXNO_SC)
void type_of_call qnext();
#endif

#if defined(CERNLIB_QXCAPT)
void type_of_call QNEXT();
#endif

#if defined(CERNLIB_QX_SC)
void type_of_call qnexte_()
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call qnexte()
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call QNEXTE()
#endif
#if defined(CERNLIB_QSIGJMP)
{     static sigjmp_buf  myenv;
      static int ireent = 0;

      if (ireent)  siglongjmp (myenv, 7);

      ireent = 77;
      sigsetjmp (myenv,7);
#endif
#if !defined(CERNLIB_QSIGJMP)
{     static jmp_buf  myenv;
      static int ireent = 0;

      if (ireent)  longjmp (myenv, 7);

      ireent = 77;
      setjmp (myenv);
#endif
#if defined(CERNLIB_QX_SC)
      qnext_();
#endif
#if defined(CERNLIB_QXNO_SC)
      qnext();
#endif
#if defined(CERNLIB_QXCAPT)
      QNEXT();
#endif
}
/*> END <----------------------------------------------------------*/
#endif
