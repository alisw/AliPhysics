/*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:35  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:25  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* cf#clos.inc
*/
#if defined(CERNLIB_QMAPO)||defined(CERNLIB_QMOS9)
#elif defined(CERNLIB_QMVAX)
#include <file.h>            /*  VAX/VMS                    */
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
  #ifdef WIN32
#   ifdef __STDC__
#    undef __STDC__
#   endif
    #include <io.h>
  #endif
#else
#include <unistd.h>          /*  default Posix              */
#endif
