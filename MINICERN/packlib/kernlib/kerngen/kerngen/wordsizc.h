/*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*
* wordsizc.h
*/
#if defined(CERNLIB_QMIRTD)
#define NBITPW 64      /* Number of bits  per word */
#define NBYTPW 8       /* Number of bytes per word */
#else
#define NBYTPW 4       /* Number of bytes per word */
#endif
