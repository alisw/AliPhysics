*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:50:55  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:11  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:49  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:18  mclareni
* Kernlib
*
*
*
* wordsize.inc
*
*    # of bits/word
*    # of characters/word
#if defined(CERNLIB_QMIRTD)
      PARAMETER   (NBITPW=64, NCHAPW=8, NBITPC=8)
#else
      PARAMETER   (NBITPW=32)
      PARAMETER   (NCHAPW=4)
#endif
