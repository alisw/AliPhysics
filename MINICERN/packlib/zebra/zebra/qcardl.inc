*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QCARDL_INC
#define CERNLIB_ZEBRA_QCARDL_INC
*
*
* qcardl.inc
*
#if (defined(CERNLIB_A4)||defined(CERNLIB_A5)||defined(CERNLIB_A6)||defined(CERNLIB_A8)||defined(CERNLIB_A10))&&(defined(CERNLIB_A10)||defined(CERNLIB_A8)||defined(CERNLIB_A6)||defined(CERNLIB_A5))
#endif
#if (defined(CERNLIB_A5M)||defined(CERNLIB_A6M)||defined(CERNLIB_A8M))&&(defined(CERNLIB_A6M)||defined(CERNLIB_A8M))&&(defined(CERNLIB_A8M))
#endif
#if (defined(CERNLIB_B32)||defined(CERNLIB_B36)||defined(CERNLIB_B48)||defined(CERNLIB_B60)||defined(CERNLIB_B64))&&(defined(CERNLIB_B64)||defined(CERNLIB_B60)||defined(CERNLIB_B48)||defined(CERNLIB_B36))
#endif
#if (defined(CERNLIB_B36M)||defined(CERNLIB_B48M)||defined(CERNLIB_B60M))&&(defined(CERNLIB_B48M)||defined(CERNLIB_B60M))&&(defined(CERNLIB_B60M))
#endif

#endif
