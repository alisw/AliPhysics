*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:22:25  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:34:16  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_LOCF_INC
#define CERNLIB_ZEBRA_Q_LOCF_INC
*
*
* q_locf.inc
*
#if defined(CERNLIB_CRAY)
      LOCF (IZLOC) = LOC (IZLOC)
#elif defined(CERNLIB_QFMSOFT_TEST)
*    This is very dangerous since the DWORD address may be odd too  !!!
      LOCF (IZLOC) = LOC (IZLOC)/4
#endif

#endif
