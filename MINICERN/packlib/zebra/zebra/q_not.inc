*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:12  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_NOT_INC
#define CERNLIB_ZEBRA_Q_NOT_INC
*
* q_not.inc
*
#if defined(CERNLIB_APOLLO)
*    NOT  is standard on APOLLO

#elif defined(CERNLIB_CDC)
      NOT (IZW)      = .NOT. IZW

#elif defined(CERNLIB_MSDOS)&&defined(CERNLIB_QF_F2C)

#elif defined(CERNLIB_QMIBMSI)
      NOT (IZW)      = LCOMPL (IZW)

#elif defined(CERNLIB_LNX)

#elif defined(CERNLIB_SUN)
*    NOT  is standard on SUN
 
#elif 1

#endif

#endif
