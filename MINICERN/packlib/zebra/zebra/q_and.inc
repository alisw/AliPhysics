*
* $Id$
*
* $Log$
* Revision 1.1  1996/04/18 16:15:05  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_AND_INC
#define CERNLIB_ZEBRA_Q_AND_INC
*
* q_and.inc
*
#if defined(CERNLIB_APOLLO)
      IAND (IZV,IZW) = AND (IZV, IZW)

#elif defined(CERNLIB_CDC)
      IAND (IZV,IZW) = IZV .AND. IZW

#elif defined(CERNLIB_MSDOS)&&defined(CERNLIB_QF_F2C)
      IAND(IZV,IZW)  = AND(IZV,IZW)

#elif defined(CERNLIB_QMIBMSI)
      IAND (IZV,IZW) = LAND (IZV, IZW)

#elif defined(CERNLIB_LNX)
      IAND(IZV,IZW)  = AND(IZV,IZW)

#elif defined(CERNLIB_SUN)
      IAND (IZV,IZW) = AND (IZV, IZW)
 
#elif 1

#endif

#endif
