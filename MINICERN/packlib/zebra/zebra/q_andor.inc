*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_ANDOR_INC
#define CERNLIB_ZEBRA_Q_ANDOR_INC
*
*    IAND / IOR / IEOR / NOT  standard
*
* q_andor.inc
*
#if defined(CERNLIB_APOLLO)
      IAND (IZV,IZW) = AND (IZV, IZW)
      IOR  (IZV,IZW) =  OR (IZV, IZW)
      IEOR (IZV,IZW) = XOR (IZV, IZW)
*    NOT  is standard on APOLLO

#elif defined(CERNLIB_CDC)
      IAND (IZV,IZW) = IZV .AND. IZW
      IOR  (IZV,IZW) = IZV .OR.  IZW
      IEOR (IZV,IZW) = IZV .XOR. IZW
      NOT (IZW)      = .NOT. IZW

#elif defined(CERNLIB_MSDOS)&&defined(CERNLIB_QF_F2C)
      IAND(IZV,IZW)  = AND(IZV,IZW)
       IOR(IZV,IZW)  =  OR(IZV,IZW)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)

#elif defined(CERNLIB_QMIBMSI)
      IAND (IZV,IZW) = LAND (IZV, IZW)
      IOR  (IZV,IZW) =  LOR (IZV, IZW)
      IEOR (IZV,IZW) = LXOR (IZV, IZW)
      NOT (IZW)      = LCOMPL (IZW)

#elif defined(CERNLIB_LNX)
      IAND(IZV,IZW)  = AND(IZV,IZW)
       IOR(IZV,IZW)  =  OR(IZV,IZW)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)

#elif defined(CERNLIB_SUN)
      IAND (IZV,IZW) = AND (IZV, IZW)
      IOR  (IZV,IZW) =  OR (IZV, IZW)
      IEOR (IZV,IZW) = XOR (IZV, IZW)
*    NOT  is standard on SUN
 
#elif 1
*  IAND / IOR / IEOR / NOT  standard

#endif

#endif
