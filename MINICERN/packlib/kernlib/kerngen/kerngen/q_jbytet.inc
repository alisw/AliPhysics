*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:40  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:29  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* q_jbytet.inc
*
#if defined(CERNLIB_QMSUN)&&defined(CERNLIB_BUGLRSHFT)
#elif defined(CERNLIB_QMSUN)&&(!defined(CERNLIB_BUGLRSHFT))
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                       lrshft (LSHIFT(IZW,33-IZP-NZB), 32-NZB))
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                       lrshft (LSHIFT(IZW,33-IZP-NZB), 32-NZB))
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW, NOT(lrshft (
     +                          LSHIFT(NOT(MZ),32-NZB), 33-IZP-NZB)))
      MBYTOR (MZ,IZW,IZP,NZB) = OR (IZW, lrshft (
     +                               LSHIFT(MZ,32-NZB), 33-IZP-NZB))
#if defined(CERNLIB_QMLNX)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                        ISHFTR (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                        ISHFTR (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW,
     +           NOT(ISHFTR (LSHIFT(NOT(MZ),32-NZB),33-IZP-NZB)) )
      MBYTOR (MZ,IZW,IZP,NZB) = OR (IZW,
     +       ISHFTR (LSHIFT(MZ,32-NZB),33-IZP-NZB))
#if defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
      JBYTET (MZ,IZW,IZP,NZB) = AND (MZ,
     +                        ishftr (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
      JBYTOR (MZ,IZW,IZP,NZB) = OR  (MZ,
     +                        ishftr (LSHIFT(IZW,33-IZP-NZB),32-NZB) )
      MBYTET (MZ,IZW,IZP,NZB) = AND (IZW,
     +           NOT(ishftr (LSHIFT(NOT(MZ),32-NZB),33-IZP-NZB)) )
      MBYTOR (MZ,IZW,IZP,NZB) = OR (IZW,
     +                    ishftr (LSHIFT(MZ,32-NZB),33-IZP-NZB) )
#elif defined(CERNLIB_QMOS9)||defined(CERNLIB_QMHPX)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ, IBITS(IZW,IZP-1,NZB))
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ, IBITS(IZW,IZP-1,NZB))
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW, NOT(ishft (
     +                          ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )
      MBYTOR (MZ,IZW,IZP,NZB) = IOR (IZW, ishft (
     +                               ISHFT(MZ,32-NZB),-33+IZP+NZB) )
#elif defined(CERNLIB_QMDOS)&&(!defined(CERNLIB_QF_F2C)) || defined(CERNLIB_QFMSOFT)
      JBYTET (MZ,IZW,IZP,NZB) = IAND (MZ,
     +                        ISHFT (ISHFT(IZW,33-IZP-NZB),-32+NZB) )
      JBYTOR (MZ,IZW,IZP,NZB) = IOR  (MZ,
     +                        ISHFT (ISHFT(IZW,33-IZP-NZB),-32+NZB) )
      MBYTET (MZ,IZW,IZP,NZB) = IAND (IZW,
     +           NOT(ISHFT (ISHFT(NOT(MZ),32-NZB),-33+IZP+NZB)) )
      MBYTOR (MZ,IZW,IZP,NZB) = IOR (IZW,
     +                    ISHFT (ISHFT(MZ,32-NZB),-33+IZP+NZB) )
#endif
