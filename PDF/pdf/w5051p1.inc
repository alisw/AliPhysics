*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:20  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:48  plothow
* Version 7.01
*
*
*
* w5051p1.inc
*
C. OLD default PDF was MRS set S0 (MODE = 45 : NPTYPE=1, NGROUP=3, NSET = 26)
      PARAMETER (LPTYPO = 1, LGROPO = 3, LNSETO = 26)
C. NEW default PDF is GRV set HO (NPTYPE=1, NGROUP=5, NSET = 3)
C.    PARAMETER (LPTYPE = 1, LGROUP = 5, LNSET = 3)
C.    PARAMETER (IFL = -5, LORD = 2, TM = 100.D0)
C. NEW default (1995) PDF is MRS set (G) (NPTYPE=1, NGROUP=3, NSET = 41)
      PARAMETER (LPTYPE = 1, LGROUP = 3, LNSET = 41)
      PARAMETER (IFL = -5, LORD = 2, TM = 180.D0)
