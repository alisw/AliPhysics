*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:27:43  mclareni
* Import version 8.04
*
* Revision 1.3  2000/04/15 14:10:42  plothow
* Version 8.03: Nuclear PDFs introduced
*
* Revision 1.2  2000/04/09 10:06:18  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
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
C.      PARAMETER (LPTYPE = 1, LGROUP = 3, LNSET = 41)
C.      PARAMETER (IFL = -5, LORD = 2, TM = 180.D0)
C. NEW default (2000) PDF is MRS set (c-g) (NPTYPE=1, NGROUP=3, NSET = 89)
      PARAMETER (LPTYPE = 1, LGROUP = 3, LNSET = 89)
      PARAMETER (IFL = -5, LORD = 2, TM = 180.D0)
C. NEW default (2000) NPDF is EKS98 set 1 (NATYPE=4, NAGROUP=1, NASET = 1)
      PARAMETER (LATYPE = 4, LAGROUP = 1, LNASET = 1)
