/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:10  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1997/02/04 17:34:23  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:31  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:24  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"

#if defined(CERNLIB_QMVAOS)
#include "vaogs/jumptn.c"
#else
/*>    ROUTINE JUMPTN
  CERN PROGLIB# Z043    JUMPTN          .VERSION KERNFOR  4.40  940929
  ORIG. 21/04/88 JZ+FCA
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    two steps are needed :

C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)

C- 3) CALL JUMPT2 (IADR,par1,par2)      to transfer
*/
#if defined(CERNLIB_QCCINDAD)
#define IADR *iadr
#endif
#if !defined(CERNLIB_QCCINDAD)
#define IADR iadr
#endif
#if defined(CERNLIB_QX_SC)
#define jumpt0 type_of_call jumpt0_
#define jumpt1 type_of_call jumpt1_
#define jumpt2 type_of_call jumpt2_
#define jumpt3 type_of_call jumpt3_
#define jumpt4 type_of_call jumpt4_
#define jumpt5 type_of_call jumpt5_
#define jumpt6 type_of_call jumpt6_
#define jumpt7 type_of_call jumpt7_
#define jumpt8 type_of_call jumpt8_
#define jumpt9 type_of_call jumpt9_
#endif
#if defined(CERNLIB_QXCAPT)
#define jumpt0 type_of_call JUMPT0
#define jumpt1 type_of_call JUMPT1
#define jumpt2 type_of_call JUMPT2
#define jumpt3 type_of_call JUMPT3
#define jumpt4 type_of_call JUMPT4
#define jumpt5 type_of_call JUMPT5
#define jumpt6 type_of_call JUMPT6
#define jumpt7 type_of_call JUMPT7
#define jumpt8 type_of_call JUMPT8
#define jumpt9 type_of_call JUMPT9
#endif
void jumpt0(iadr)
     void (type_of_call **IADR)();
{
    (**IADR)();
    return;
}

void jumpt1(iadr,ixa)
     void (type_of_call **IADR)();
     char *ixa;
{
    (**IADR)(ixa);
    return;
}

void jumpt2(iadr, ixa, ixb)
     void (type_of_call **IADR)();
     char *ixa, *ixb;
{
    (**IADR)(ixa, ixb);
    return;
}
void jumpt3(iadr, ixa, ixb, ixc)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc;
{
    (**IADR)(ixa, ixb, ixc);
    return;
}
void jumpt4(iadr, ixa, ixb, ixc, ixd)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd;
{
    (**IADR)(ixa, ixb, ixc, ixd);
    return;
}
void jumpt5(iadr, ixa, ixb, ixc, ixd, ixe)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd, *ixe;
{
    (**IADR)(ixa, ixb, ixc, ixd, ixe);
    return;
}
void jumpt6(iadr, ixa, ixb, ixc, ixd, ixe, ixf)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf;
{
    (**IADR)(ixa, ixb, ixc, ixd, ixe, ixf);
    return;
}
void jumpt7(iadr, ixa, ixb, ixc, ixd, ixe, ixf, ixg)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg;
{
    (**IADR)(ixa, ixb, ixc, ixd, ixe, ixf, ixg);
    return;
}
void jumpt8(iadr, ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg, *ixh;
{
    (**IADR)(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh);
    return;
}
void jumpt9(iadr, ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh, ixi)
     void (type_of_call **IADR)();
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg, *ixh, *ixi;
{
    (**IADR)(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh, ixi);
    return;
}
/*> END <----------------------------------------------------------*/
#endif
