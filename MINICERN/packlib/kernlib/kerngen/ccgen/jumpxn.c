/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:11  mclareni
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
#include "vaogs/jumpxn.c"
#else
/*>    ROUTINE JUMPXN
  CERN PROGLIB# Z042    JUMPXN          .VERSION KERNFOR  4.40  940929
  ORIG. 21/04/88 JZ+FCA
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    three steps are needed :

C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)

C- 2) CALL JUMPST (IADR)           to set the tranfer address

C- 3) CALL JUMPX2 (par1,par2)      to transfer
*/
#if defined(CERNLIB_QCCINDAD)
#define IADR *iadr
#endif
#if !defined(CERNLIB_QCCINDAD)
#define IADR iadr
#endif
#if defined(CERNLIB_QX_SC)
#define jumpad type_of_call jumpad_
#define jumpst type_of_call jumpst_
#define jumpx0 type_of_call jumpx0_
#define jumpx1 type_of_call jumpx1_
#define jumpx2 type_of_call jumpx2_
#define jumpx3 type_of_call jumpx3_
#define jumpx4 type_of_call jumpx4_
#define jumpx5 type_of_call jumpx5_
#define jumpx6 type_of_call jumpx6_
#define jumpx7 type_of_call jumpx7_
#define jumpx8 type_of_call jumpx8_
#define jumpx9 type_of_call jumpx9_
#endif
#if defined(CERNLIB_QXCAPT)
#define jumpad type_of_call JUMPAD
#define jumpst type_of_call JUMPST
#define jumpx0 type_of_call JUMPX0
#define jumpx1 type_of_call JUMPX1
#define jumpx2 type_of_call JUMPX2
#define jumpx3 type_of_call JUMPX3
#define jumpx4 type_of_call JUMPX4
#define jumpx5 type_of_call JUMPX5
#define jumpx6 type_of_call JUMPX6
#define jumpx7 type_of_call JUMPX7
#define jumpx8 type_of_call JUMPX8
#define jumpx9 type_of_call JUMPX9
#endif

static void  (type_of_call *tarsub)();

/* ----   jumpad   ---------------------------------------------  */
int jumpad(ifun)
    char *ifun;
{
    return (int) ifun;
}
/* ----   jumpst   ---------------------------------------------  */
void jumpst(iadr)
     void (type_of_call **IADR)();
{
    tarsub = *IADR;
}
/* ----   jumpxn   ---------------------------------------------  */
void jumpx0()
{
    (*tarsub)();
    return;
}

void jumpx1(ixa)
     char *ixa;
{
    (*tarsub)(ixa);
    return;
}

void jumpx2(ixa, ixb)
     char *ixa, *ixb;
{
    (*tarsub)(ixa, ixb);
    return;
}
void jumpx3(ixa, ixb, ixc)
     char *ixa, *ixb, *ixc;
{
    (*tarsub)(ixa, ixb, ixc);
    return;
}
void jumpx4(ixa, ixb, ixc, ixd)
     char *ixa, *ixb, *ixc, *ixd;
{
    (*tarsub)(ixa, ixb, ixc, ixd);
    return;
}
void jumpx5(ixa, ixb, ixc, ixd, ixe)
     char *ixa, *ixb, *ixc, *ixd, *ixe;
{
    (*tarsub)(ixa, ixb, ixc, ixd, ixe);
    return;
}
void jumpx6(ixa, ixb, ixc, ixd, ixe, ixf)
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf;
{
    (*tarsub)(ixa, ixb, ixc, ixd, ixe, ixf);
    return;
}
void jumpx7(ixa, ixb, ixc, ixd, ixe, ixf, ixg)
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg;
{
    (*tarsub)(ixa, ixb, ixc, ixd, ixe, ixf, ixg);
    return;
}
void jumpx8(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh)
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg, *ixh;
{
    (*tarsub)(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh);
    return;
}
void jumpx9(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh, ixi)
     char *ixa, *ixb, *ixc, *ixd, *ixe, *ixf, *ixg, *ixh, *ixi;
{
    (*tarsub)(ixa, ixb, ixc, ixd, ixe, ixf, ixg, ixh, ixi);
    return;
}
/*> END <----------------------------------------------------------*/
#endif
