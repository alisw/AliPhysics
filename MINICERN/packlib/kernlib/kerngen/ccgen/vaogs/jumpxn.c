/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:32  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE JUMPXN
  CERN PROGLIB# Z042    JUMPXN          .VERSION KERNVMI  1.08  930527
  ORIG. 21/04/88 JZ+FCA, adapted 11/05/93 AP+JZ
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    three steps are needed :

C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)

C- 2) CALL JUMPST (IADR)           to set the tranfer address

C- 3) CALL JUMPX2 (par1,par2)      to transfer
*/

static void (*tarsub)();

/* ----   jumpad   ---------------------------------------------  */
int jumpad_(ifun)
    char *ifun;
{
    long temp;

    temp = (long)ifun - (long)jumpad_;
    return (int) temp;
}

/* ----   jumpst   ---------------------------------------------  */
void jumpst_(iadr)
    int  *iadr;
{
    long true;

    true = (long)jumpad_;
    true = true + *iadr;
    tarsub = (void (*)())true;
}

/* ----   jumpxn   ---------------------------------------------  */
jumpx0_()
{
    (*tarsub)();
    return;
}

jumpx1_(ipara)
    char *ipara;
{
    (*tarsub)(ipara);
    return;
}

jumpx2_(ipara, iparb)
    char *ipara, *iparb;
{
    (*tarsub)(ipara, iparb);
    return;
}
jumpx3_(ipara, iparb, iparc)
    char *ipara, *iparb, *iparc;
{
    (*tarsub)(ipara, iparb, iparc);
    return;
}
jumpx4_(ipara, iparb, iparc, ipard)
    char *ipara, *iparb, *iparc, *ipard;
{
    (*tarsub)(ipara, iparb, iparc, ipard);
    return;
}
/*> END <----------------------------------------------------------*/
