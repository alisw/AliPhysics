/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:32  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE JUMPTN
  CERN PROGLIB# Z043    JUMPTN          .VERSION KERNVMI  1.09  940531
  ORIG. 21/04/88 JZ+FCA
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    two steps are needed :

C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)

C- 3) CALL JUMPT2 (IADR,par1,par2)      to transfer
*/
#define jumpt0 jumpt0_
#define jumpt1 jumpt1_
#define jumpt2 jumpt2_
#define jumpt3 jumpt3_
#define jumpt4 jumpt4_

extern int jumpad_();

static void (*jumpto)();

void jumpt0(iadr)
     int *iadr;
{
    long func;
    func = *iadr + (long)jumpad_;
    jumpto = (void(*)()) func;
    jumpto();
    return;
}

void jumpt1(iadr,ipara)
     int *iadr;
     char *ipara;
{
    long func;
    func = *iadr + (long)jumpad_;
    jumpto = (void(*)()) func;
    jumpto (ipara);
    return;
}

void jumpt2(iadr, ipara, iparb)
     int *iadr;
     char *ipara, *iparb;
{
    long func;
    func = *iadr + (long)jumpad_;
    jumpto = (void(*)()) func;
    jumpto (ipara, iparb);
    return;
}
void jumpt3(iadr, ipara, iparb, iparc)
     int *iadr;
     char *ipara, *iparb, *iparc;
{
    long func;
    func = *iadr + (long)jumpad_;
    jumpto = (void(*)()) func;
    jumpto (ipara, iparb, iparc);
    return;
}
void jumpt4(iadr, ipara, iparb, iparc, ipard)
     int *iadr;
     char *ipara, *iparb, *iparc, *ipard;
{
    long func;
    func = *iadr + (long)jumpad_;
    jumpto = (void(*)()) func;
    jumpto (ipara, iparb, iparc, ipard);
    return;
}
/*> END <----------------------------------------------------------*/
