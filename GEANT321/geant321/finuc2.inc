*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:39  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FINUC2_INC
#define CERNLIB_GEANT321_FINUC2_INC
*
*
* finuc2.inc
*
*
*=== finuc2 ===========================================================*
*
*----------------------------------------------------------------------*
*     Last change on 5 april 1990 by   Alfredo Ferrari, INFN - Milan   *
*     include file: finuc2 copy                  created 26/11/86 by pa*
*     note: if changed change finuc copy also                          *
*     changes: none                                                    *
*     included in the following subroutines or functions: not updated  *
*                                                                      *
*     description of the common block(s) and variable(s)               *
*                                                                      *
*     /finuc/ is the storage for secondaries created in event          *
*        irn       = number of secondaries                             *
*        itrn(ip)  = type of the secondary ip                          *
*        cxrn(ip)  = direction cosine of the secondary ip              *
*                    with respect to x-axis                            *
*        cyrn(ip)  = direction cosine of the secondary ip              *
*                    with respect to y-axis                            *
*        czrn(ip)  = direction cosine of the secondary ip              *
*                    with respect to z-axis                            *
*        elr(ip)   = kinetic (or total) energy of secondary ip         *
*        plr(ip)   = momentum of the secondary ip                      *
*        wei(ip)   = weight of the secondary ip                        *
*        tv        = excitation energy                                 *
*        tvcms     = actual excitation energy of the residual nucleus  *
*        tvrecl    = recoil kinetic energy of the residual nucleus     *
*        tvheav    = recoil kinetic energies of heavy (2-H, 3-H, 3-He, *
*                    4-He) fragments after evaporation                 *
*        tvbind    = approximate energy wasted in nuclear binding      *
*                    effects
*----------------------------------------------------------------------*
*
      PARAMETER (MXP=999)
      COMMON / FKFINU / CXRN (MXP), CYRN (MXP), CZRN (MXP), ELR (MXP),
     &                 PLR  (MXP), WEI  (MXP), TV, TVCMS, TVRECL,
     &                 TVHEAV, TVBIND, NP0, IRN, ITRN (MXP)
 

#endif
