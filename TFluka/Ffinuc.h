#ifndef FFINUC
#define FFINUC_H 1
                                                                                
#include "Rtypes.h"
#include "cfortran.h"
                                                                                
#include "Fdimpar.h"
extern "C" {
//*$ create finuc.add
//*copy finuc
//*
//*=== finuc ============================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     include file: finuc (new version of old finuc of fluka86)        *
//*                                                                      *
//*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      *
//*     !!!!     s e e   a l s o   i n c l u d e   f i l e     !!!!      *
//*     !!!!                 f i n u c 2                       !!!!      *
//*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      *
//*                                                                      *
//*     created on  20 january 1996  by    alfredo ferrari & paola sala  *
//*                                                   infn - milan       *
//*                                                                      *
//*     last change on 26-jul-97     by    alfredo ferrari               *
//*                                                                      *
//*     included in the following subroutines or functions: not updated  *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*     /finuc/ is the storage for secondaries created in event          *
//*        np        = number of secondaries                             *
//*       kpart (ip) = type of the secondary ip                          *
//*         cxr (ip) = direction cosine of the secondary ip              *
//*                    with respect to x-axis                            *
//*         cyr (ip) = direction cosine of the secondary ip              *
//*                    with respect to y-axis                            *
//*         czr (ip) = direction cosine of the secondary ip              *
//*                    with respect to z-axis                            *
//*      cxrpol (ip) = direction cosine of the secondary ip polarization *
//*                    with respect to x-axis                            *
//*      cyrpol (ip) = direction cosine of the secondary ip polarization *
//*                    with respect to y-axis                            *
//*      czrpol (ip) = direction cosine of the secondary ip polarization *
//*                    with respect to z-axis                            *
//*         tki (ip) = kinetic energy of secondary ip                    *
//*         plr (ip) = momentum of the secondary ip                      *
//*         wei (ip) = weight of the secondary ip                        *
//*      agesec (ip) = "age" of the secondary ip with respect to the     *
//*                    interaction time                                  *
//*        tv        = excitation energy                                 *
//*        tvcms     = actual excitation energy of the residual nucleus  *
//*        tvrecl    = recoil kinetic energy of the residual nucleus     *
//*        tvheav    = recoil kinetic energies of heavy (2-h, 3-h, 3-he, *
//*                    4-he) fragments after evaporation                 *
//*        tvbind    = approximate energy wasted in nuclear binding      *
//*                    effects (not yet operational)                     *
//*                                                                      *
//*----------------------------------------------------------------------*
//*
const Int_t mxp = mxpscs;
//*

typedef struct {
   Double_t cxr[mxp];
   Double_t cyr[mxp];
   Double_t czr[mxp];
   Double_t cxrpol[mxp];
   Double_t cyrpol[mxp];
   Double_t czrpol[mxp];
   Double_t tki[mxp];
   Double_t plr[mxp];
   Double_t wei[mxp];
   Double_t agesec[mxp];
   Double_t tv;
   Double_t tvcms;
   Double_t tvrecl;
   Double_t tvheav;
   Double_t tvbind;
   Int_t    np0;
   Int_t    np;
   Int_t    kpart[mxp];
} finucCommon;
#define FINUC COMMON_BLOCK(FINUC,finuc)
COMMON_BLOCK_DEF(finucCommon,FINUC);
}
#endif
