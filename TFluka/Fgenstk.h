#ifndef FGENSTK
#define FGENSTK_H 1
                                                                                
#include "Rtypes.h"
#include "cfortran.h"
                                                                                
#include "Fdimpar.h"
//*$ CREATE GENSTK.ADD
//*COPY GENSTK
//*
//*=== Genstk ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     hadron GENerator STacK for FLUKA: (new version of old Finuc of   *
//*     FLUKA86 by J.Ranft)                                              *
//*                                                                      *
//*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      *
//*     !!!!     S E E   A L S O   I N C L U D E   F I L E     !!!!      *
//*     !!!!                 G E N S T K 2                     !!!!      *
//*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      *
//*                                                                      *
//*     Created on  20 january 1996  by    Alfredo Ferrari & Paola Sala  *
//*                                                   Infn - Milan       *
//*                                                                      *
//*     Last change on 15-may-05     by    Alfredo Ferrari               *
//*                                                                      *
//*                                                                      *
//*     /Genstk/ is the storage for secondaries created in hadronic      *
//*              events                                                  *
//*        Np        = total number of secondaries                       *
//*       Kpart (ip) = (Paprop) id of the ip_th secondary                *
//*         Cxr (ip) = x-axis direction cosine of the ip_th secondary    *
//*         Cyr (ip) = y-axis direction cosine of the ip_th secondary    *
//*         Czr (ip) = z-axis direction cosine of the ip_th secondary    *
//*      Cxrpol (ip) = x-axis direction cosine of the ip_th secondary    *
//*                    polarization vector (rest frame when applicable)  *
//*      Cyrpol (ip) = y-axis direction cosine of the ip_th secondary    *
//*                    polarization vector (rest frame when applicable)  *
//*      Czrpol (ip) = z-axis direction cosine of the ip_th secondary    *
//*                    polarization vector (rest frame when applicable)  *
//*         Tki (ip) = laboratory kinetic energy of ip_th secondary (GeV)*
//*         Plr (ip) = laboratory momentum of the ip_th secondary (GeV/c)*
//*         Wei (ip) = statistical weight of the ip_th secondary         *
//*      Agesec (ip) = "age" of the ip_th secondary with respect to the  *
//*                    interaction time                                  *
//*        Tv        = excitation energy (GeV)                           *
//*        Tvcms     = actual excitation energy of the residual nucleus  *
//*        Tvrecl    = recoil kinetic energy of the residual nucleus     *
//*        Tvheav    = recoil kinetic energies of heavy (2-H, 3-H, 3-He, *
//*                    4-He) fragments after evaporation                 *
//*        Tvbind    = approximate energy wasted in nuclear binding      *
//*                    effects (not yet operational)                     *
//*      Infext (ip) = possible extra infos for the ip_th secondary      * 2006.3
//*                                                                      *
//*----------------------------------------------------------------------*

extern "C" {
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
    Int_t    kpart[mxp];
    Int_t    infext[mxpscs];
    Int_t    np0;
    Int_t    np;
} genstkCommon;
#define GENSTK COMMON_BLOCK(GENSTK,genstk)
COMMON_BLOCK_DEF(genstkCommon,GENSTK);
}
#endif
