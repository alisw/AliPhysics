#ifndef FPART_H
#define FPART_H 1

#include "Rtypes.h"
#include "cfortran.h"
#include "Fdimpar.h" //For some constants
extern "C" {
//*$ create part.add
//*copy part
//*
//*=== part =============================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     include file: part copy        revised on 20-08-96 by a. ferrari *
//*                                                                      *
//*     last change   on   14-oct-00     by       alfredo ferrari        *
//*                                                                      *
//*     included in the following subroutines or functions:              *
//*                                                                      *
//*     w a r n i n g !!!! check also part2 and part3 for any change!!!  *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*             am = particle mass (gev/c^2)                             *
//*             ga = particle width (gev)                                *
//*            tau = particle mean life (s)                              *
//*         amdisc = "effective" particle mass for energy balance (gev)  *
//*         zmnabs = lower width (adimensional unit) to be used during   *
//*                  particle decay to assure that at least one decay    *
//*                  channel is physically open                          *
//*         atnmna = atan (zmnabs)                                       *
//*            ich = particle electric charge                            *
//*           ibar = particle baryon number                              *
//*         isosym = index of the isospin reversed (t_z --> -t_z)        *
//*                  particle (if any, if 0 no such particle is available*
//*                  in the part listing)                                *
//*         ichcon = index of the charge conjugated (antiparticle)       *
//*                  particle (if any, if 0 no such particle is available*
//*                  in the part listing)                                *
//*             k1 = index of first decay channel                        *
//*             k2 = index of last  decay channel                        *
//*         kptoip = conversion from part to paprop numbering            *
//*         iptokp = conversion from paprop to part numbering            *
//*         kptoia = conversion from part to abltis numbering            *
//*         iatokp = conversion from abltis to part numbering            *
//*         idcflg = decay flag                                          *
//*         iptype = particle type                                       *
//*                  -1: heavy fragments                                 *
//*                   0: unknown particle or lepton                      *
//*                   1: nucleon                                         *
//*                   2: antinucleon                                     *
//*                   3: pion                                            *
//*                   4: k+/k0                                           *
//*                  -4: kshrt/klong                                     *
//*                   5: k-/k0bar                                        *
//*                   6: lamda/sigma   (strangeness -1 hyperon)          *
//*                   7: xsi           (strangeness -2 hyperon)          *
//*                   8: omega         (strangeness -3 hyperon)          *
//*                   9: alamda/asigma (strangeness +1 antihyperon)      *
//*                  10: axsi          (strangeness +2 antihyperon)      *
//*                  11: aomega        (strangeness +3 antihyperon)      *
//*                  12: d+/d0                                           *
//*                  13: d-/d0bar                                        *
//*                  14: d_s+/d_s-                                       *
//*                  15: lambda_c+                                       *
//*                  16: xsi_c+/xsi_c0                                   *
//*                  17: xsi'_c+/xsi'_c0                                 *
//*                  18: omega_c                                         *
//*                  19: alambda_c+                                      *
//*                  20: axsi_c-/axsi_c0                                 *
//*                  21: axsi'_c-/axsi'_c0                               *
//*                  22: aomega_c                                        *
//*          aname = particle literal name                               *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Double_t am[idmaxp+7];
   Double_t ga[idmaxp+7];
   Double_t tau[idmaxp+7];
   Double_t amdisc[idmaxp+7];
   Double_t zmnabs[idmaxp+7];
   Double_t atnmna[idmaxp+7];
   Int_t    ich[idmaxp+7];
   Int_t    ibar[idmaxp+7];
   Int_t    isosym[idmaxp+7];
   Int_t    ichcon[idmaxp+7];
   Int_t    k1[idmaxp+7];
   Int_t    k2[idmaxp+7];
   Int_t    kptoip[idmaxp+7];
   Int_t    iptokp[nallwp+7];
   Int_t    kptoia[idmaxp+7];
   Int_t    iatokp[mxpabl+7];
   Int_t    idcflg[nallwp+7];
   Int_t    iptype[nallwp+7];
} partCommon;
#define PART COMMON_BLOCK(PART,part)
COMMON_BLOCK_DEF(partCommon,PART);

typedef struct {
   Char_t   aname[idmaxp+7][8];
} chpartCommon;
#define CHPART COMMON_BLOCK(CHPART,chpart)
COMMON_BLOCK_DEF(chpartCommon,CHPART);
}

//Get functions
inline Double_t GetFlukaAM(unsigned int i) {return PART.am[i+6];}
inline Double_t GetFlukaGA(unsigned int i) {return PART.ga[i+6];}
inline Double_t GetFlukaTAU(unsigned int i) {return PART.tau[i+6];}
inline Double_t GetFlukaAMDISC(unsigned int i) {return PART.amdisc[i+6];}
inline Double_t GetFlukaZMNABS(unsigned int i) {return PART.zmnabs[i+6];}
inline Double_t GetFlukaATNMNA(unsigned int i) {return PART.atnmna[i+6];}
inline Int_t    GetFlukaICH(unsigned int i) {return PART.ich[i+6];}
inline Int_t    GetFlukaIBAR(unsigned int i) {return PART.ibar[i+6];}
inline Int_t    GetFlukaISOSYM(unsigned int i) {return PART.isosym[i+6];}
inline Int_t    GetFlukaICHCON(unsigned int i) {return PART.ichcon[i+6];}
inline Int_t    GetFlukaK1(unsigned int i) {return PART.k1[i+6];}
inline Int_t    GetFlukaK2(unsigned int i) {return PART.k2[i+6];}
inline Int_t    GetFlukaKPTOIP(unsigned int i) {return PART.kptoip[i+6];}
inline Int_t    GetFlukaIPTOKP(unsigned int i) {return PART.iptokp[i+6];}
inline Int_t    GetFlukaKPTOIA(unsigned int i) {return PART.kptoia[i+6];}
inline Int_t    GetFlukaIATOKP(unsigned int i) {return PART.iatokp[i+6];}
inline Int_t    GetFlukaIDCFLG(unsigned int i) {return PART.idcflg[i+6];}
inline Int_t    GetFlukaIPTYPE(unsigned int i) {return PART.iptype[i+6];}

#endif
