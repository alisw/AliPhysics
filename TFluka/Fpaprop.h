#ifndef FPAPROP_H
#define FPAPROP_H 1

#include "Rtypes.h"
#include "cfortran.h"

#include "Fdimpar.h"

extern "C" {
//*$ create paprop.add
//*copy paprop
//*
//*=== paprop ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     include file: paprop copy                                        *
//*                                                                      *
//*     !!!!    n e w   v e r s i o n   !!!!                             *
//*                                                                      *
//*     created on    07 may 1991    by    alfredo ferrari & paola sala  *
//*                                                   infn - milan       *
//*                                                                      *
//*     last change on 03-jul-97     by    alfredo ferrari               *
//*                                                                      *
//*     included in the following subroutines or functions: not updated  *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*     /paprop/ contains particle properties                            *
//*        btype  = literal name of the particle                         *
//*        am     = particle mass in gev                                 *
//*        ichrge = electric charge of the particle                      *
//*        ibarch = baryonic charge of the particle                      *
//*        iscore = explanations for the scored distribution             *
//*        genpar = names of the generalized particles                   *
//*        ijdisc = list of the particle types to be discarded           *
//*        thalf  = half life of the particle in sec                     *
//*        biasdc = decay biasing factors                                *
//*        biasin = inelastic interaction biasing factors                *
//*        lhadro = flag for hadrons                                     *
//*        jspinp = particle spin (in units of 1/2)                      *
//*        iparty = particle parity (when defined)                       *
//*        iparid = flag used to identify particle types                 *
//*        lbsdcy = logical flag for biased decay: if .true. the biasing *
//*                 factor is used as an upper limit to the decay length *
//*        lprbsd = logical flag for biased decay: if .true. the biasing *
//*                 factor is applied only to primaries                  *
//*        lprbsi = logical flag for inelastic interaction biasing: if   *
//*                 .true. the biasing factor is applied only to prima-  *
//*                 ries                                                 *
//*        lsclwf = logical flag for low energy neutron fission scoring  *
//*        lscnbl = logical flag for neutron balance scoring             *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

const Int_t mxgnpr =  33;
typedef struct {
   Double_t am[nallwp+7];         //(-6:NALLWP)
   Double_t amdisc[nallwp+7];     //(-6:NALLWP)
   Double_t thalf[nallwp+7];      //(-6:NALLWP)
   Double_t biasdc[nallwp+7];     //(-6:NALLWP)
   Double_t biasin[nallwp+7];     //(-6:NALLWP)
   Int_t    ichrge[nallwp+7];     //(-6:NALLWP)
   Int_t    ibarch[nallwp+7];     //(-6:NALLWP)
   Int_t    ijdisc[nallwp+7];     //(-6:NALLWP)
   Int_t    jspinp[nallwp+7];     //(-6:NALLWP)
   Int_t    iparty[nallwp+7];     //(-6:NALLWP)
   Int_t    iparid[nallwp+7];     //(-6:NALLWP)
   Int_t    lhadro[nallwp+7];     //(-6:NALLWP)
   Int_t    lbsdcy[nallwp+7];     //(-6:NALLWP)
   Int_t    iscore[10];
   Int_t    lprbsd;
   Int_t    lprbsi;
   Int_t    lsclwf;
   Int_t    lscnbl;
} papropCommon;
#define PAPROP COMMON_BLOCK(PAPROP,paprop)
COMMON_BLOCK_DEF(papropCommon,PAPROP);

typedef struct {
   Char_t   btype[nallwp+7][8];     //(-6:NALLWP)
   Char_t   genpar[mxgnpr][8];          //(30)
} chpprpCommon;
#define CHPPRP COMMON_BLOCK(CHPPRP,chpprp)
COMMON_BLOCK_DEF(chpprpCommon,CHPPRP);
}

//Get functions
//inline Double_t GetFlukaAM(unsigned int i) {return PAPROP.am[i+6];}
//inline Double_t GetFlukaAMDISC(unsigned int i) {return PAPROP.amdisc[i+6];}
//inline Double_t GetFlukaTHALH(unsigned int i) {return PAPROP.thalf[i+6];}
//inline Double_t GetFlukaBIASDC(unsigned int i) {return PAPROP.biasdc[i+6];}
//inline Double_t GetFlukaBIASIN(unsigned int i) {return PAPROP.biasin[i+6];}
//inline Int_t    GetFlukaICHRGE(unsigned int i) {return PAPROP.ichrge[i+6];}
//inline Int_t    GetFlukaIBARCH(unsigned int i) {return PAPROP.ibarch[i+6];}
//inline Int_t    GetFlukaIJDISC(unsigned int i) {return PAPROP.ijdisc[i+6];}
//inline Int_t    GetFlukaJSPINP(unsigned int i) {return PAPROP.jspinp[i+6];}
//inline Int_t    GetFlukaIPARTY(unsigned int i) {return PAPROP.iparty[i+6];}
//inline Int_t    GetFlukaIPARID(unsigned int i) {return PAPROP.iparid[i+6];}
//inline Int_t    GetFlukaLHADRO(unsigned int i) {return PAPROP.lhadro[i+6];}
//inline Int_t    GetFlukaLBSDCY(unsigned int i) {return PAPROP.lbsdcy[i+6];}
#endif
