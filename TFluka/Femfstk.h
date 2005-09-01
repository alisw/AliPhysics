#ifndef FEMFSTK
#define FEMFSTK_H 1

#include "Rtypes.h"
#include "cfortran.h"

#include "Fdimpar.h"
extern "C" {
//*$ create emfstk.add
//*copy emfstk
//*
//*=== emfstk ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     common emfstk (emf stack) for emf                                *
//*                                                                      *
//*     last change on  08-oct-97    by    alfredo ferrari               *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Double_t etemf[mestck]; // total energy in MeV
   Double_t pmemf[mestck];
   Double_t xemf[mestck]; // particle x-coordinate
   Double_t yemf[mestck]; // particle y-coordinate
   Double_t zemf[mestck]; // particle z-coordinate
   Double_t uemf[mestck]; // x direction cosine
   Double_t vemf[mestck]; // y direction cosine
   Double_t wemf[mestck]; // z direction cosine
   Double_t dnear[mestck]; // equivalent to GEANT "safety"
   Double_t upol[mestck]; // polarisation in x direction
   Double_t vpol[mestck]; // polarisation in y direction
   Double_t wpol[mestck]; // polarisation in z direction
   Double_t usnrml[mestck];
   Double_t vsnrml[mestck];
   Double_t wsnrml[mestck];
   Double_t wtemf[mestck]; // weight
   Double_t agemf[mestck]; // age
   Double_t cmpemf[mestck];
   Double_t espark[mestck][mkbmx1];
   Double_t rdlyem[mestck];
   Int_t    iespak[mestck][mkbmx2];
   Int_t    ichemf[mestck]; // charge
   Int_t    iremf[mestck];  // region
   Int_t    irlatt[mestck]; // lattice cell
   Int_t    nhpemf[mestck];
   Int_t    lloemf[mestck]; // generation number
   Int_t    louemf[mestck];
   Int_t    lrdemf[mestck];
   Int_t    npemf;  // number of particles in stack
   Int_t    npstrt; // EMF stack index before the interaction (since
                    // the projectile disappears it is also the starting
                    // index of secondaries)
//*d === obsolete variable names === *
//*d     parameter ( idmemf = mestck )
//*d     dimension e (idmemf), wt (idmemf), iq (idmemf), ir (idmemf)
//*d     equivalence ( e   (1), etemf  (1) )
//*d     equivalence ( wt  (1), wtemf  (1) )
//*d     equivalence ( ir  (1), iremf  (1) )
//*d     equivalence ( iq  (1), ichemf (1) )
//*d     equivalence ( np, npemf )
} emfstkCommon;
#define EMFSTK COMMON_BLOCK(EMFSTK,emfstk)
COMMON_BLOCK_DEF(emfstkCommon,EMFSTK);
}
#endif
