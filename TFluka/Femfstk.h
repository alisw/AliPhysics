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
const Int_t idmemf = mestck;

typedef struct {
   Double_t e[idmemf]; // total energy in MeV
   Double_t x[idmemf]; // particle x-coordinate
   Double_t y[idmemf]; // particle y-coordinate
   Double_t z[idmemf]; // particle z-coordinate
   Double_t u[idmemf]; // x direction cosine
   Double_t v[idmemf]; // y direction cosine
   Double_t w[idmemf]; // z direction cosine
   Double_t dnear[idmemf]; // equivalent to GEANT "safety"
   Double_t upol[idmemf]; // polarisation in x direction
   Double_t vpol[idmemf]; // polarisation in y direction
   Double_t wpol[idmemf]; // polarisation in z direction
   Double_t usnrml[idmemf];
   Double_t vsnrml[idmemf];
   Double_t wsnrml[idmemf];
   Double_t wt[idmemf]; // weight
   Double_t agemf[idmemf]; // age
   Double_t espark[idmemf][mkbmx1];
   Int_t    iespak[idmemf][mkbmx2];
   Int_t    iq[idmemf]; // charge
   Int_t    ir[idmemf]; // region
   Int_t    irlatt[idmemf]; // lattice cell
   Int_t    nhpemf[idmemf];
   Int_t    lloemf[idmemf]; // generation number
   Int_t    louemf[idmemf];
   Int_t    np; // number of particles in stack
   Int_t    npstrt; // EMF stack index before the interaction (since
                    // the projectile disappears it is also the starting
                    // index of secondaries)
} emfstkCommon;
#define EMFSTK COMMON_BLOCK(EMFSTK,emfstk)
COMMON_BLOCK_DEF(emfstkCommon,EMFSTK);
}
#endif
