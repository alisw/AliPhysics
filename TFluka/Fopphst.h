#ifndef OPPHST
#define OPPHST_H 1

#include "Rtypes.h"
#include "cfortran.h"
#include "Fdimpar.h"

extern "C" {

//*$ CREATE OPPHST.ADD
//*COPY OPPHST
//*
//*=== Opphst ===========================================================//*
//*
//*----------------------------------------------------------------------//*
//*                                                                      //*
//*     OPtical PHoton STack:                                            //*
//*                                                                      //*
//*     Created on 19 september 1997 by    Alfredo Ferrari & Paola Sala  //*
//*                                                   Infn - Milan       //*
//*                                                                      //*
//*     Last change on 13-oct-98     by    Alfredo Ferrari               //*
//*                                                                      //*
//*        wtopph = weight of the photon                                 //*
//*        poptph = laboratory momentum of the photon in GeV/c           //*
//*        xoptph = x-coordinate of the photon                           //*
//*        yoptph = y-coordinate of the photon                           //*
//*        zoptph = z-coordinate of the photon                           //*
//*        txopph = direction cosine of the photon                       //*
//*                 with respect to x-axis                               //*
//*        tyopph = direction cosine of the photon                       //*
//*                 with respect to y-axis                               //*
//*        tzopph = direction cosine of the photon                       //*
//*                 with respect to z-axis                               //*
//*        txpopp = direction cosine of the photon polarization          //*
//*        typopp = direction cosine of the photon polarization          //*
//*        tzpopp = direction cosine of the photon polarization          //*
//*        donear = distance to the nearest boundary                     //*
//*        agopph = age of the photon (seconds)                          //*
//*        cmpopp = total path length of the photon (cm)                 //*
//*        loopph = generation of the photon                             //*
//*        louopp = user flag                                            //*
//*        nregop = number of the region of the photon                   //*
//*        nlatop = number of the lattice cell of the photon             //*
//*        tpropp = kinetic energy of parent particle of the photon      //*
//*        apropp = age of the parent particle of the photon (seconds)   //*
//*        ipropp = id (paprop) of the parent particle of the photon     //*
//*        lpropp = generation of the parent particle of the photon      //*
//*        npropp = # of the primary track which generated the photon    //*
//*                 (not used for the moment)                            //*
//*        lstopp = stack pointer                                        //*
//*        lmxopp = highest value of the stack pointer encountered       //*
//*                 in the run                                           //*
//*                                                                      //*
//*----------------------------------------------------------------------//*
//*
typedef struct {
    Double_t wtopph [mostck];
    Double_t poptph [mostck];
    Double_t xoptph [mostck];
    Double_t yoptph [mostck];
    Double_t zoptph [mostck];
    Double_t txopph [mostck];
    Double_t tyopph [mostck];
    Double_t tzopph [mostck];
    Double_t txpopp [mostck];
    Double_t typopp [mostck];
    Double_t tzpopp [mostck];
    Double_t donear [mostck];
    Double_t agopph [mostck];
    Double_t tpropp [mostck];
    Double_t apropp [mostck];
    Double_t cmpopp [mostck];
    Double_t sparok [mostck][mkbmx1];
    Int_t    ispork [mostck][mkbmx2];
    Int_t    loopph [mostck];
    Int_t    louopp [mostck];
    Int_t    nregop [mostck];
    Int_t    nlatop [mostck];
    Int_t    ipropp [mostck];
    Int_t    lpropp [mostck];
    Int_t    npropp [mostck];
    Int_t    lstopp;
    Int_t    lmxopp;
} opphstCommon;

#define OPPHST COMMON_BLOCK(OPPHST,opphst)
COMMON_BLOCK_DEF(opphstCommon, OPPHST);

}
#endif

