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
    Double_t WTOPPH [mostck];
    Double_t POPTPH [mostck];
    Double_t XOPTPH [mostck];
    Double_t YOPTPH [mostck];
    Double_t ZOPTPH [mostck];
    Double_t TXOPPH [mostck];
    Double_t TYOPPH [mostck];
    Double_t TZOPPH [mostck];
    Double_t TXPOPP [mostck];
    Double_t TYPOPP [mostck];
    Double_t TZPOPP [mostck];
    Double_t DONEAR [mostck];
    Double_t AGOPPH [mostck];
    Double_t TPROPP [mostck];
    Double_t APROPP [mostck];
    Double_t CMPOPP [mostck];
    Double_t SPAROK [mkbmx1][mostck];
    Int_t    ISPORK [mkbmx2][mostck];
    Int_t    LOOPPH [mostck];
    Int_t    LOUOPP [mostck];
    Int_t    NREGOP [mostck];
    Int_t    NLATOP [mostck];
    Int_t    IPROPP [mostck];
    Int_t    LPROPP [mostck];
    Int_t    NPROPP [mostck];
    Int_t    LSTOPP;
    Int_t    LMXOPP;
} opphstCommon;

#define OPPHST COMMON_BLOCK(OPPHST,opphst)
COMMON_BLOCK_DEF(opphstCommon, OPPHST);

}
#endif

