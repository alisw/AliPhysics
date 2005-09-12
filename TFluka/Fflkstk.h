#ifndef FFLKSTK_H
#define FFLKSTK_H 1

#include "cfortran.h"
#include "Rtypes.h"

#include "Fdimpar.h"

extern "C" {
//*$ create flkstk.add
//*copy flkstk
//*
//*=== flkstk ============================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     FLUKA90-200x particle stack:                                     *
//*                                                                      *
//*     Changes: last change on 15-may-2005   by    Alfredo Ferrari      *
//*                                                  INFN, Milan         *
//*                                                                      *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*     /Flkstk/ stack for the primaries                                  *
//*        Wtflk  = particle statistical weight                          *
//*        Pmoflk = particle (laboratory) momentum (GeV/c)               *
//*        Tkeflk = particle (laboratory) kinetic energy (GeV)           *
//*        Xflk   = particle position  x-coordinate                      *
//*        Yflk   = particle position  y-coordinate                      *
//*        Zflk   = particle position  z-coordinate                      *
//*        Txflk  = particle direction x-coordinate                      *
//*        Tyflk  = particle direction y-coordinate                      *
//*        Tzflk  = particle direction z-coordinate                      *
//*        Txpol  = x direction cosine of the particle polarization      *
//*        Typol  = y direction cosine of the particle polarization      *
//*        Tzpol  = z direction cosine of the particle polarization      *
//*        Txnor  = x direction cosine of a (possible) surface normal    *
//*        Tynor  = y direction cosine of a (possible) surface normal    *
//*        Tznor  = z direction cosine of a (possible) surface normal    *
//*        Dfnear = distance to the nearest boundary                     *
//*        Agestk = age of the particle (seconds)                        *
//*        Aknshr = Kshort component of K0/K0bar                         *
//*        Frcphn = cross section for force photonuclear interaction (if *
//*                 < 0), distance to a forced photonuclear interaction  *
//*                (if > 0)                                              *
//*        Lfrphn = flag for forced photonuclear interaction             *
//*        Raddly = delay (s) in production wrt the nominal primary "0"  *
//*                 time for particle produced in radioactive decays     *
//*                (i.e. those coming from decays of daughter isotopes), *
//*                 when in analogue mode, flag for position in the      *
//*                 activr array when in non-analogue mode               *
//*        Cmpath = cumulative path travelled by the particle since it   *
//*                 was produced (cm)                                    *
//*        Sparek = spare real variables available for K.W.Burn          *
//*        Ispark = spare integer variables available for K.W.Burn       *
//*        Iloflk = particle identity (Paprop numbering)                 *
//*        Igroup = energy group for low energy neutrons                 *
//*        Loflk  = particle generation                                  *
//*        Louse  = user flag                                            *
//*        Nrgflk = particle region number                               *
//*        Nlattc = particle lattice cell number                         *
//*        Nhspnt = pointer to the history object (Geant4 geometry)      *
//*        Nevent = number of the event which created the particle       *
//*        Numpar = particle number                                      *
//*        Lraddc = flag for particles generated in radioactive decays   *
//*        Nparma = largest particle number ever reached                 *
//*        Nstmax = highest value of the stack pointer ever reached      *
//*                 in the run                                           *
//*        Npflka = Fluka stack pointer                                  *
//*        Nstaol = stack pointer of the last processed particle         *
//*        Igroun = energy group number of the last processed particle   *
//*                 if it is a low energy neutron                        *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Double_t xflk[mfstck+1];           //(0:MFSTCK)
   Double_t yflk[mfstck+1];           //(0:MFSTCK)
   Double_t zflk[mfstck+1];           //(0:MFSTCK)
   Double_t txflk[mfstck+1];          //(0:MFSTCK)
   Double_t tyflk[mfstck+1];          //(0:MFSTCK)
   Double_t tzflk[mfstck+1];          //(0:MFSTCK)
   Double_t txpol[mfstck+1];          //(0:MFSTCK)
   Double_t typol[mfstck+1];          //(0:MFSTCK)
   Double_t tzpol[mfstck+1];          //(0:MFSTCK)
   Double_t txnor[mfstck+1];          //(0:MFSTCK)
   Double_t tynor[mfstck+1];          //(0:MFSTCK)
   Double_t tznor[mfstck+1];          //(0:MFSTCK)
   Double_t wtflk[mfstck+1];          //(0:MFSTCK)
   Double_t pmoflk[mfstck+1];         //(0:MFSTCK)
   Double_t tkeflk[mfstck+1];         //(0:MFSTCK)
   Double_t dfnear[mfstck+1];         //(0:MFSTCK)
   Double_t agestk[mfstck+1];         //(0:MFSTCK)
   Double_t aknshr[mfstck+1];         //(0:MFSTCK)
   Double_t raddly[mfstck+1];         //(0:MFSTCK)
   Double_t cmpath[mfstck+1];         //(0:MFSTCK)
   Double_t frcphn[mfstck+1];         //(0:MFSTCK)
   Double_t sparek[mfstck+1][mkbmx1]; //(MKBMX1,0:MFSTCK)
   Int_t    ispark[mfstck+1][mkbmx2]; //(MKBMX2,0:MFSTCK)
   Int_t    iloflk[mfstck+1];         //(0:MFSTCK)
   Int_t    igroup[mfstck+1];         //(0:MFSTCK)
   Int_t    loflk[mfstck+1];          //(0:MFSTCK)
   Int_t    louse[mfstck+1];          //(0:MFSTCK)
   Int_t    nrgflk[mfstck+1];         //(0:MFSTCK)
   Int_t    nlattc[mfstck+1];         //(0:MFSTCK)
   Int_t    nhspnt[mfstck+1];         //(0:MFSTCK)
   Int_t    nevent[mfstck+1];         //(0:MFSTCK)
   Int_t    numpar[mfstck+1];         //(0:MFSTCK)
   Int_t    lraddc[mfstck+1];         //(0:MFSTCK)
   Int_t    lfrphn[mfstck+1];         //(0:MFSTCK)
   Int_t    nparma;
   Int_t    nstmax;
   Int_t    npflka;
   Int_t    nstaol;
   Int_t    igroun;
} flkstkCommon;
#define FLKSTK COMMON_BLOCK(FLKSTK,flkstk)
COMMON_BLOCK_DEF(flkstkCommon,FLKSTK);
}

#endif
