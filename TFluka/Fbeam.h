#ifndef FBEAM_H
#define FBEAM_H 1

#include "cfortran.h"
#include "Rtypes.h"
extern "C" {
//*$ create beam.add
//*copy beam
//*
//*=== beam =============================================================*
//*
//*----------------------------------------------------------------------*
//*     include file: beam copy                    created 26/11/86 by pa*
//*                                                                      *
//*     changes: on 22-oct-1993     by             alfredo ferrari       *
//*                                                                      *
//*     included in the following subroutines or functions: not updated  *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*                                                                      *
//*     /beam/ contains properties of the beam of primary particles      *
//*        pbeam  = average momentum of the beam particles in gev/c      *
//*        dpbeam = momentum spread of the beam in gev/c                 *
//*        divbm  = angular divergense of the beam in mrad               *
//*        xspot  = beam width in x-direction in cm                      *
//*        yspot  = beam width in y-direction in cm                      *
//*        xina   = x-coordinate of the centre of the beam spot          *
//*        yina   = y-coordinate of the centre of the beam spot          *
//*        zina   = z-coordinate of the centre of the beam spot          *
//*        tinx   = direction cosine of the beam with respect to         *
//*                 x-axis                                               *
//*        tiny   = direction cosine of the beam with respect to         *
//*                 y-axis                                               *
//*        tinz   = direction cosine of the beam with respect to         *
//*                 z-axis                                               *
//*        tinpx  = direction cosine of the beam polariz. with respect to*
//*                 x-axis                                               *
//*        tinpy  = direction cosine of the beam polariz. with respect to*
//*                 y-axis                                               *
//*        tinpz  = direction cosine of the beam polariz. with respect to*
//*                 z-axis                                               *
//*        polfra = polarization fraction                                *
//*        nforce = number of the region of forced interaction           *
//*        xfor   = x-coord. of the starting point of the region nforce  *
//*        yfor   = y-coord. of the starting point of the region nforce  *
//*        zfor   = z-coord. of the starting point of the region nforce  *
//*        disfor = thickness of the region nforce in cm                 *
//*        wfor   = relative weight of the particle due to forcing       *
//*        ijbeam = beam particle type (see btype in /paprop/)           *
//*        ijhion = heavy ion type if ijbeam = -2                        *
//*        ipbite = flag describing the shape of the momentum            *
//*                 distribution of the beam                             *
//*                 0=rectangular, 1=gaussian                            *
//*        idiv   = flag describing the shape of the angular             *
//*                 divergence distribution of the beam                  *
//*                 0=rectangular, 1=gaussian                            *
//*        ixspot = flag describing the shape of the spatial             *
//*                 distribution of the beam spot in x-direction         *
//*                 0=rectangular, 1=gaussian                            *
//*        iyspot = flag describing the shape of the spatial             *
//*                 distribution of the beam spot in y-direction         *
//*                 0=rectangular, 1=gaussian                            *
//*        beawei = weight of the beam particles                         *
//*        lbeamc = flag for an annular beam                             *
//*        lpperp = flag for polar. perp. to the beam direction          *
//*        lpfrac = flag for interpreting the polar. fraction            *
//*                                                                      *
//*----------------------------------------------------------------------*
//      LOGICAL LBEAMC, LPPERP, LPFRAC

typedef struct {
   Double_t pbeam;
   Double_t dpbeam;
   Double_t divbm;
   Double_t xspot;
   Double_t yspot;
   Double_t xina;
   Double_t yina;
   Double_t zina;
   Double_t tinx;
   Double_t tiny;
   Double_t tinz;
   Double_t tinpx;
   Double_t tinpy;
   Double_t tinpz;
   Double_t polfra;
   Double_t beawei;
   Double_t xfor;
   Double_t yfor;
   Double_t zfor;
   Double_t disfor;
   Double_t wfor;
   Int_t    ijbeam;
   Int_t    ijhion;
   Int_t    ipbite;
   Int_t    idiv;
   Int_t    ixspot;
   Int_t    iyspot;
   Int_t    nforce;
   Int_t    lbeamc;
   Int_t    lpperp;
   Int_t    lpfrac;
} beamCommon;
#define BEAM COMMON_BLOCK(BEAM,beam)
COMMON_BLOCK_DEF(beamCommon,BEAM);
}

#endif
