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
//*        beawei = weight of the beam particles                         *
//*      bmaxis(j,i) = j_th component of the i_th axis used to define the*
//*                 conventional x,y,z beam reference frame              *
//*!!!!! ATTENTION: in C++ it is the component bmaxis(i,j) !!!!!         *
//*        ijbeam = beam particle type (see btype in /paprop/)           *
//*        ijhion = heavy ion type if ijbeam = -2                        *
//*        ldpgss = true for a gaussian momentum distribution of the     *
//*                 beam particles, false for a rectangular one          *
//*        ldvgss = true for a gaussian angular divergence distribution  *
//*                 of the beam particles, false for a rectangular one   *
//*        ldxgss = true for a gaussian spatial distribution of the beam *
//*                 spot in the x-direction, false for a rectangular one *
//*        ldygss = true for a gaussian spatial distribution of the beam *
//*                 spot in the y-direction, false for a rectangular one *
//*        lbeamc = flag for an annular beam                             *
//*        lpperp = flag for polar. perp. to the beam direction          *
//*        lpfrac = flag for interpreting the polar. fraction
//*        lbaxis = logical flag for using a beam axis frame different   *
//*                 from the standard one                                *
//*                                                                      *
//*----------------------------------------------------------------------*
//      LBEAMC, LPPERP, LPFRAC, LDPGSS, LDVGSS, LDXGSS, LDYGSS,
//     &        LBAXIS

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
   Double_t bmaxis[3][3];
   Int_t    ijbeam;
   Int_t    ijhion;
   Int_t    ldpgss;
   Int_t    ldvgss;
   Int_t    ldxgss;
   Int_t    ldygss;
   Int_t    lbeamc;
   Int_t    lpperp;
   Int_t    lpfrac;
   Int_t    lbaxis;
} beamCommon;
#define BEAM COMMON_BLOCK(BEAM,beam)
COMMON_BLOCK_DEF(beamCommon,BEAM);
}

#endif
