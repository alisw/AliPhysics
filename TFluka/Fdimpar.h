#ifndef FDIMPAR_H
#define FDIMPAR_H

#include "cfortran.h"
#include "Rtypes.h"
extern "C" {
//*$ create dimpar.add
//*copy dimpar
//*                                                                      *
//*=== dimpar ===========================================================*
//*                                                                      *
//*----------------------------------------------------------------------*
//*                                                                      *
//*      dimpar: included in any routine                                 *
//*                                                                      *
//*          mxxrgn = maximum number of regions                          *
//*          mxxmdf = maximum number of media in fluka                   *
//*          mxxmde = maximum number of media in emf                     *
//*          mfstck = stack dimension in fluka                           *
//*          mestck = stack dimension in emf                             *
//*          mostck = stack dimension for optical photons                *
//*          mxprsn = secondary stack dimension for resonance generator  *
//*          mxpdpm = secondary stack dimension for dpm generators       *
//*          Mxpscs = secondary stack dimension overall                  *
//*          mxoutu = maximum number of output units                     *
//*          nallwp = number of allowed particles                        *
//*          nelemx = number of maximum allowed elements of a compound   *
//*          mpdpdx = number of particle types for which em de/dx pro-   *
//*                   cesses (ion,pair,bremss) have to be computed       *
//*          mxhttr = maximum number of (hit) target nucleons for a      *
//*                   given collision generation                         *
//*          icomax = maximum number of materials for compounds (equal   *
//*                   to the sum of the number of materials for every    *
//*                   compound )                                         *
//*          ichmax = maximum number of harmonic oscillator levels for   *
//*                   compounds (equal to the sum of the number of har-  *
//*                   monic oscillator levels for every compound )       *
//*          nstbis = number of stable isotopes recorded in common iso-  *
//*                   top                                                *
//*          mxpabl = number of resonances inside hadrin common blocks   *
//*          idmaxp = number of particles/resonances defined in common   *
//*                   part                                               *
//*          idmxdc = number of particles/resonances decay channels      *
//*                   defined in common decayc                           *
//*          ihypmx = maximum number of hyperons in a hypernucleus       *
//*          mkbmx1 = dimension for kwb real spare array in fluka stack  *
//*          mkbmx2 = dimension for kwb int. spare array in fluka stack  *
//*                                                                      *
//*----------------------------------------------------------------------*
//*                                                                      *
const Int_t mxxrgn = 10000;
const Int_t mxxmdf = 510;
const Int_t mxxmde = 502;
const Int_t mfstck = 5500;
const Int_t mestck = 100;
const Int_t mostck = 2000;
const Int_t mxprsn = 100;
const Int_t mxpdpm = 800;
const Int_t mxpscs = 3999;
const Int_t mxoutu = 50;
const Int_t nallwp = 64;
const Int_t nelemx = 80;
const Int_t mpdpdx = 18;
const Int_t mxhttr = 260;
const Int_t icomax = 700;
const Int_t ichmax = icomax+mxxmdf;
const Int_t nstbis = 304;
//* till 3-aug-99:
//*     const Int_t mxpabl =  110;
const Int_t mxpabl = 120;
const Int_t idmaxp = 450;
const Int_t idmxdc = 2000;
const Int_t mxmcin = 410;
const Int_t ihypmx = 4;
//* till 19-jul-2000:
//*     const Int_t mkbmx1 =    9;
//*     const Int_t mkbmx2 =    3;
const Int_t mkbmx1 = 11;
const Int_t mkbmx2 = 11;
}

#endif
