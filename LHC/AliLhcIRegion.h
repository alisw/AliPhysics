#ifndef ALILHCIREGION_H
#define ALILHCIREGION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Realisation of AliLhcMonitor simulating an LHC interaction region.
// The interaction region is described by the two LHC beams,
// by the beta* and the crossing angle. 
// As a monitor it records the luminosity, average luminosity and beta*
// time evolution.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//
#include <TNamed.h>
#include "AliLhcMonitor.h"

class AliLHC;
class AliLhcBeam;

class AliLhcIRegion : public TNamed, public AliLhcMonitor
{
 public:
    AliLhcIRegion(AliLHC* lhc, const char* name, const char* tile);
    AliLhcIRegion(const AliLhcIRegion &region);
    virtual ~AliLhcIRegion();
    AliLhcIRegion & operator=(const AliLhcIRegion & rhs);
    
    virtual void Init();
    virtual Float_t Luminosity();
    virtual Float_t InitialLumi()   const {return fLuminosity0;}
    virtual Float_t BetaStar()      const {return fBetaStar;}
    virtual void  SetBetaStar(Float_t beta) {fBetaStar = beta;}
    virtual void  SetCrossingAngle(Float_t angle) {fCrossingAngle = angle;}	    
    virtual void  SetAccelerator(AliLHC* acc) {fAccelerator = acc;}
    virtual void  Update();
    virtual void  SetMonitor(Int_t n);
    virtual void  Record();
    virtual void  DrawPlots();

 protected:
    AliLHC*     fAccelerator;         // Accelerator
    AliLhcBeam* fBeam1;               // Beam1
    AliLhcBeam* fBeam2;               // Beam2       
    Float_t     fLuminosity;          // Luminosity
    Float_t     fLuminosity0;         // Initial Luminosity
    Float_t     fAverageLumi;         // Average Initial Luminosity
    Float_t     fBetaStar;            // Beta*
    Float_t     fBetaStar0;           // Initial Beta*
    Float_t     fCrossingAngle;       // Crossing Angle
    Float_t     fFrequency;           // Frequency
    Float_t*    fLumiArray;           // [fNmax] Luminosity(t)
    Float_t*    fAverageLumiArray;    // [fNmax] Average Luminosity(t)
    Float_t*    fBetaStarArray;       // [fNmax] Beta*(t)
//
    ClassDef(AliLhcIRegion,1) // LHC Interaction Region
};

#endif





