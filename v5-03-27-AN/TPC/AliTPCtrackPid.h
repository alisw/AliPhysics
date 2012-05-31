#ifndef ALITPCTRACKPID_H
#define ALITPCTRACKPID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class to determine the particle ID for TPC tracks
//

#include <TMath.h>
#include <TObject.h>
#include "Riostream.h"

//_____________________________________________________________________________
class AliTPCtrackPid : public TObject {
public:
    AliTPCtrackPid();
    virtual ~AliTPCtrackPid(){}

private:
    Float_t fWpi;     // Particle ID
    Float_t fWk;      // Particle k
    Float_t fWp;      // Particle momentum
    Float_t fSignal;  // Signal
    Float_t fMom;     // Momentum
    Float_t fPhi;     // Phi
    Float_t fLam;     // Lambda
    Int_t   fPcode;   // Particle code
    Int_t   fLabel;   // Particle Label
    
    Float_t fGSignal; // Signal G
    Float_t fGMom;    // Moment G
    Float_t fGpx;     // px G
    Float_t fGpy;     // py G
    Float_t fGpz;     // pz G
    Float_t fGx;      // x G
    Float_t fGy;      // y G
    Float_t fGz;      // z G
    Int_t   fGcode;   // Code G
    Int_t   fGlab;    // Lab G
    //

  ClassDef(AliTPCtrackPid,1)  // TPC track PID
};

#endif


