#ifndef ALIFASTMUONTRACKINGRES_H
#define ALIFASTMUONTRACKINGRES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of AliFastResponse for the Muon Spectrometer resolution.
// Author: andreas.morsch@cern.ch
//
#include "AliFastResponse.h"
class AliMUONFastTracking;

class AliFastMuonTrackingRes :  public AliFastResponse {
 public:
    AliFastMuonTrackingRes();
    virtual ~AliFastMuonTrackingRes(){;}
    void SetBackground(Float_t bg = 1.) {fBackground = bg;}
    void SetCharge(Float_t charge = 1.) {fCharge     = charge;}
    virtual void Init();
    virtual void Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			  Float_t& pS,  Float_t& thetaS, Float_t&  phiS);
    virtual Float_t Evaluate(AliFastParticle* part) {
      return AliFastResponse::Evaluate(part);
    }
    virtual Float_t Evaluate(Float_t  pt,  Float_t  theta , Float_t   phi) {
      return AliFastResponse::Evaluate(pt,theta,phi);
    }
 protected:
    Float_t              fBackground;   // Background level
    Float_t              fCharge;       // Current charge
    
    AliMUONFastTracking* fFastTracking; //!Pointer to Fast Tracking Data Handler
    ClassDef(AliFastMuonTrackingRes,1)  // Fast MUON Tracking 
};

#endif


