#ifndef ALIFASTMUONTRACKINGACC
#define ALIFASTMUONTRACKINGACC
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliFastResponse.h"
class AliMUONFastTracking;

class AliFastMuonTrackingAcc :  public AliFastResponse {
 public:
    AliFastMuonTrackingAcc();
    ~AliFastMuonTrackingAcc(){;}
    void SetBackground(Float_t bg = 1.) {fBackground = bg;}
    void SetCharge(Float_t charge = 1.) {fCharge     = charge;}
    virtual void Init();
    virtual Float_t Evaluate(Float_t pt, Float_t theta, Float_t phi);
    virtual void    Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS) {
      AliFastResponse::Evaluate(p,  theta,  phi, pS, thetaS, phiS);
    }
    virtual Float_t Evaluate(AliFastParticle* part) {
      return AliFastResponse::Evaluate(part);
    }
 protected:
    Float_t              fBackground;   // Background level
    Float_t              fCharge;       // Current charge
    
    AliMUONFastTracking* fFastTracking; //!Pointer to Fast Tracking Data Handler
    ClassDef(AliFastMuonTrackingAcc,1)  // Fast MUON Tracking Acceptance
};

#endif





