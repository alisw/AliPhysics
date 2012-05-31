#ifndef ALIFASTMUONTRACKINGACC_H
#define ALIFASTMUONTRACKINGACC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Realisation of AliFastResponse for the
// fast simulation of the muon spectrometer acceptance.
// The acceptance depends on the muon 3-vector which can be passed as (pt, theta, phi), 
// where pt is the transverse momentum, theta the polar angle and phi the azimuthal angle.
// Author: Andreas Morsch
// andreas.morsch@cern.ch 

#include "AliFastResponse.h"
class AliMUONFastTracking;

class AliFastMuonTrackingAcc :  public AliFastResponse {
 public:
    AliFastMuonTrackingAcc();
    AliFastMuonTrackingAcc(const AliFastMuonTrackingAcc& acc);
    virtual ~AliFastMuonTrackingAcc(){;}
    void SetBackground(Float_t bg = 1.) {fBackground = bg;}
    void SetCharge(Float_t charge = 1.) {fCharge     = charge;}
    virtual void    Init();
    virtual Float_t Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi);
    virtual void    Evaluate(Float_t charge, Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
	{AliFastResponse::Evaluate(charge, p, theta, phi, pS, thetaS, phiS);}
    virtual void    Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
	{AliFastResponse::Evaluate(p, theta, phi, pS, thetaS, phiS);}
    
    // Copy
    AliFastMuonTrackingAcc& operator=(const AliFastMuonTrackingAcc& rhs);
 protected:
    Float_t              fBackground;   // Background level
    Float_t              fCharge;       // Current charge
    
    AliMUONFastTracking* fFastTracking; //!Pointer to Fast Tracking Data Handler
    ClassDef(AliFastMuonTrackingAcc,1)  // Fast MUON Tracking Acceptance
};

#endif





