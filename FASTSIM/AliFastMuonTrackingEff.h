#ifndef ALIFASTMUONTRACKINGEFF_H
#define ALIFASTMUONTRACKINGEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//
// Class for fast simulation of the ALICE Muon Spectrometer
// Tracking Efficiency.
// The efficiency depends on transverse momentum pt, polar angle theta and azimuthal angle phi.
//
// Author: Alessandro de Falco 
// alessandro.de.falco@ca.infn.it

#include "AliFastResponse.h"
class AliMUONFastTracking;

class AliFastMuonTrackingEff :  public AliFastResponse {
 public:
    AliFastMuonTrackingEff();
    virtual ~AliFastMuonTrackingEff(){;}
    AliFastMuonTrackingEff(const AliFastMuonTrackingEff& eff);
    void SetBackground(Float_t bg = 1.) {fBackground = bg;}
    void SetCharge(Float_t charge = 1.) {fCharge     = charge;}
    virtual void Init();
    virtual Float_t Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi);
    virtual void    Evaluate(Float_t charge, Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
	{AliFastResponse::Evaluate(charge, p, theta, phi, pS, thetaS, phiS);}
    virtual void    Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
	{AliFastResponse::Evaluate(p, theta, phi, pS, thetaS, phiS);}
    
    // Copy
    AliFastMuonTrackingEff& operator=(const AliFastMuonTrackingEff& rhs);
 protected:
    Float_t              fBackground;   // Background level
    Float_t              fCharge;       // Current charge
    
    AliMUONFastTracking* fFastTracking; //!Pointer to Fast Tracking Data Handler
    ClassDef(AliFastMuonTrackingEff,1)  // Fast MUON Tracking Efficiency
};

#endif





