#ifndef ALIFASTMUONTRIGGEREFF_H
#define ALIFASTMUONTRIGGEREFF_H
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <AliFastResponse.h>
enum CutTupe {kLow, kHigh};

class AliFastMuonTriggerEff : public AliFastResponse {
    
 public:
    AliFastMuonTriggerEff();
    AliFastMuonTriggerEff(char* Name, char* Title):AliFastResponse(Name, Title) {;}    
    virtual ~AliFastMuonTriggerEff(){;}
    virtual void    Init();
    virtual void    Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi,
			     Float_t& effLow, Float_t& effHigh, Float_t& eff);
    virtual Float_t Evaluate(Float_t charge, Float_t pt, Float_t theta, Float_t phi);
    virtual void    Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			     Float_t& pS,  Float_t& thetaS, Float_t&  phiS) {
      AliFastResponse::Evaluate(p, theta , phi, pS, thetaS, phiS);
    }
    virtual Float_t Evaluate(Float_t pt, Float_t theta, Float_t phi) {
      return AliFastResponse::Evaluate(pt, theta, phi);
    }
    virtual Float_t Evaluate(AliFastParticle* part) {
      return AliFastResponse::Evaluate(part);
    }
    virtual void    SetCut(Int_t cut = kLow) {fCut = cut;}
    virtual Float_t Cut() {return fCut;}
  protected:
    virtual void InitTree();
  protected:
    Int_t fLook[2][10][20];       // Look up table for bkg=0
    Float_t fDpt;                 // Delta_pt
    Float_t fPhiMin;              // lower limit for phi 
    Float_t fPhiMax;              // upper limit for phi
    Float_t fDphi;                // Delta_phi
    Float_t fThetaMin;            // lower limit for theta
    Float_t fThetaMax;            // upper limit for theta
    Float_t fDtheta;              // Delta_theta
    Int_t   fCut;                 // Cut type (low/high)
    Int_t   fZones;               // Total number of zones
    static const Int_t   fSim=2;  // Number of pt extentions (internal use)
    Float_t** fEffLow;            // Table for low-pt  cut bkg=0
    Float_t** fEffHigh;           // Table for high-pt cut bkg=0
    
    ClassDef(AliFastMuonTriggerEff,1)    // Fast Muon Trigger response
};

#endif 



