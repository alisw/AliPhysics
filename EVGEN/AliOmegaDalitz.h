#ifndef ALIOMEGADALITZ_H
#define ALIOMEGADALITZ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------------
//
// Generate lepton-pair mass distributions for Dalitz decays according
// to the Kroll-Wada parametrization: N. Kroll, W. Wada: Phys. Rev 98(1955)1355
//
// For the electromagnetic form factor the parameterization from
// Lepton-G is used: L.G. Landsberg et al.: Phys. Rep. 128(1985)301
//
//-----------------------------------------------------------------------------


#include "AliDecayer.h"
#include <TLorentzVector.h>

class TH1F;
class TClonesArray;

class AliOmegaDalitz : public AliDecayer
{
 public:
    AliOmegaDalitz();
    virtual void    Init();
    virtual void    Decay(Int_t idpart, TLorentzVector* p);
    virtual Int_t   ImportParticles(TClonesArray *particles);
    virtual void    SetForceDecay(Int_t)                      {;}
    virtual void    ForceDecay()                              {;}
    virtual Float_t GetPartialBranchingRatio(Int_t /*ipart*/) {return -1;}
    virtual Float_t GetLifetime(Int_t /*kf*/)                 {return -1;}
    virtual void    ReadDecayTable()                          {;}
    virtual TH1F*   LeptonPairMassHisto()                     {return  fLPMass;}
	    
 private:
    virtual void    Rot(Double_t pin[3], Double_t pout[3],
			Double_t costheta, Double_t sintheta,
			Double_t cosphi, Double_t sinphi);
    
 protected:
    TH1F*           fLPMass;       // Histogram for lepton pair mass
    TLorentzVector  fProducts[3];  // Decay products
    ClassDef(AliOmegaDalitz, 1) // AliDecayer implementation for omega Dalitz
};

#endif
