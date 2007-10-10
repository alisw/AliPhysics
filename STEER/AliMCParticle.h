#ifndef AliMCParticle_H
#define AliMCParticle_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AliVParticle realisation for MC Particles
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <Rtypes.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRefArray.h>

#include "AliTrackReference.h"
#include "AliVParticle.h"

class AliMCParticle: public AliVParticle {
public:
    AliMCParticle();
    AliMCParticle(TParticle* part, TRefArray* rarray = 0);
    virtual ~AliMCParticle() {}
    AliMCParticle(const AliMCParticle& mcPart); 
    AliMCParticle& operator=(const AliMCParticle& mcPart);
    
    // Kinematics
    virtual Double_t Px()        const;
    virtual Double_t Py()        const;
    virtual Double_t Pz()        const;
    virtual Double_t Pt()        const;
    virtual Double_t P()         const;
    
    virtual Double_t OneOverPt() const;
    virtual Double_t Phi()       const;
    virtual Double_t Theta()     const;
    
    
    virtual Double_t E()         const;
    virtual Double_t M()         const;
    
    virtual Double_t Eta()       const;
    virtual Double_t Y()         const;
    
    virtual Short_t Charge()     const;
    
    // PID
    virtual const Double_t *PID() const {return 0;} // return PID object (to be defined, still)

    // Track References
    Int_t              GetNumberOfTrackReferences() {return fNTrackRef;}
    AliTrackReference* GetTrackReference(Int_t i)
	{return dynamic_cast<AliTrackReference*>((*fTrackReferences)[i]);}
 private:
    TParticle *fParticle;             // The wrapped TParticle
    TRefArray *fTrackReferences;      // Reference array to track references
    Int_t      fNTrackRef;            // Number of track references
    
  ClassDef(AliMCParticle,0)  // AliVParticle realisation for MCParticles
};

inline Double_t AliMCParticle::Px()        const {return fParticle->Px();}
inline Double_t AliMCParticle::Py()        const {return fParticle->Py();}
inline Double_t AliMCParticle::Pz()        const {return fParticle->Pz();}
inline Double_t AliMCParticle::Pt()        const {return fParticle->Pt();}
inline Double_t AliMCParticle::P()         const {return fParticle->P(); }
inline Double_t AliMCParticle::OneOverPt() const {return 1. / fParticle->Pt();}
inline Double_t AliMCParticle::Phi()       const {return fParticle->Phi();}
inline Double_t AliMCParticle::Theta()     const {return fParticle->Theta();}
inline Double_t AliMCParticle::E()         const {return fParticle->Energy();}
inline Double_t AliMCParticle::Eta()       const {return fParticle->Eta();}

inline Double_t AliMCParticle::M()         const
{
    TParticlePDG* pdg = fParticle->GetPDG();
    if (pdg) {
	return (pdg->Mass());
    } else {
	return (fParticle->GetCalcMass());
    }
}


inline Double_t AliMCParticle::Y()         const 
{
    Double_t e  = E();
    Double_t pz = TMath::Abs(Pz());
    if (e != pz) { 
	return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { 
	return -999.;
    }
}

inline Short_t AliMCParticle::Charge()     const
{
    TParticlePDG* pdg = fParticle->GetPDG();
    if (pdg) {
	return (Short_t (pdg->Charge()));
    } else {
	return -99;
    }
}

#endif
