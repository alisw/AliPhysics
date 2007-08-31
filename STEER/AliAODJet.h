#ifndef AliAODJet_H
#define AliAODJet_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD jet class
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliVParticle.h"
#include <TArrayI.h>
#include "AliAODVertex.h"


class AliAODJet : public AliVParticle {

 public:
    AliAODJet();
    AliAODJet(Double_t px, Double_t py, Double_t pz, Double_t e);
    AliAODJet(TLorentzVector & p);  
    virtual ~AliAODJet();
    AliAODJet(const AliAODJet& jet); 
    AliAODJet& operator=(const AliAODJet& jet);
// AliVParticle methods
    virtual Double_t Px()         const { return fMomentum->Px();      }
    virtual Double_t Py()         const { return fMomentum->Py();      }
    virtual Double_t Pz()         const { return fMomentum->Pz();      }
    virtual Double_t Pt()         const { return fMomentum->Pt();      }
    virtual Double_t P()          const { return fMomentum->P();       }
    virtual Double_t OneOverPt()  const { return 1. / fMomentum->Pt(); }
    virtual Double_t Phi()        const;
    virtual Double_t Theta()      const { return fMomentum->Theta();   }
    virtual Double_t E()          const { return fMomentum->E();       }
    virtual Double_t M()          const { return fMomentum->M();       }
    virtual Double_t Eta()        const { return fMomentum->Eta();     }
    virtual Double_t Y()          const { return fMomentum->Rapidity();}
    

//
    virtual void     AddTrack(TObject *tr) {fRefTracks->Add(tr);}
    TObject* GetTrack(Int_t i) {return fRefTracks->At(i);}
    virtual void     SetBgEnergy(Double_t bgEnCh, Double_t bgEnNe)
	{fBackgEnergy[0] = bgEnCh; fBackgEnergy[1] = bgEnNe;}
    virtual void     SetEffArea(Double_t effACh, Double_t effANe)
	{fEffectiveArea[0] = effACh; fEffectiveArea[1] = effANe;}
    
    virtual TRefArray* GetRefTracks()           const { return  fRefTracks;}
    virtual Double_t   ChargedBgEnergy()        const { return  fBackgEnergy[0];}
    virtual Double_t   NeutralBgEnergy()        const { return  fBackgEnergy[1];}
    virtual Double_t   TotalBgEnergy()          const { return (fBackgEnergy[0] + fBackgEnergy[1]);}

    virtual Double_t   EffectiveAreaCharged()   const { return  fEffectiveArea[0];}
    virtual Double_t   EffectiveAreaNeutral()   const { return  fEffectiveArea[1];}

    virtual void     Print(Option_t* /*option*/) const;
    
// Dummy  
    virtual Short_t Charge()      const { return 0;}
    virtual const Double_t* PID() const { return NULL;}
//
    
    
 private:
    Double32_t      fBackgEnergy[2];     // Subtracted background energy
    Double32_t      fEffectiveArea[2];   // Effective jet area used for background subtraction

    TLorentzVector* fMomentum;           // Jet 4-momentum vector
    TRefArray*      fRefTracks;          // array of references to the tracks belonging to the jet

    ClassDef(AliAODJet,3);

};

inline Double_t AliAODJet::Phi() const
{
    // Return phi
    Double_t phi = fMomentum->Phi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}

#endif
