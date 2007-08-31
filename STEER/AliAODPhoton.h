#ifndef AliAODPhoton_H
#define AliAODPhoton_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD photon class
//     Author: Yves Schutz, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliVParticle.h"
#include "AliAODVertex.h"


class AliAODPhoton : public AliVParticle {

 public:
    AliAODPhoton();
    AliAODPhoton(Double_t px, Double_t py, Double_t pz, Double_t e);
    AliAODPhoton(TLorentzVector & p);  
    virtual ~AliAODPhoton();
    AliAODPhoton(const AliAODPhoton& photon); 
    AliAODPhoton& operator=(const AliAODPhoton& photon);
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

    virtual void     Print(Option_t* /*option*/) const;
    
// Dummy  
    virtual Short_t Charge()      const { return 0;}
    virtual const Double_t* PID() const { return NULL;}
//
    
    
 private:
    TLorentzVector* fMomentum;           // Photon 4-momentum vector
    ClassDef(AliAODPhoton,1);
};

inline Double_t AliAODPhoton::Phi() const
{
    // Return phi
    Double_t phi = fMomentum->Phi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}

#endif
