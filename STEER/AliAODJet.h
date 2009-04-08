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
    virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }
    virtual Double_t Phi()        const;
    virtual Double_t Theta()      const { return fMomentum->Theta();   }
    virtual Double_t E()          const { return fMomentum->E();       }
    virtual Double_t M()          const { return fMomentum->M();       }
    virtual Double_t Eta()        const { return fMomentum->Eta();     }
    virtual Double_t Y()          const { return fMomentum->Rapidity();}
    virtual Double_t Xv()         const {return -999.;} // put reasonable values here
    virtual Double_t Yv()         const {return -999.;} //
    virtual Double_t Zv()         const {return -999.;} //
    virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }  

//
    virtual void     AddTrack(TObject *tr) {fRefTracks->Add(tr);}
    TObject* GetTrack(Int_t i) {return fRefTracks->At(i);}
    virtual void     SetBgEnergy(Double_t bgEnCh, Double_t bgEnNe)
	{fBackgEnergy[0] = bgEnCh; fBackgEnergy[1] = bgEnNe;}
    virtual void     SetEffArea(Double_t effACh, Double_t effANe)
	{fEffectiveArea[0] = effACh; fEffectiveArea[1] = effANe;}
    virtual void     SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e);

    virtual TRefArray* GetRefTracks()           const { return  fRefTracks;}
    virtual Double_t   ChargedBgEnergy()        const { return  fBackgEnergy[0];}
    virtual Double_t   NeutralBgEnergy()        const { return  fBackgEnergy[1];}
    virtual Double_t   TotalBgEnergy()          const { return (fBackgEnergy[0] + fBackgEnergy[1]);}

    virtual Double_t   EffectiveAreaCharged()   const { return  fEffectiveArea[0];}
    virtual Double_t   EffectiveAreaNeutral()   const { return  fEffectiveArea[1];}
    virtual Double_t   DeltaR(const AliVParticle* part);

    
    TLorentzVector*    MomentumVector()         const {return fMomentum;}
    virtual void       Print(Option_t* /*option*/) const;
    
// Dummy  
    virtual Short_t Charge()      const { return 0;}
    virtual const Double_t* PID() const { return NULL;}
    virtual Int_t   GetLabel()    const { return -1;}
//

    /** Compare this class with an other instance of this class
     *  used in a TClonesArray::Sort()
     *  @param   obj  ptr to other instance
     *  @return  Returns 0 when equal, 1 when this is smaller
     *  and -1 when bigger -- sorts descending
     */
    Int_t Compare( const TObject* obj) const;
    
    
    /** Defines this class as being sortable in a TClonesArray
     *  @return     always kTRUE;
     */
    Bool_t IsSortable() const  { return kTRUE; }

 private:
    Double32_t      fBackgEnergy[2];     // Subtracted background energy
    Double32_t      fEffectiveArea[2];   // Effective jet area used for background subtraction

    TLorentzVector* fMomentum;           // Jet 4-momentum vector
    TRefArray*      fRefTracks;          // array of references to the tracks belonging to the jet

    ClassDef(AliAODJet,4);

};

inline Double_t AliAODJet::Phi() const
{
    // Return phi
    Double_t phi = fMomentum->Phi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}

#endif
