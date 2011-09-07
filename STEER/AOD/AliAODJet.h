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
    virtual Bool_t   IsTriggeredEMCAL(){return (fTrigger&kEMCALTriggered)==kEMCALTriggered;}
    virtual Bool_t   IsTriggeredTRD(){return (fTrigger&kTRDTriggered)==kTRDTriggered;}
    virtual UChar_t  Trigger(){return fTrigger;}

    virtual void     AddTrack(TObject *tr);
    
    TObject* GetTrack(Int_t i) {return fRefTracks->At(i);}
    virtual void     SetBgEnergy(Double_t bgEnCh, Double_t bgEnNe)
	{fBackgEnergy[0] = bgEnCh; fBackgEnergy[1] = bgEnNe;}
    virtual void     SetEffArea(Double_t effACh, Double_t effANe, Double_t effAErrCh = 0, Double_t effAErrNe = 0)
	{
	  fEffectiveArea[0] = effACh; fEffectiveArea[1] = effANe;
	  fEffectiveAreaError[0] = effAErrCh;
	  fEffectiveAreaError[1] = effAErrNe;
	}
    virtual void     SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e);
    virtual void     SetPtEtaPhiM(Double_t pt, Double_t eta, Double_t phi, Double_t m);
    virtual void     SetTrigger(UChar_t f){fTrigger |= f;}
    virtual void     ResetTrigger(UChar_t f){fTrigger &= ~f;}
    virtual void     SetNEF(Double_t nef) {fNeutralFraction=nef;}
    virtual Double_t GetNEF() const {return fNeutralFraction;}

    virtual TRefArray* GetRefTracks()           const { return  fRefTracks;}
    virtual Double_t   ChargedBgEnergy()        const { return  fBackgEnergy[0];}
    virtual Double_t   NeutralBgEnergy()        const { return  fBackgEnergy[1];}
    virtual Double_t   TotalBgEnergy()          const { return (fBackgEnergy[0] + fBackgEnergy[1]);}

    virtual Double_t   EffectiveAreaCharged()   const { return  fEffectiveArea[0];}
    virtual Double_t   EffectiveAreaNeutral()   const { return  fEffectiveArea[1];}
    virtual void SetVectorAreaCharged(TLorentzVector *effVACh){
      if(!fVectorAreaCharged)fVectorAreaCharged= new TLorentzVector(*effVACh);
      else *fVectorAreaCharged = *effVACh;
    }
    virtual TLorentzVector*  VectorAreaCharged()   const {return fVectorAreaCharged;}



    virtual Double_t   ErrorEffectiveAreaCharged()   const { return  fEffectiveAreaError[0];}
    virtual Double_t   ErrorEffectiveAreaNeutral()   const { return  fEffectiveAreaError[1];}
    virtual Double_t   DeltaR(const AliVParticle* part);

    
    TLorentzVector*    MomentumVector()         const {return fMomentum;}
    virtual void       Print(Option_t* /*option*/) const;
    
// Dummy  
    virtual Short_t Charge()      const { return 0;}
    virtual const Double_t* PID() const { return NULL;}
    virtual Int_t   GetLabel()    const { return -1;}
  // Dummy
    virtual Int_t    PdgCode()    const {return 0;}

//

    // first only one bit for EMCAL and TRD, leave space for more
    // trigger types and/or other detectors
    enum {kEMCALTriggered = 1<<0,
	  kTRDTriggered =   1<<2,
	  kHighTrackPtTriggered = 1<<7};


 private:
    Double32_t      fBackgEnergy[2];         // Subtracted background energy
    Double32_t      fEffectiveArea[2];       // Effective jet area used for background subtraction
    Double32_t      fEffectiveAreaError[2];  //[0,1,10] relative error of jet areas, 10 bit precision
    Double32_t      fNeutralFraction;        //[0,1,12] Neutral fraction between 0 and 1 12 bit precision;
    UChar_t         fTrigger;                // Bit mask to flag jets triggered by a certain detector  
    TLorentzVector* fMomentum;               // Jet 4-momentum vector
    TLorentzVector* fVectorAreaCharged;      // jet area four momentum 
    TRefArray*      fRefTracks;              // array of references to the tracks belonging to the jet

    ClassDef(AliAODJet,9);

};

inline Double_t AliAODJet::Phi() const
{
    // Return phi
    Double_t phi = fMomentum->Phi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}

#endif
