#ifndef AliAODParticleCorrelation_H
#define AliAODParticleCorrelation_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODParticleCorrelation.h  $ */

//-------------------------------------------------------------------------
//     Copy of AOD photon class, adapted for particle identification
//     and correlations analysis
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --
#include <TLorentzVector.h>
class TString;

//-- Analysis system
#include "AliAODJet.h"
#include "AliVParticle.h"

class AliAODParticleCorrelation : public AliVParticle {

 public:
    AliAODParticleCorrelation();
    AliAODParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e);
    AliAODParticleCorrelation(TLorentzVector & p);  
    virtual ~AliAODParticleCorrelation();
    AliAODParticleCorrelation(const AliAODParticleCorrelation& photon); 
    AliAODParticleCorrelation& operator=(const AliAODParticleCorrelation& photon);

    // AliVParticle methods
    virtual Double_t Px()         const { return fMomentum->Px();      }
    virtual Double_t Py()         const { return fMomentum->Py();      }
    virtual Double_t Pz()         const { return fMomentum->Pz();      }
    virtual Double_t Pt()         const { return fMomentum->Pt();      }
    virtual Double_t P()          const { return fMomentum->P();       }
    virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }
    virtual Double_t OneOverPt()  const { return 1. / fMomentum->Pt(); }
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
    virtual void     Print(Option_t* /*option*/) const;

    //
    //Dummy
    virtual Short_t Charge()      const { return 0;}
    virtual const Double_t* PID() const { return NULL;}
    //

    virtual Float_t GetR()   const {return fR ; }
    virtual Int_t GetPdg()   const {return fPdg ; }
    virtual Int_t GetTag()   const {return fTag ; }
    virtual Int_t GetLabel()   const {return fLabel ; }
    virtual TString GetDetector()   const {return fDetector ; }

    virtual void SetR(Float_t r)   {fR = r ; }
    virtual void SetPdg(Int_t pdg)   {fPdg = pdg ; }
    virtual void SetTag(Int_t tag)   {fTag = tag ; }
    virtual void SetLabel(Int_t l)   {fLabel = l ; }
    virtual void SetDetector(TString d)   {fDetector = d ; }

    virtual TRefArray* GetRefTracks()           const { return  fRefTracks;}
    virtual void     AddTrack(TObject *tr) {fRefTracks->Add(tr);}
    TObject* GetTrack(Int_t i) {return fRefTracks->At(i);}

    virtual TRefArray* GetRefClusters()           const { return  fRefClusters;}
    virtual void     AddCluster(TObject *tr) {fRefClusters->Add(tr);}
    TObject* GetCluster(Int_t i) {return fRefClusters->At(i);}

    virtual TRefArray* GetRefIsolationConeTracks()           const { return  fRefIsolationConeTracks;}
    virtual void     AddIsolationConeTrack(TObject *tr) {fRefIsolationConeTracks->Add(tr);}
    TObject* GetIsolationConeTrack(Int_t i) {return fRefIsolationConeTracks->At(i);}

    virtual TRefArray* GetRefIsolationConeClusters()           const { return  fRefIsolationConeClusters;}
    virtual void     AddIsolationConeCluster(TObject *tr) {fRefIsolationConeClusters->Add(tr);}
    TObject* GetIsolationConeCluster(Int_t i) {return fRefIsolationConeClusters->At(i);}

    virtual TRefArray* GetRefBackgroundTracks()           const { return  fRefBackgroundTracks;}
    virtual void     AddBackgroundTrack(TObject *tr) {fRefBackgroundTracks->Add(tr);}
    TObject* GetBackgroundTrack(Int_t i) {return fRefBackgroundTracks->At(i);}

    virtual TRefArray* GetRefBackgroundClusters()           const { return  fRefBackgroundClusters;}
    virtual void     AddBackgroundCluster(TObject *tr) {fRefBackgroundClusters->Add(tr);}
    TObject* GetBackgroundCluster(Int_t i) {return fRefBackgroundClusters->At(i);}

    virtual void SetLeadingDetector(TString d)   {fLeadingDetector = d ; }
    virtual TString GetLeadingDetector()   const {return fLeadingDetector ; }
    
    virtual TLorentzVector  GetLeading()           const { return  fLeading;}
    virtual void  SetLeading(TLorentzVector lead) {fLeading = lead;}

    virtual TLorentzVector  GetCorrelatedJet()           const { return  fCorrJet;}
    virtual void  SetCorrelatedJet(TLorentzVector jet) {fCorrJet = jet;}

    virtual TLorentzVector  GetCorrelatedBackground()           const { return  fCorrBkg;}
    virtual void  SetCorrelatedBackground(TLorentzVector bkg) {fCorrBkg = bkg;}

    void SetRefJet(AliAODJet* jet)  { fRefJet = jet;}
    AliAODJet* GetJet() {return ((AliAODJet*) fRefJet.GetObject());}
    TRef GetRefJet() {return fRefJet;}

 private:
    TLorentzVector* fMomentum;           // Photon 4-momentum vector
    Int_t           fPdg; // id of particle
    Int_t           fTag; // tag of particle (decay, fragment, prompt photon)
    Int_t           fLabel; // MC label
    TString         fDetector; // Detector where particle was measured.
    Float_t         fR ; // Isolation cone size
    TRefArray*     fRefTracks;  // array of references to the tracks belonging to the jet / all selected hadrons  
    TRefArray*     fRefClusters; // array of references to the clusters belonging to the jet / all selected hadrons  
    
    TRefArray*     fRefIsolationConeTracks;  // array of references to the tracks belonging to the cone around direct particle candidate  
    TRefArray*     fRefIsolationConeClusters; // array of references to the clusters belonging to the  cone around direct particle candidate  

    TRefArray*     fRefBackgroundTracks;  // array of references to the tracks for background stimation
    TRefArray*     fRefBackgroundClusters; // array of references to the clusters for background stimation 

    TString        fLeadingDetector; // Detector where leading particle was measured.
    TLorentzVector fLeading;     // Leading Particle 4-momentum vector
    
    TLorentzVector fCorrJet;     // Jet  4-momentum vector
    TLorentzVector fCorrBkg;     // Background  4-momentum vector

    TRef           fRefJet; // Rerence to jet found with JETAN and correlated with particle

    ClassDef(AliAODParticleCorrelation,1);
};

inline Double_t AliAODParticleCorrelation::Phi() const
{
  // Return phi
  Double_t phi = fMomentum->Phi();
  if (phi < 0.) phi += 2. * TMath::Pi();
  return phi;
}

#endif
