#ifndef ALIAODPWG4PARTICLECORRELATION_H
#define ALIAODPWG4PARTICLECORRELATION_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODPWG4ParticleCorrelation.h  $ */

//-------------------------------------------------------------------------
//     Copy of AOD photon class, adapted for particle identification
//     and correlations analysis
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --

//-- Analysis system
#include "AliAODJet.h" 
#include "AliAODPWG4Particle.h"

class AliAODPWG4ParticleCorrelation : public AliAODPWG4Particle {

 public:
    AliAODPWG4ParticleCorrelation();
    AliAODPWG4ParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e);
    AliAODPWG4ParticleCorrelation(TLorentzVector & p);  
    AliAODPWG4ParticleCorrelation(AliAODPWG4Particle & p);  

    virtual ~AliAODPWG4ParticleCorrelation();
    AliAODPWG4ParticleCorrelation(const AliAODPWG4ParticleCorrelation& photon); 
    AliAODPWG4ParticleCorrelation& operator=(const AliAODPWG4ParticleCorrelation& photon);

    virtual Float_t IsIsolated() const { return fIsolated ;}
    virtual void    SetIsolated(Bool_t iso) {fIsolated = iso ;}

    virtual TRefArray* GetRefTracks()    const { return  fRefTracks;}
    virtual void       AddTrack(TObject *tr) {fRefTracks->Add(tr);}
    virtual TObject*   GetTrack(Int_t i) const { return fRefTracks->At(i);}

    virtual TRefArray* GetRefClusters()    const { return  fRefClusters;}
    virtual void       AddCluster(TObject *tr) {fRefClusters->Add(tr);}
    virtual TObject*   GetCluster(Int_t i) const { return fRefClusters->At(i);}

    virtual TRefArray* GetRefIsolationConeTracks()  const { return  fRefIsolationConeTracks;}
    virtual void     AddIsolationConeTrack(TObject *tr) {fRefIsolationConeTracks->Add(tr);}
    virtual TObject* GetIsolationConeTrack(Int_t i) const { return fRefIsolationConeTracks->At(i);}

    virtual TRefArray* GetRefIsolationConeClusters()  const { return  fRefIsolationConeClusters;}
    virtual void     AddIsolationConeCluster(TObject *tr) {fRefIsolationConeClusters->Add(tr);}
    virtual TObject* GetIsolationConeCluster(Int_t i) const { return fRefIsolationConeClusters->At(i);}

    virtual TRefArray* GetRefBackgroundTracks()  const { return  fRefBackgroundTracks;}
    virtual void     AddBackgroundTrack(TObject *tr) {fRefBackgroundTracks->Add(tr);}
    virtual TObject* GetBackgroundTrack(Int_t i) const { return fRefBackgroundTracks->At(i);}

    virtual TRefArray* GetRefBackgroundClusters()    const { return  fRefBackgroundClusters;}
    virtual void       AddBackgroundCluster(TObject *tr) {fRefBackgroundClusters->Add(tr);}
    virtual TObject*   GetBackgroundCluster(Int_t i) const { return fRefBackgroundClusters->At(i);}

    virtual TString GetLeadingDetector()   const {return fLeadingDetector ; }
    virtual void SetLeadingDetector(TString d)   {fLeadingDetector = d ; }
    
    virtual TLorentzVector  GetLeading()               const { return  fLeading;}
    virtual void  SetLeading(TLorentzVector lead) {fLeading = lead;}

    virtual TLorentzVector  GetCorrelatedJet()         const { return  fCorrJet;}
    virtual void  SetCorrelatedJet(TLorentzVector jet) {fCorrJet = jet;}

    virtual TLorentzVector  GetCorrelatedBackground()  const { return  fCorrBkg;}
    virtual void  SetCorrelatedBackground(TLorentzVector bkg) {fCorrBkg = bkg;}

    virtual void SetRefJet(AliAODJet* jet)  { fRefJet = jet;}
    virtual      AliAODJet* GetJet() const {return ((AliAODJet*) fRefJet.GetObject());}
    virtual TRef GetRefJet()         const {return fRefJet;}

    virtual void   Print(Option_t* /*option*/) const;

 private:

	Float_t		   fIsolated ; //Particle is isolated or not
	
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

    ClassDef(AliAODPWG4ParticleCorrelation,1);
};


#endif //ALIAODPWG4PARTICLECORRELATION_H
