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
#include "TList.h"

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
  
  virtual TRefArray* GetRefArray(TString refname)  const { return   (TRefArray*) fListOfRefArrays->FindObject(refname); } 
  virtual TList*     GetRefArrayList()             const { return  fListOfRefArrays; } 
  virtual void       AddRefArray(TRefArray * refarray)  { fListOfRefArrays->Add(refarray); }

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
  
  virtual Bool_t IsIsolated() const      { return fIsolated ; }
  virtual void   SetIsolated(Bool_t iso) { fIsolated = iso ; }

  virtual void   Print(Option_t* /*option*/) const;
  
 private:
  Bool_t         fIsolated ;        //Particle is isolated or not  
  TString        fLeadingDetector;  // Detector where leading particle was measured.
  TLorentzVector fLeading;          // Leading Particle 4-momentum vector
  TLorentzVector fCorrJet;          // Jet  4-momentum vector
  TLorentzVector fCorrBkg;          // Background  4-momentum vector
  TRef           fRefJet;           // Reference to jet found with JETAN and correlated with particle
  TList   *      fListOfRefArrays ; // List with correlation reference arrays
  
  ClassDef(AliAODPWG4ParticleCorrelation, 2);
};


#endif //ALIAODPWG4PARTICLECORRELATION_H
