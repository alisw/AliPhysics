#ifndef ALIAODPWG4PARTICLECORRELATION_H
#define ALIAODPWG4PARTICLECORRELATION_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
/// \class AliAODPWG4ParticleCorrelation
/// \brief Daughter of AliAODPWG4Particle that includes correlation part
///
///  AOD objects class in use in the CaloTrackCorrelations
///  analysis pacackge ($ALICE_PHYSICS/PWGGA/CaloTrackCorrelations). 
///  Common format for selected tracks or calorimeter clusters to give as input
///  for different analysis. Basically it contains the particle kinematics
///  and some detailed parameters of the calorimeter cluster and of the intermediate
///  steps of the analysis. 
///
///  Daughter of AliAODPWG4Particle, it includes correlations with respect 
///  this (trigger) object: jets, tracks.
///
///  More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
///  \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//-------------------------------------------------------------------------

//-- ROOT system --
#include "TList.h"
#include "AliAODJet.h"

//-- Analysis system
#include "AliAODPWG4Particle.h"

class AliAODPWG4ParticleCorrelation : public AliAODPWG4Particle {

 public:
  
  AliAODPWG4ParticleCorrelation();
  AliAODPWG4ParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e);
  AliAODPWG4ParticleCorrelation(TLorentzVector & p);  
  AliAODPWG4ParticleCorrelation(AliAODPWG4Particle & p);  
  AliAODPWG4ParticleCorrelation(const AliAODPWG4ParticleCorrelation& p);

  virtual ~AliAODPWG4ParticleCorrelation();
  virtual void   Clear(const Option_t* /*opt*/);
  
  // Deal with the arrays of correlated objects
  virtual TObjArray* GetObjArray(TString refname)      const
  { if(fListOfObjArrays) return (TObjArray*) fListOfObjArrays->FindObject(refname);
    else                 return 0x0 ; }
  
  virtual TList* GetObjArrayList()                     const { return  fListOfObjArrays        ; }
  virtual void   AddObjArray(TObjArray * refarray)           { fListOfObjArrays->Add(refarray) ; }

  // General info on leading particle opposite to trigger
  virtual Int_t  GetLeadingDetector()                  const { return fLeadingDetector ; }
  virtual void   SetLeadingDetector(Int_t d)                 { fLeadingDetector = d    ; }
  
  virtual TLorentzVector  GetLeading()                 const { return  fLeading ; }
  virtual void   SetLeading(TLorentzVector lead)             { fLeading = lead  ; }
  
  // trigger-Jet correlation
  virtual TLorentzVector  GetCorrelatedJet()           const { return  fCorrJet ; }
  virtual void   SetCorrelatedJet(TLorentzVector jet)        { fCorrJet = jet   ; }
  
  virtual TLorentzVector  GetCorrelatedBackground()    const { return  fCorrBkg ; }
  virtual void   SetCorrelatedBackground(TLorentzVector bkg) { fCorrBkg = bkg   ; }
  
  virtual void   SetRefJet(AliAODJet* jet)                   { fRefJet = jet    ; }
  virtual        AliAODJet* GetJet()                   const { return ((AliAODJet*) fRefJet.GetObject()) ; }
  virtual TRef   GetRefJet()                           const { return fRefJet   ; }
  
  virtual void   Print(Option_t* /*option*/) const;
  
 private:
  
  Int_t          fLeadingDetector;  ///< Detector where leading particle was measured.
  TLorentzVector fLeading;          ///< Leading Particle 4-momentum vector on opposite side of trigger particle
  TLorentzVector fCorrJet;          ///< Jet  4-momentum vector
  TLorentzVector fCorrBkg;          ///< Background  4-momentum vector
  TRef           fRefJet;           ///< Reference to jet found with JETAN and correlated with particle
  TList   *      fListOfObjArrays ; ///< List with correlation reference arrays
  
  /// Assignment operator not implemented.
  AliAODPWG4ParticleCorrelation& operator=(const AliAODPWG4ParticleCorrelation& p);
    
  /// \cond CLASSIMP
  ClassDef(AliAODPWG4ParticleCorrelation, 6);
  /// \endcond

};


#endif //ALIAODPWG4PARTICLECORRELATION_H
