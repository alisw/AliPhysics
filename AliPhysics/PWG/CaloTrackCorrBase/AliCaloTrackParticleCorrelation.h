#ifndef ALICALOTRACKPARTICLECORRELATION_H
#define ALICALOTRACKPARTICLECORRELATION_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
/// \class AliCaloTrackParticleCorrelation
/// \ingroup CaloTrackCorrelationsBase
/// \brief Daughter of AliCaloTrackParticle that includes correlation part
///
///  AOD objects class in use in the CaloTrackCorrelations
///  analysis pacackge ($ALICE_PHYSICS/PWGGA/CaloTrackCorrelations). 
///  Common format for selected tracks or calorimeter clusters to give as input
///  for different analysis. Basically it contains the particle kinematics
///  and some detailed parameters of the calorimeter cluster and of the intermediate
///  steps of the analysis. 
///
///  Daughter of AliCaloTrackParticle, it includes correlations with respect 
///  this (trigger) object: jets, tracks.
///
///  First version in use before september 2017 in $ALICE_ROOT/STEEER/AOD/AliAODPWG4Particle
///
///  More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
///  \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//-------------------------------------------------------------------------

//-- ROOT system --
#include "TList.h"
#include "AliAODJet.h"

//-- Analysis system
#include "AliCaloTrackParticle.h"

class AliCaloTrackParticleCorrelation : public AliCaloTrackParticle {

 public:
  
  AliCaloTrackParticleCorrelation();
  AliCaloTrackParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e);
  AliCaloTrackParticleCorrelation(TLorentzVector & p);  
  AliCaloTrackParticleCorrelation(AliCaloTrackParticle & p);  
  AliCaloTrackParticleCorrelation(const AliCaloTrackParticleCorrelation& p);

  virtual ~AliCaloTrackParticleCorrelation();
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
  AliCaloTrackParticleCorrelation& operator=(const AliCaloTrackParticleCorrelation& p);
    
  /// \cond CLASSIMP
  ClassDef(AliCaloTrackParticleCorrelation, 1);
  /// \endcond

};


#endif //ALICALOTRACKPARTICLECORRELATION_H
