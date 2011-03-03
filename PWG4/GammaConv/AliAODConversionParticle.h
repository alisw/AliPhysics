#ifndef ALIAODCONVERSIONPARTICLE_H
#define ALIAODCONVERSIONPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class containing the aod information from conversions
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---

class AliStack;
class AliESDEvent;
#include "TMath.h"
#include "AliAODPhoton.h"
#include "AliGammaConversionAODObject.h"
#include "AliKFParticle.h"

class AliAODConversionParticle : public AliAODPhoton {

 public: 

  //Constructors
  AliAODConversionParticle();    
  AliAODConversionParticle(TLorentzVector& momentum); 
  AliAODConversionParticle(AliGammaConversionAODObject *obj);
  AliAODConversionParticle(AliKFParticle * gammakf, Int_t label1, Int_t label2);

  //Constructor Mother Particle
  AliAODConversionParticle(AliAODConversionParticle *y1,AliAODConversionParticle *y2);

  //Copy Constructor
  AliAODConversionParticle(const AliAODConversionParticle & g);           
  //assignment operator
  AliAODConversionParticle & operator = (const AliAODConversionParticle & g);

  //Destructor
  virtual ~AliAODConversionParticle() {;}

  ///Set the Chi2 of reconstructed conversion gamma
  void SetChi2(Float_t chi2) {fChi2 = chi2;}

  ///Set track or MC labels
  void SetLabel1(Int_t label){fLabel[0] = label;}
  void SetLabel2(Int_t label){fLabel[1] = label;}
  void SetTrackLabels(Int_t label1, Int_t label2){fLabel[0] = label1; fLabel[1] = label2;}
  
  ///Set Invariant mass
  void SetIMass(Float_t im) { fIMass = im; }

  ///Set the tag for decay meson
  void SetTag( Bool_t tagged ) { fTagged = tagged; }

  ///Set pointer to MC stack
  void SetStack(AliStack* const stack){fMCStack=stack;}
  
  ///Set pointer to ESD event
  void SetESDEvent(AliESDEvent* const esdEvent){fESDEvent = esdEvent;}

  //Get the Chi2 of particle
  Float_t Chi2() const {return fChi2;}

  ///Get Invariant mass of particle
  Float_t IMass() const { return fIMass; }

  ///Get track or MC labels
  Int_t GetLabel1() const{return fLabel[0];}
  Int_t GetLabel2() const {return fLabel[1];}
  Int_t GetTrackLabel(Int_t i) const {return fLabel[i];}
  void  GetTrackLabels(Int_t * labels) { labels[0] = fLabel[0]; labels[1] = fLabel[1];} 

  Int_t GetMCLabel(Int_t Label) const;

  /* This function returns the Gamma MC label */
  Int_t GetGammaMCLabel() const;
  
  /* This function returns the unique id  of the electrons (if they have the same mother and unique id) */
  Int_t GetElectronUniqueID() const;
  Int_t GetElectronUniqueID(Int_t mcLabel1, Int_t mcLabel2) const;
  
  /// Get MC labels of electrons
  Int_t GetElectronMCLabel1() const;
  Int_t GetElectronMCLabel2() const;

  Bool_t IsMySpawn(const Int_t trackId, const Int_t nSpawn, const Int_t * const spawn) const;

 private:

  Int_t fLabel[2];
  Float_t fChi2; // Chi sq of reconstructed mother
  Float_t fIMass; //Invariant mass, 0 for photons
  Bool_t fTagged; // Is it tagged as decay pion (only for gammas)
  AliStack* fMCStack; //!transient pointer to the mc stack
  AliESDEvent * fESDEvent; //!transient pointer to the esdevent

  ClassDef(AliAODConversionParticle,2)
};


#endif



