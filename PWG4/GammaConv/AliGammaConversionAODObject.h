#ifndef ALIGAMMACONVERSIONAODOBJECT_H
#define ALIGAMMACONVERSIONAODOBJECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class containing the aod information from conversions
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include "TObject.h" 
class AliStack;
class AliESDEvent;
#include "TMath.h"
class AliGammaConversionAODObject : public TObject {

 public: 

  AliGammaConversionAODObject();                                        //constructor
  AliGammaConversionAODObject(const AliGammaConversionAODObject & g);                   //copy constructor
  AliGammaConversionAODObject & operator = (const AliGammaConversionAODObject & g);     //assignment operator
  virtual ~AliGammaConversionAODObject() {;}                            //virtual destructor


  /* This function sets the Px */
  void SetPx(Float_t px){fPx = px;}

  /* This function sets the Py */
  void SetPy(Float_t py){fPy = py;}
  
  /* This function sets the Pz */
  void SetPz(Float_t pz){fPz = pz;}

  ///Set E
  void SetE(Float_t e) { fE = e; }

  ///Set the Chi2 of reconstructed conversion gamma
  void SetChi2(Float_t chi2) {fChi2 = chi2;}

  /* This function sets the esd label of the first decay particle */
  void SetLabel1(Int_t label){fLabel1 = label;}

  /* This function sets the esd label of the second decay particle */
  void SetLabel2(Int_t label){fLabel2 = label;}
  
  ///Set Invariant mass
  void SetIMass(Float_t im) { fIMass = im; }

  ///Set the tag for decay meson
  void SetTag( Bool_t tagged ) { fTagged = tagged; }

  /* This function sets the MC stack */
  void SetStack(AliStack* const stack){fMCStack=stack;}
  
  /* This function sets the ESD event*/
  void SetESDEvent(AliESDEvent* const esdEvent){fESDEvent = esdEvent;}

  /// Get Px
  Float_t Px() const{return fPx;}

  /* This function returns the Py*/
  Float_t Py() const{return fPy;}

  /* This function returns the Pz*/
  Float_t Pz() const{return fPz;}

  //Get the Chi2 of particle
  Float_t Chi2() const {return fChi2;}

  //Get Pt of particle
  Float_t Pt() const {return TMath::Sqrt(fPx* fPx + fPy*fPy);}
  
  ///Get Energy of particle
  Float_t E() const {return fE;}

  ///Get Phi of particle
  Float_t Phi() const {    return fPx == 0.0 && fPy == 0.0 ? 0.0 : TMath::ATan2(fPy,fPx); }

  ///Get Invariant mass of particle
  Float_t IMass() const { return fIMass; }

  /* This function returns the esd label of the first electron */
  Int_t GetLabel1() const{return fLabel1;}
  
  /* This function returns the esd label of the second electron */
  Int_t GetLabel2()const {return fLabel2;}

  /// Get gamma tag (if tagged as decay photon)
  Bool_t IsTagged() const { return fTagged; }

  /* This function returns the Gamma MC label */
  Int_t GetGammaMCLabel() const;
  
  /* This function returns the unique id  of the electrons (if they have the same mother and unique id) */
  Int_t GetElectronUniqueID() const;
  
  /* This function returns the MC label of the first electron*/
  Int_t GetElectronMCLabel1() const;
  
  /* This function returns the MC label of the second electron*/
  Int_t GetElectronMCLabel2() const;

 private:

  Float_t fPx; // Px 
  Float_t fPy; // Py 
  Float_t fPz; // Pz 
  Float_t fE; // Pz 
  Int_t fLabel1; //mc label or TClonesArray index of first decay particle
  Int_t fLabel2; //mc label or TClonesArray index of second decay particle
  Float_t fChi2; // Chi sq of reconstructed mother
  Float_t fIMass; //Invariant mass, 0 for photons
  Bool_t fTagged; // Is it tagged as decay pion (only for gammas)
  AliStack* fMCStack; //!transient pointer to the mc stack
  AliESDEvent * fESDEvent; //!transient pointer to the esdevent

  ClassDef(AliGammaConversionAODObject,2)
};


#endif



