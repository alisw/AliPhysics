#ifndef ALIEMCALJET_H
#define ALIEMCALJET_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
#include <TObject.h>
//*-- Author: Andreas Morsch (CERN)


class AliEMCALJet : public TObject {
 public:
  AliEMCALJet();
  AliEMCALJet(Float_t energy, Float_t phi, Float_t eta); 
  virtual ~AliEMCALJet();
  void SetEnergy(Float_t val) {fEnergy = val;}
  void SetEMCALEnergy(Float_t val) {fEMCALEnergy = val;}
  void SetEMCALEnergyBGSub(Float_t val){fEMCALEnergyBGSub = val;}
  void SetTrackEnergy(Float_t val) {fTrackEnergy = val;}
  void SetTrackEnergyPtCut(Float_t val){fTrackEnergyPtCut = val;}
  void SetHCEnergy(Float_t val) {fHCEnergy = val;}
  void SetPhi(Float_t val)    {fPhi    = val;}  
  void SetEta(Float_t val)    {fEta    = val;}    
  void SetIsWeightedEnergy(Bool_t flag)    {fIsWeightedEnergy    = flag;}    
  void SetTrackList(Int_t val, Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg);
  Float_t Energy()  {return fEnergy;}
  Float_t EMCALEnergy()  {return fEMCALEnergy;}
  Float_t EMCALEnergyBGSub()  {return fEMCALEnergyBGSub;}
  Float_t TrackEnergy()  {return fTrackEnergy;}
  Float_t TrackEnergyPtCut()  {return fTrackEnergyPtCut;}
  Float_t HCEnergy()  {return fHCEnergy;}
  Float_t Phi()     {return fPhi;}
  Float_t Eta()     {return fEta;}
  Int_t   TrackList(Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg);
  Int_t   NTracks() {return fNt;} 
  
protected:
  Float_t  fEnergy;      // Jet Energy
  Float_t  fEMCALEnergy; // EMCAL component of Energy inside Jet cone before BG subtraction
  Float_t  fEMCALEnergyBGSub; // EMCAL component of Energy inside Jet cone after BG subtraction
  Float_t  fTrackEnergy; // Charge tracks component of Energy inside Jet cone with no pT cut
  Float_t  fTrackEnergyPtCut; // Charge tracks component of Energy inside Jet cone after pT cut
  Float_t  fHCEnergy;    // HC  component of Energy inside Jet cone
  Bool_t   fIsWeightedEnergy; // Store flag regarding energy calculation
  Float_t  fEta;      // Jet Phi
  Float_t  fPhi;      // Jet Eta
  Int_t    fNt;       // Number of associated tracks
  Float_t  fPtT [1000]; // Track pt 
  Float_t  fEtaT[1000]; // Track eta
  Float_t  fPhiT[1000]; // Track phi
  Int_t    fPdgT[1000]; // Track pdg code
  ClassDef(AliEMCALJet,6) // Jet for EMCAL

} ;

#endif // ALIEMCALJet_H
