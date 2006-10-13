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

  AliEMCALJet (const AliEMCALJet&);
  AliEMCALJet & operator = (const AliEMCALJet & ) {
    Fatal("operator =", "not implemented") ;
    return *this ;
  }

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
  Float_t Energy() const  {return fEnergy;}
  Float_t EMCALEnergy() const  {return fEMCALEnergy;}
  Float_t EMCALEnergyBGSub()  const {return fEMCALEnergyBGSub;}
  Float_t TrackEnergy() const  {return fTrackEnergy;}
  Float_t TrackEnergyPtCut() const  {return fTrackEnergyPtCut;}
  Float_t HCEnergy() const  {return fHCEnergy;}
  Float_t Phi() const     {return fPhi;}
  Float_t Eta()  const    {return fEta;}
  Int_t   TrackList(Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg)const;
  Int_t   NTracks() const {return fNt;} 
  
protected:
  Float_t  fEnergy;      // Jet Energy
  Float_t  fEMCALEnergy; // EMCAL component of Energy inside Jet cone before BG subtraction
  Float_t  fEMCALEnergyBGSub; // EMCAL component of Energy inside Jet cone after BG subtraction
  Float_t  fTrackEnergy; // Charge tracks component of Energy inside Jet cone with no pT cut
  Float_t  fTrackEnergyPtCut; // Charge tracks component of Energy inside Jet cone after pT cut
  Float_t  fHCEnergy;    // HC  component of Energy inside Jet cone
  Bool_t   fIsWeightedEnergy; // Store flag regarding energy calculation
  Float_t  fEta;      // Jet Eta
  Float_t  fPhi;      // Jet Phi
  Int_t    fNt;       // Number of associated tracks
  Float_t  *fPtT;      //[fNt] Track pt 
  Float_t  *fEtaT;     //[fNt] Track eta
  Float_t  *fPhiT;     //[fNt] Track phi
  Int_t    *fPdgT;     //[fNt] Track pdg code
  ClassDef(AliEMCALJet,8) // Jet for EMCAL

} ;

#endif // ALIEMCALJet_H
