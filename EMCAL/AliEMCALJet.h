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
  void SetPhi(Float_t val)    {fPhi    = val;}  
  void SetEta(Float_t val)    {fEta    = val;}    
  void SetTrackList(Int_t val, Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg);
  Float_t Energy()  {return fEnergy;}
  Float_t Phi()     {return fPhi;}
  Float_t Eta()     {return fEta;}
  Int_t   TrackList(Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg);
  Int_t   NTracks() {return fNt;} 
  
protected:
  Float_t  fEnergy;   // Jet Energy
  Float_t  fEta;      // Jet Phi
  Float_t  fPhi;      // Jet Eta
  Int_t    fNt;       // Number of associated tracks
  Float_t  fPtT [1000]; // Track pt 
  Float_t  fEtaT[1000]; // Track eta
  Float_t  fPhiT[1000]; // Track phi
  Int_t    fPdgT[1000]; // Track pdg code
  ClassDef(AliEMCALJet,4) // Jet for EMCAL

} ;

#endif // ALIEMCALJet_H
