#ifndef ALIEMCALPARTON_H
#define ALIEMCALPARTON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include <TObject.h>
//*-- Author: Renan Cabrera (Creighton U.)


class AliEMCALParton : public TObject {
public:
  AliEMCALParton();
  AliEMCALParton(Float_t energy, Float_t phi, Float_t eta);
  virtual ~AliEMCALParton();
  void SetEnergy(Float_t val) {fEnergy = val;}
  void SetPhi(Float_t val)    {fPhi    = val;}
  void SetEta(Float_t val)    {fEta    = val;}
  void SetTrackList(Int_t num ,Float_t* ar1,Float_t* ar2,Float_t* ar3,Int_t* ar4);
  void GetTrackList(Float_t* ar1,Float_t* ar2,Float_t* ar3,Int_t* ar4)const;
  
  Int_t GetNTracks() const {return fNTracks;}
  Float_t Energy()  const {return fEnergy;}
  Float_t Phi()    const {return fPhi;}
  Float_t Eta()    const {return fEta;}
  void SetPartonCode(Int_t code){fPartonCode = code;} 
  Int_t GetPartonCode()const{return fPartonCode;}
protected:
  Float_t   fEnergy;   // Jet Energy
  Float_t   fEta;      // Jet Phi
  Float_t   fPhi;      // Jet Eta
  Int_t     fNTracks;      // Number of tracks 
  Float_t  fTrackEnergy[1000];  // Jet Tracks Energy
  Float_t  fTrackEta[1000];     // Jet Tracks Eta
  Float_t  fTrackPhi[1000];     // Jet Tracks Phi
  Int_t    fTrackPDG[1000];     // Jet Tracks PDG code
  Int_t	   fPartonCode;		// Store the type of parton
  ClassDef(AliEMCALParton,7) // Jet for EMCAL
    
} ;

#endif // ALIEMCALParton_H
