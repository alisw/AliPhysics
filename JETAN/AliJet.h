#ifndef ALIJET_H
#define ALIJET_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Jet class 
// Stores the output of a jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
 
#include <TObject.h>
#include <TArrayI.h>


class TClonesArray;
class TLorentzVector;
 
class AliJet : public TObject
{
 public:
 
  AliJet();
  ~AliJet();

  // Getters
  Int_t GetNinput() const { return fNInput; }
  Int_t GetNJets() const {return fNJets;}
  TClonesArray* GetJets() const {return fJets;}
  TArrayI GetInJet() const {return fInJet;}
  TArrayI GetMultiplicities() const {return fMultiplicities;}

  TLorentzVector* GetJet(Int_t i);
  Int_t GetMultiplicity(Int_t i);
  Double_t GetPx(Int_t i);
  Double_t GetPy(Int_t i);
  Double_t GetPz(Int_t i);
  Double_t GetP(Int_t i);
  Double_t GetE(Int_t i);
  Double_t GetPt(Int_t i);
  Double_t GetEta(Int_t i);
  Double_t GetPhi(Int_t i);
  Double_t GetTheta(Int_t i);
  Double_t GetMass(Int_t i);
   
  // Setters
  void SetNinput(Int_t i) {fNInput = i;}
  void AddJet(Double_t px, Double_t py, Double_t pz, Double_t e);
  void SetInJet(Int_t* j);
  void SetMultiplicities(Int_t* m);

  // others
  Bool_t OutOfRange(Int_t i, const char *s) const;
  void ClearJets(Option_t *option="");
  void PrintJets();

 protected:

  Int_t fNInput;               // number of input objects
  Int_t fNJets;                // number of jets found
  TArrayI fInJet;              // i-input object belongs to k-jet 
  TArrayI fMultiplicities;     // Multiplicity of each jet
  TClonesArray* fJets;         // 4-momenta of jets

  ClassDef(AliJet,1)
};
 
#endif
