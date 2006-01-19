#ifndef ALILEADING_H
#define ALILEADING_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Class to find and store the leading particle in event and
// store its correlation to associated particles
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TObject.h>
#include <TArrayI.h>

class TLorentzVector;
class AliJetReader;

class AliLeading : public TObject
{
 public:

  AliLeading();
  ~AliLeading();

  // Getters
  Int_t GetNassoc() const {return fNassoc;}
  TLorentzVector* GetLeading() const {return fLeading;}
  TArrayI GetCorr() const {return fCorr;}
  Bool_t LeadingFound() const {return fFound;}
  Double_t GetE();
  Double_t GetPt();
  Double_t GetEta();
  Double_t GetPhi();
   
  // Setters
  // Others
  void FindLeading(AliJetReader* reader);
  void PrintLeading();
  void Reset();

 protected:

  Int_t fNassoc;            // number of associated particles
  TLorentzVector* fLeading; // leading particle
  TArrayI fCorr;            // array to store azimuthal correlation
                            // between leading and assoc particles
  Int_t fnBin;              // number of bins in array
  Double_t fLow;            // value corresponding to lower bound of bin 0
  Bool_t  fFound;           // leading found
  
  ClassDef(AliLeading,1);
};

#endif
