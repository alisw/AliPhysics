#ifndef ALIPROFILESUBSET__H
#define ALIPROFILESUBSET__H
//Helper function to select a subrange of a TProfile
#include "TProfile.h"
#include "TArrayD.h"
class AliProfileSubset : public TProfile {
public:
  AliProfileSubset();
  AliProfileSubset(TProfile *inpf);
  ~AliProfileSubset();
  Int_t GetNCells() { return fNcells; };
  TArrayD Get_fBinEntries() { return fBinEntries; };
  TArrayD Get_fBinSumw2() { return fBinSumw2; };
  TArrayD Get_fSumw2() { return fSumw2; };
  Double_t* Get_fArray() { return fArray; };
  EErrorType GetErrorMode() { return fErrorMode; };
  Double_t GetSumwy() { return fTsumwy; };
  Double_t GetSumwy2() { return fTsumwy2; };
  Double_t GetYMax() { return fYmax; };
  Double_t GetYMin() { return fYmin; };
  Bool_t CopyFromProfile(Int_t StartI, Int_t StopI, TProfile *source);
  ClassDef(AliProfileSubset,1);
};
#endif
