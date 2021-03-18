#ifndef ALIPROFILEBS__H
#define ALIPROFILEBS__H
//Container to store bootstrapping profiles
//Vytautas Vislavicius
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TCollection.h"


class AliProfileBS: public TProfile {
public:
  AliProfileBS();
  ~AliProfileBS();
  AliProfileBS(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins);
  AliProfileBS(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  TList *fListOfEntries;
  void InitializeSubsamples(Int_t nSub);
  void FillProfile(const Double_t &xv, const Double_t &yv, const Double_t &w, const Double_t &rn);
  void FillProfile(const Double_t &xv, const Double_t &yv, const Double_t &w);
  Long64_t Merge(TCollection *collist);
  ClassDef(AliProfileBS,1);
protected:
  Bool_t fProfInitialized;
  Int_t fNSubs;
};
#endif
