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
  void RebinMulti(Int_t nbins);
  void RebinMulti(Int_t nbins, Double_t *binedges);
  TH1 *getHist(Int_t ind=-1);
  ClassDef(AliProfileBS,1);
protected:
  TH1* getHistRebinned(TProfile *inpf); //Performs rebinning, if required, and returns a projection of profile
  Bool_t fProfInitialized;
  Int_t fNSubs;
  Int_t fMultiRebin; //! externaly set runtime, no need to store
  Double_t *fMultiRebinEdges; //! externaly set runtime, no need to store
};
#endif
