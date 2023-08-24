/*
Author: Vytautas Vislavicius
Container to store multiple profiles used for bootstrapping.
Extra layer to calculate the moments. Current implementation supports only the first and second moments.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#ifndef ALIPROFILEBS__H
#define ALIPROFILEBS__H
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TCollection.h"
#include "TMath.h"


class AliProfileBS: public TProfile {
public:
  AliProfileBS();
  ~AliProfileBS();
  AliProfileBS(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins);
  AliProfileBS(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup);
  TList *fListOfEntries;
  void MergeBS(AliProfileBS *target);
  void InitializeSubsamples(Int_t nSub);
  void FillProfile(const Double_t &xv, const Double_t &yv, const Double_t &w, const Double_t &rn);
  void FillProfile(const Double_t &xv, const Double_t &yv, const Double_t &w);
  Long64_t Merge(TCollection *collist);
  void RebinMulti(Int_t nbins);
  void RebinMulti(Int_t nbins, Double_t *binedges);
  TH1 *getHist(Int_t ind=-1);
  TProfile *getProfile(Int_t ind=-1);
  TProfile *getSummedProfiles();
  void OverrideMainWithSub();
  Int_t getNSubs() { return fListOfEntries->GetEntries(); };
  void PresetWeights(AliProfileBS *targetBS) { fPresetWeights = targetBS; };
  void ResetBin(Int_t nbin) { ResetBin((TProfile*)this,nbin); for(Int_t i=0;i<fListOfEntries->GetEntries(); i++) ResetBin((TProfile*)fListOfEntries->At(i),nbin); };
  ClassDef(AliProfileBS,2);
protected:
  TH1* getHistRebinned(TProfile *inpf); //Performs rebinning, if required, and returns a projection of profile
  TH1* getWeightBasedRebin(Int_t ind=-1);
  Bool_t fProfInitialized;
  Int_t fNSubs;
  Int_t fMultiRebin; //! externaly set runtime, no need to store
  Double_t *fMultiRebinEdges; //! externaly set runtime, no need to store
  AliProfileBS *fPresetWeights; //! AliProfileBS whose weights we should copy
  void ResetBin(TProfile *tpf, Int_t nbin) {tpf->SetBinEntries(nbin,0); tpf->SetBinContent(nbin,0); tpf->SetBinError(nbin,0); };
};
#endif
