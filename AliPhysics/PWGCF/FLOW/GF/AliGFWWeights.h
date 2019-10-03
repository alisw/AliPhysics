#ifndef ALIGFWWEIGHTS__H
#define ALIGFWWEIGHTS__H
#include "TObjArray.h"
#include "TNamed.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"

class AliGFWWeights: public TNamed
{
 public:
  AliGFWWeights();
  ~AliGFWWeights();
  void Init(Bool_t AddData=kTRUE, Bool_t AddM=kTRUE);
  void Fill(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype); //htype: 0 for data, 1 for mc rec, 2 for mc gen
  Double_t GetWeight(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype); //htype: 0 for data, 1 for mc rec, 2 for mc gen
  Double_t GetNUA(Double_t phi, Double_t eta, Double_t vz); //This just fetches correction from integrated NUA, should speed up
  Bool_t IsDataFilled() { return fDataFilled; };
  Bool_t IsMCFilled() { return fMCFilled; };
  Double_t FindMax(TH3D *inh, Int_t &ix, Int_t &iy, Int_t &iz);
  void MCToEfficiency();
  TObjArray *GetRecArray() { return fW_mcrec; };
  TObjArray *GetGenArray() { return fW_mcgen; };
  TObjArray *GetDataArray(){ return fW_data; };
  void CreateNUA(Bool_t IntegrateOverCentAndPt=kTRUE);
  void CreateNUE(Bool_t IntegrateOverCentrality=kTRUE);
  TH3D *GetIntegratedEfficiency() { return fEffInt; };
  void SetDataFilled(Bool_t newval) { fDataFilled=newval; };
  void SetMCFilled(Bool_t newval) { fMCFilled = newval; };
  void ReadAndMerge(const char *filelinks);
  void SetPtBins(Int_t Nbins, Double_t *bins);
  Long64_t Merge(TCollection *collist);
  void RebinNUA(Int_t nX=1, Int_t nY=2, Int_t nZ=5);
  void OverwriteNUA();
  private:
  Bool_t fDataFilled;
  Bool_t fMCFilled;
  TObjArray *fW_data;
  TObjArray *fW_mcrec;
  TObjArray *fW_mcgen;
  TH3D *fEffInt; //!
  TH3D *fAccInt; //!
  Int_t fNbinsPt; //! do not store
  Double_t *fbinsPt; //! do not store
  void AddArray(TObjArray *targ, TObjArray *sour);
  const char *GetBinName(Double_t ptv, Double_t v0mv,const char *pf="") {
    Int_t ptind = 0;//GetPtBin(ptv);
    Int_t v0mind = 0;//GetV0MBin(v0mv);
    return Form("Bin%s_weights%i_%i",pf,ptind,v0mind);
  };

  ClassDef(AliGFWWeights,1);
};


#endif
