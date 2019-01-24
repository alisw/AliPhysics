#ifndef ALIGFWWEIGHTS__H
#define ALIGFWWEIGHTS__H
#include "TObjArray.h"
#include "TNamed.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

class AliGFWWeights: public TNamed
{
 public:
  AliGFWWeights();
  ~AliGFWWeights();
  void Init(Bool_t AddData=kTRUE, Bool_t AddM=kTRUE);
  void Fill(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype); //htype: 0 for data, 1 for mc rec, 2 for mc gen
  Double_t GetWeight(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype); //htype: 0 for data, 1 for mc rec, 2 for mc gen
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
  private:
  Bool_t fDataFilled;
  Bool_t fMCFilled;
  TObjArray *fW_data;
  TObjArray *fW_mcrec;
  TObjArray *fW_mcgen;
  TH3D *fEffInt; //!
  TH3D *fAccInt; //!
  void AddArray(TObjArray *targ, TObjArray *sour);
  /*TODO: implement external setters:*/
  static const Int_t fNbinsPtDefault=48;
  Double_t *fbinsPtDefault;
  /*const Int_t fgV0MNBins=10;//9;
  Double_t fgV0MBins[fgV0MNBins+1] = {0,1,5,10,15,20,30,40,50,70,100};
  TAxis *fg_PtAxis = new TAxis(fNbinsPtDefault,fbinsPtDefault);
  TAxis *fg_V0MAxis = new TAxis(fgV0MNBins, fgV0MBins);*/
  //Int_t GetPtBin(Double_t ptval) { return fg_PtAxis->FindBin(ptval); };
  //Int_t GetV0MBin(Double_t V0MVal) { return fg_V0MAxis->FindBin(V0MVal); };
  /*const char *GetBinName(Double_t ptv, Double_t v0mv,const char *pf="") {
    Int_t ptind = GetPtBin(ptv);
    Int_t v0mind = GetV0MBin(v0mv);
    return Form("Bin%s%i_%i",pf,ptind,v0mind);
    };*/
  
  ClassDef(AliGFWWeights,1);
};


#endif
