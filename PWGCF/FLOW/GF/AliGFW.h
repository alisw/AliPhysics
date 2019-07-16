#ifndef AliGFW__H
#define AliGFW__H
#include "AliGFWCumulant.h"
#include <vector>
#include <utility>
#include <algorithm>
#include "TString.h"
#include "TObjArray.h"
using std::vector;

class AliGFW {
 public:
  struct Region {
    Int_t Nhar, Npar, NpT;
    vector<Int_t> NparVec;
    Double_t EtaMin=-999;
    Double_t EtaMax=-999;
    Int_t BitMask=1;
    TString rName="";
    bool operator<(const Region& a) const {
      return EtaMin < a.EtaMin;
    };
    Region operator=(const Region& a) {
      Nhar=a.Nhar;
      Npar=a.Npar;
      NparVec=a.NparVec;
      NpT =a.NpT;
      EtaMin=a.EtaMin;
      EtaMax=a.EtaMax;
      rName=a.rName;
      BitMask=a.BitMask;
      return *this;
    };
    void PrintStructure() {printf("%s: eta [%f.. %f].",rName.Data(),EtaMin,EtaMax); };
  };
  struct CorrConfig {
    vector<Int_t> Regs {};
    vector<Int_t> Hars {};
    vector<Int_t> Regs2 {};
    vector<Int_t> Hars2 {};
    Bool_t pTDif=kFALSE;
    TString Head="";
  };
  AliGFW();
  ~AliGFW();
  vector<Region> fRegions;
  vector<AliGFWCumulant> fCumulants;
  vector<Int_t> fEmptyInt;
  void AddRegion(TString refName, Int_t lNhar, Int_t lNpar, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT=1, Int_t BitMask=1);
  void AddRegion(TString refName, Int_t lNhar, Int_t *lNparVec, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT=1, Int_t BitMask=1);
  Int_t CreateRegions();
  void Fill(Double_t eta, Int_t ptin, Double_t phi, Double_t weight, Int_t mask);
  void Clear();// { for(auto ptr = fCumulants.begin(); ptr!=fCumulants.end(); ++ptr) ptr->ResetQs(); };
  AliGFWCumulant GetCumulant(Int_t index) { return fCumulants.at(index); };
  TComplex Calculate(TString config, Bool_t SetHarmsToZero=kFALSE);
  CorrConfig GetCorrelatorConfig(TString config, TString head = "", Bool_t ptdif=kFALSE);
  TComplex Calculate(CorrConfig corconf, Int_t ptbin, Bool_t SetHarmsToZero, Bool_t DisableOverlap=kFALSE);
 private:
  Bool_t fInitialized;
  void SplitRegions();
  AliGFWCumulant fEmptyCumulant;
  TComplex TwoRec(Int_t n1, Int_t n2, Int_t p1, Int_t p2, Int_t ptbin, AliGFWCumulant*, AliGFWCumulant*, AliGFWCumulant*);
  TComplex RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> hars, vector<Int_t> pows={}); //POI, Ref. flow, overlapping region
  //Deprecated and not used (for now):
  void AddRegion(Region inreg) { fRegions.push_back(inreg); };
  Region GetRegion(Int_t index) { return fRegions.at(index); };
  Int_t FindRegionByName(TString refName);
  vector<TString> fCalculatedNames;
  vector<TComplex> fCalculatedQs;
  Int_t FindCalculated(TString identifier);
  //Calculateing functions:
  TComplex Calculate(Int_t poi, Int_t ref, vector<Int_t> hars, Int_t ptbin=0); //For differential, need POI and reference
  TComplex Calculate(Int_t poi, vector<Int_t> hars); //For integrated case
  //Process one string (= one region)
  TComplex CalculateSingle(TString config);

  Bool_t SetHarmonicsToZero(TString &instr);

};
#endif
