/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please aknowledge the author of this code.
*/
#ifndef AliGFW__H
#define AliGFW__H
#include "AliGFWCumulant.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <complex>
using std::vector;
using std::complex;
using std::string;
class AliGFW {
 public:
  struct Region {
    Int_t Nhar, Npar, NpT;
    vector<Int_t> NparVec;
    Double_t EtaMin=-999;
    Double_t EtaMax=-999;
    Int_t BitMask=1;
    string rName="";
    bool operator<(const Region& a) const {
      return EtaMin < a.EtaMin;
    };
    void PrintStructure() {printf("%s: eta [%f.. %f].",rName.c_str(),EtaMin,EtaMax); };
  };
  struct CorrConfig {
    vector<vector<Int_t>> Regs {};
    vector<vector<Int_t>> Hars {};
    vector<Int_t> Overlap;
    vector<Int_t> ptInd;
    Bool_t pTDif=kFALSE;
    string Head="";
  };
  AliGFW();
  ~AliGFW();
  vector<Region> fRegions;
  vector<AliGFWCumulant> fCumulants;
  void AddRegion(string refName, Int_t lNhar, Int_t lNpar, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT=1, Int_t BitMask=1);
  void AddRegion(string refName, Int_t lNhar, Int_t *lNparVec, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT=1, Int_t BitMask=1);
  Int_t CreateRegions();
  void Fill(Double_t eta, Int_t ptin, Double_t phi, Double_t weight, Int_t mask, Double_t secondWeight=-1);
  void Clear();
  AliGFWCumulant GetCumulant(Int_t index) { return fCumulants.at(index); };
  CorrConfig GetCorrelatorConfig(string config, string head = "", Bool_t ptdif=kFALSE);
  complex<Double_t> Calculate(CorrConfig corconf, Int_t ptbin, Bool_t SetHarmsToZero, Bool_t DisableOverlap=kFALSE);
public:
  Bool_t fInitialized;
  complex<Double_t> TwoRec(Int_t n1, Int_t n2, Int_t p1, Int_t p2, Int_t ptbin, AliGFWCumulant*, AliGFWCumulant*, AliGFWCumulant*);
  complex<Double_t> RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> &hars, vector<Int_t> &pows); //POI, Ref. flow, overlapping region
  complex<Double_t> RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> &hars); //POI, Ref. flow, overlapping region
  void AddRegion(Region inreg) { fRegions.push_back(inreg); };
  Region GetRegion(Int_t index) { return fRegions.at(index); };
  Int_t FindRegionByName(string refName);
  //Calculating functions:
  complex<Double_t> Calculate(Int_t poi, Int_t ref, vector<Int_t> hars, Int_t ptbin=0); //For differential, need POI and reference
  complex<Double_t> Calculate(Int_t poi, vector<Int_t> hars); //For integrated case
  //Operations on strings. Equivalent to TString operations, but one to rid of root dependence
  Int_t s_index(string &instr, const string &pattern, const Int_t &spos=0);
  Bool_t s_contains(string &instr, const string &pattern);
  void s_replace(string &instr, const string &pattern1, const string &pattern2, const Int_t &spos=0);
  void s_replace_all(string &instr, const string &pattern1, const string &pattern2);
  Bool_t s_tokenize(string &instr, string &substr, Int_t &spos, const string &delim);
};
#endif
