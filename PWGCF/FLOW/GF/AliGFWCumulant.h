/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
A part of <AliGFW.cxx/h>
A container to store Q vectors for one subevent with an extra layer to recursively calculate particle correlations.
If used, modified, or distributed, please aknowledge the author of this code.
*/
#ifndef ALIGFWCUMULANT__H
#define ALIGFWCUMULANT__H
#include "RtypesCore.h" //needed for root types only. Irrelevant, if running c++-SA code
#include <cmath>
#include <complex>
#include <vector>
using std::vector;
using std::complex;
class AliGFWCumulant {
 public:
  AliGFWCumulant();
  ~AliGFWCumulant();
  void ResetQs();
  void FillArray(Int_t ptin, Double_t phi, Double_t weight=1, Double_t SecondWeight=-1);
  enum UsedFlags_t {kBlank = 0, kFull=1, kPt=2};
  void SetType(UInt_t infl) { DestroyComplexVectorArray(); fUsed = infl; };
  void Inc() { fNEntries++; };
  Int_t GetN() { return fNEntries; };
  Bool_t IsPtBinFilled(Int_t ptb);
  void CreateComplexVectorArray(Int_t N=1, Int_t P=1, Int_t Pt=1);
  void CreateComplexVectorArrayVarPower(Int_t N=1, vector<Int_t> Pvec={1}, Int_t Pt=1);
  Int_t PW(Int_t ind) { return fPowVec.at(ind); }; //No checks to speed up, be carefull!!!
  void DestroyComplexVectorArray();
  complex<Double_t> Vec(Int_t, Int_t, Int_t ptbin=0); //envelope class to summarize pt-dif. Q-vec getter
 protected:
  complex<Double_t> ***fQvector;
  UInt_t fUsed;
  Int_t fNEntries;
  //Q-vectors. Could be done recursively, but maybe defining each one of them explicitly is easier to read
  Int_t fN; //! Harmonics
  Int_t fPow; //! Power
  vector<Int_t> fPowVec; //! Powers array
  Int_t fPt; //!fPt bins
  Bool_t *fFilledPts;
  Bool_t fInitialized; //Arrays are initialized
  complex<Double_t> fNullQ = complex<Double_t>(0.,0.);
};

#endif
