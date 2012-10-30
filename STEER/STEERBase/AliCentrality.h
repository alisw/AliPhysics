//-*- Mode: C++ -*-
#ifndef ALICentrality_H
#define ALICentrality_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentrality
//   author: Alberica Toia
//*****************************************************

#include "TNamed.h"

class AliCentrality : public TNamed
{
 public:

  AliCentrality();  /// constructor
  ~AliCentrality();  /// destructor
  AliCentrality(const AliCentrality& cnt); /// copy constructor
  AliCentrality& operator=(const AliCentrality& cnt);   /// assignment operator

  /// set centrality result
  void SetQuality(Int_t quality) {fQuality = quality;} 
  void SetCentralityV0M(Float_t cent) {fCentralityV0M = cent;} 
  void SetCentralityV0A(Float_t cent) {fCentralityV0A = cent;} 
  void SetCentralityV0C(Float_t cent) {fCentralityV0C = cent;} 
  void SetCentralityFMD(Float_t cent) {fCentralityFMD = cent;}
  void SetCentralityTRK(Float_t cent) {fCentralityTRK = cent;}
  void SetCentralityTKL(Float_t cent) {fCentralityTKL = cent;}
  void SetCentralityCL0(Float_t cent) {fCentralityCL0 = cent;}
  void SetCentralityCL1(Float_t cent) {fCentralityCL1 = cent;}
  void SetCentralityCND(Float_t cent) {fCentralityCND = cent;}
  void SetCentralityNPA(Float_t cent) {fCentralityNPA = cent;}
  void SetCentralityZNA(Float_t cent) {fCentralityZNA = cent;}
  void SetCentralityV0MvsFMD(Float_t cent) {fCentralityV0MvsFMD = cent;}
  void SetCentralityTKLvsV0M(Float_t cent) {fCentralityTKLvsV0M = cent;}
  void SetCentralityZEMvsZDC(Float_t cent) {fCentralityZEMvsZDC = cent;}

  /// get centrality result
  Float_t GetCentralityPercentile(const char *method) const;
  Int_t   GetCentralityClass10(const char *method) const;
  Int_t   GetCentralityClass5(const char *method) const;
  Bool_t  IsEventInCentralityClass(Float_t a, Float_t b, const char *method) const;

  Float_t GetCentralityPercentileUnchecked(const char *method) const;
  Int_t   GetCentralityClass10Unchecked(const char *method) const;
  Int_t   GetCentralityClass5Unchecked(const char *method) const;
  Bool_t  IsEventInCentralityClassUnchecked(Float_t a, Float_t b, const char *method) const;

  Int_t GetQuality() const;
  void  Reset();

 private:
  Int_t   fQuality; // Quality of centrality determination
  Float_t fCentralityV0M;   // Centrality from V0A+V0C
  Float_t fCentralityV0A;   // Centrality from V0A
  Float_t fCentralityV0C;   // Centrality from V0C
  Float_t fCentralityFMD;   // Centrality from FMD
  Float_t fCentralityTRK;   // Centrality from tracks
  Float_t fCentralityTKL;   // Centrality from tracklets
  Float_t fCentralityCL0;   // Centrality from Clusters in layer 0
  Float_t fCentralityCL1;   // Centrality from Clusters in layer 1
  Float_t fCentralityCND;   // Centrality from tracks (candle condition)
  Float_t fCentralityNPA;   // Centrality from Npart (MC)
  Float_t fCentralityZNA;   // Centrality from ZNA)
  Float_t fCentralityV0MvsFMD;   // Centrality from V0 vs FMD
  Float_t fCentralityTKLvsV0M;   // Centrality from tracklets vs V0
  Float_t fCentralityZEMvsZDC;   // Centrality from ZEM vs ZDC

  ClassDef(AliCentrality, 5)
};
#endif //ALICENTRALITY_H
