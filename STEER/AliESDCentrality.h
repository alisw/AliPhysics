//-*- Mode: C++ -*-
#ifndef ALIESDCentrality_H
#define ALIESDCentrality_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************

#include "TNamed.h"

class AliESDCentrality : public TNamed
{
 public:

  AliESDCentrality();  /// constructor
  ~AliESDCentrality();  /// destructor
  AliESDCentrality(const AliESDCentrality& cnt); /// copy constructor
  AliESDCentrality& operator=(const AliESDCentrality& cnt);   /// assignment operator

  /// set centrality result
  void SetCentralityV0M(Float_t cent) {fCentralityV0M = cent;} 
  void SetCentralityFMD(Float_t cent) {fCentralityFMD = cent;}
  void SetCentralityTRK(Float_t cent) {fCentralityTRK = cent;}
  void SetCentralityTKL(Float_t cent) {fCentralityTKL = cent;}
  void SetCentralityCL0(Float_t cent) {fCentralityCL0 = cent;}
  void SetCentralityCL1(Float_t cent) {fCentralityCL1 = cent;}
  void SetCentralityV0MvsFMD(Float_t cent) {fCentralityV0MvsFMD = cent;}
  void SetCentralityTKLvsV0M(Float_t cent) {fCentralityTKLvsV0M = cent;}
  void SetCentralityZEMvsZDC(Float_t cent) {fCentralityZEMvsZDC = cent;}

  /// get centrality result
  Float_t GetCentralityPercentile(const char *method);
  Int_t   GetCentralityClass10(const char *method);
  Int_t   GetCentralityClass5(const char *method);
  Bool_t  IsEventInCentralityClass(Float_t a, Float_t b, const char *method);

 private:
  Float_t fCentralityV0M;   // Centrality from V0
  Float_t fCentralityFMD;   // Centrality from FMD
  Float_t fCentralityTRK;   // Centrality from tracks
  Float_t fCentralityTKL;   // Centrality from tracklets
  Float_t fCentralityCL0;   // Centrality from Clusters in layer 0
  Float_t fCentralityCL1;   // Centrality from Clusters in layer 0
  Float_t fCentralityV0MvsFMD;   // Centrality from V0 vs FMD
  Float_t fCentralityTKLvsV0M;   // Centrality from tracklets vs V0
  Float_t fCentralityZEMvsZDC;   // Centrality from ZEM vs ZDC

  ClassDef(AliESDCentrality, 1)
};
#endif //ALIESDCENTRALITY_H
