#ifndef AliPHOSIhepAnalyze_H
#define AliPHOSIhepAnalyze_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
// Algorythm class to analyze PHOSv1 events:
// Construct histograms and displays them.
// Used the IHEP CPV/PHOS reconstruction algorithm.
//*--
//*-- Author : Boris Polichtchouk (IHEP)

// --- ROOT system ---
#include "TObject.h"

// --- Standard library ---

// --- AliRoot header files ---
class AliRunLoader;

class AliPHOSIhepAnalyze : public TObject {

 public:

  AliPHOSIhepAnalyze() ;              // ctor
  AliPHOSIhepAnalyze(Text_t * name) ; // ctor

  void AnalyzeCPV1(Int_t Nevents); // resolutions, mult and cluster lengths for CPV
  void AnalyzeEMC1(Int_t Nevents); // resolutions, mult and cluster lengths for EMC
  void AnalyzeCPV2(Int_t Nevents); // delta(gen)/delta(rec) between hits 
  void CpvSingle(Int_t Nevents); // signle particle analysis
  virtual void HitsCPV(Int_t event); 
  TString GetFileName() { return fFileName; }

 private:

  Bool_t IsCharged(Int_t pdg_code);

 private:
 
  AliRunLoader *fRunLoader;
  TString fFileName;

ClassDef(AliPHOSIhepAnalyze,1)  // PHOSv1 event analyzis algorithm

};

#endif // AliPHOSIhepAnalyze_H


