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
  AliPHOSIhepAnalyze(const AliPHOSIhepAnalyze & obj) : TObject(obj), fRunLoader(0), fFileName()
  {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
  virtual ~AliPHOSIhepAnalyze() {}  ; // dtor
  AliPHOSIhepAnalyze & operator = (const AliPHOSIhepAnalyze & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

  void AnalyzeCPV1(Int_t Nevents); // resolutions, mult and cluster lengths for CPV
  void AnalyzeEMC1(Int_t Nevents); // resolutions, mult and cluster lengths for EMC
  void AnalyzeCPV2(Int_t Nevents); // delta(gen)/delta(rec) between hits 
  void CpvSingle(Int_t Nevents); // signle particle analysis
  virtual void HitsCPV(Int_t event); 
  const TString GetFileName() const { return fFileName; }

 private:

  Bool_t IsCharged(Int_t pdgCode);

 private:
 
  AliRunLoader *fRunLoader; // run loader
  TString fFileName;        // filename with headers (e.g. galice.root)

ClassDef(AliPHOSIhepAnalyze,1)  // PHOSv1 event analyzis algorithm

};

#endif // AliPHOSIhepAnalyze_H


