#ifndef ALICOMPARISONOBJECT_H
#define ALICOMPARISONOBJECT_H

//------------------------------------------------------------------------------
// Abstract class to keep information from comparison of 
// reconstructed and MC particle tracks.   
// 
// Author: J.Otwinowski 04/14/2008 
//------------------------------------------------------------------------------

#include "TNamed.h"
#include "TFolder.h"

class AliMCInfo;
class AliESDRecInfo;

class AliComparisonObject : public TNamed {
public :
  AliComparisonObject(); 
  AliComparisonObject(const char* name="AliComparisonObject", const char* title="AliComparisonObject"); 
  virtual ~AliComparisonObject();

  // Init data members
  // call once before event loop
  virtual void Init() = 0;

  // Execute analysis
  // call in the event loop 
  virtual void Exec(AliMCInfo* const infoMC=0, AliESDRecInfo* const infoRC=0) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list=0) = 0;

  // Analyse output histograms
  virtual void Analyse() = 0;

  // Get output folder for analysed histograms
  virtual TFolder* GetAnalysisFolder() const = 0;

  // set and get analysisMode
  void SetAnalysisMode(Int_t analysisMode=0) {fAnalysisMode = analysisMode;} 
  Int_t GetAnalysisMode() {return fAnalysisMode;}

  // set and get hpt generator 
  void SetHptGenerator(Bool_t hptGenerator=kFALSE) {fHptGenerator = hptGenerator;}
  Bool_t IsHptGenerator() {return fHptGenerator;}

protected: 

 // analysis mode
 Int_t fAnalysisMode;  // 0-TPC, 1-TPCITS, 2-Constrained

 // hpt generator
 Bool_t fHptGenerator; // hpt event generator

  ClassDef(AliComparisonObject,1);
};

#endif
