#ifndef ALIPERFORMANCEOBJECT_H
#define ALIPERFORMANCEOBJECT_H

//------------------------------------------------------------------------------
// Base class to keep information from comparison of 
// reconstructed and MC particle tracks.   
// 
// Author: J.Otwinowski 04/14/2008 
//------------------------------------------------------------------------------

#include "TNamed.h"
#include "TFolder.h"

class AliMCEvent;
class AliESDEvent;
class AliRecInfoCuts;
class AliMCInfoCuts;

class AliPerformanceObject : public TNamed {
public :
  AliPerformanceObject(); 
  AliPerformanceObject(const char* name="AliPerformanceObject", const char* title="AliPerformanceObject"); 
  virtual ~AliPerformanceObject();

  // Init data members
  // call once before event loop
  virtual void Init() = 0;

  // Execute analysis
  // call in the event loop 
  virtual void Exec(AliMCEvent* const infoMC=0, AliESDEvent* const infoRC=0, const Bool_t bUseMC=kFALSE) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list=0) = 0;

  // Analyse output histograms
  virtual void Analyse() = 0;

  // Get output folder for analysed histograms
  virtual TFolder* GetAnalysisFolder() const = 0;

  // 
  virtual void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) = 0;
  virtual void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) = 0; 

  // set and get analysisMode
  void SetAnalysisMode(const Int_t analysisMode=0) {fAnalysisMode = analysisMode;} 
  Int_t GetAnalysisMode() const {return fAnalysisMode;}

  // set and get hpt generator 
  void SetHptGenerator(const Bool_t hptGenerator=kFALSE) {fHptGenerator = hptGenerator;}
  Bool_t IsHptGenerator() const {return fHptGenerator;}

  // draw all histograms from the folder
  void DrawHisto(Bool_t logz = kTRUE); 

  // create log axis 
  Double_t *CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax); 

protected: 

 // analysis mode
 Int_t fAnalysisMode;  // 0-TPC, 1-TPCITS, 2-Constrained

 // hpt generator
 Bool_t fHptGenerator; // hpt event generator

  ClassDef(AliPerformanceObject,1);
};

#endif
