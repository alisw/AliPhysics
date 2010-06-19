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
class AliESDfriend;
class AliESDVertex;

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
  virtual void Exec(AliMCEvent* const infoMC=0, AliESDEvent* const infoRC=0, AliESDfriend* const infoFriend=0, const Bool_t bUseMC=kFALSE, const Bool_t bUseESDfriend=kFALSE) = 0;

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
  void PrintHisto(Bool_t logz = kTRUE, Char_t * outFileName = "PerformanceQAHisto.ps"); 

  // create log axis 
  Double_t *CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax); 

  // trigger class selection
  void SetTriggerClass(const Char_t *triggerClass) { fTriggerClass = triggerClass; }
  const Char_t* GetTriggerClass() const { return fTriggerClass; }

  // use track vertex
  void SetUseTrackVertex(Bool_t trackVtx = kTRUE) { fUseTrackVertex = trackVtx; }
  Bool_t IsUseTrackVertex() { return fUseTrackVertex; }

protected: 

  // analysis mode
  Int_t fAnalysisMode;  // 0-TPC, 1-TPCITS, 2-Constrained, 3-TPC inner wall, 4-TPC outer wall

  // hpt generator
  Bool_t fHptGenerator; // hpt event generator

  // trigger class
  const Char_t * fTriggerClass;

  // use track vertex
  Bool_t fUseTrackVertex; // use track vertex

  AliPerformanceObject(const AliPerformanceObject&); // not implemented
  AliPerformanceObject& operator=(const AliPerformanceObject&); // not implemented

  ClassDef(AliPerformanceObject,1);
};

#endif
