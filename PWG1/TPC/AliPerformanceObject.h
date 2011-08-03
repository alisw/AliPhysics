#ifndef ALIPERFORMANCEOBJECT_H
#define ALIPERFORMANCEOBJECT_H

//------------------------------------------------------------------------------
// Base class to keep information from comparison of 
// reconstructed and MC particle tracks.   
// 
// Author: J.Otwinowski 04/14/2008 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

#include "TNamed.h"
#include "TFolder.h"
#include "THnSparse.h"

class TTree;
class AliMCEvent;
class AliESDEvent;
class AliRecInfoCuts;
class AliMCInfoCuts;
class AliESDfriend;
class AliESDVertex;

class AliPerformanceObject : public TNamed {
public :
  AliPerformanceObject(); 
  AliPerformanceObject(const char* name="AliPerformanceObject", const char* title="AliPerformanceObject", Int_t run=-1, Bool_t highMult=kFALSE); 
  virtual ~AliPerformanceObject();

  // Init data members
  // call once before event loop
  virtual void Init() = 0;
  
  // init for high multiplicity (PbPb) 
  // to be called instead of Init()
  virtual void InitHighMult();
  
  // Execute analysis
  // call in the event loop 
  virtual void Exec(AliMCEvent* const infoMC=0, AliESDEvent* const infoRC=0, AliESDfriend* const infoFriend=0, const Bool_t bUseMC=kFALSE, const Bool_t bUseESDfriend=kFALSE) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list=0) = 0;

  // project to 1d,2d,3d
  // is called from FinishTaskOuput() in AliPerformanceTask
  virtual void Analyse() = 0;

  // Get output folder for analysed histograms
  virtual TFolder* GetAnalysisFolder() const = 0;
  
  // create a summary stored in a ttree 
  // has to be implented
  virtual TTree* CreateSummary() { return 0; }
  
  // project to 1d,2d,3d
  // is called from Terminate() in AliPerformanceTask
  // final spectra calculation
  virtual void AnalyseFinal() { ; }

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
  
  Bool_t IsHighMultiplicity() { return fHighMultiplicity; }  
  
  void SetRunNumber(Int_t run) { fRunNumber = run; }
  Int_t GetRunNumber() const { return fRunNumber; }

  // use kink daughters
  void SetUseKinkDaughters(Bool_t kinkDaughters = kTRUE) { fUseKinkDaughters = kinkDaughters; }
  Bool_t IsUseKinkDaughters() { return fUseKinkDaughters; }

  // Centrality bin to be used
  void  SetUseCentralityBin(Int_t bin) { fUseCentralityBin = bin; }
  Int_t GetUseCentralityBin()          { return fUseCentralityBin; }

  // use tof bunch crossing
  void SetUseTOFBunchCrossing(Bool_t tofBunching = kTRUE) { fUseTOFBunchCrossing = tofBunching; }
  Bool_t IsUseTOFBunchCrossing() { return fUseTOFBunchCrossing; }

protected: 

  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, TString* selString = 0);
  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, Int_t yDim, TString* selString = 0);
  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, Int_t yDim, Int_t zDim, TString* selString = 0);

  // analysis mode
  Int_t fAnalysisMode;  // 0-TPC, 1-TPCITS, 2-Constrained, 3-TPC inner wall, 4-TPC outer wall

  Int_t fRunNumber;

  // hpt generator
  Bool_t fHptGenerator; // hpt event generator

  // trigger class
  const Char_t * fTriggerClass;

  // use track vertex
  Bool_t fUseTrackVertex; // use track vertex
  
  // PbPb mode?
  Bool_t fHighMultiplicity; // flag to switch between pp and PbPb  

  Bool_t fUseKinkDaughters; // use kink daughthers, default is yes

  Int_t  fUseCentralityBin;  // centrality bin to be used 

  Bool_t fUseTOFBunchCrossing; // use TOFBunchCrossing, default is yes

  AliPerformanceObject(const AliPerformanceObject&); // not implemented
  AliPerformanceObject& operator=(const AliPerformanceObject&); // not implemented

  ClassDef(AliPerformanceObject,6);
};

#endif
