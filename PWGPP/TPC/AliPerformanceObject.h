#ifndef ALIPERFORMANCEOBJECT_H
#define ALIPERFORMANCEOBJECT_H

//------------------------------------------------------------------------------
// Base class to keep information from comparison of 
// reconstructed and MC particle tracks.   
// 
// Author: J.Otwinowski 04/14/2008 
// Changes by M.Knichel 15/10/2010
// Changes by J.Salzwedel 29/9/2014
//------------------------------------------------------------------------------

#include "TNamed.h"
#include "TFolder.h"
#include "THnSparse.h"
#include "AliMergeable.h"

class TTree;
class AliMCEvent;
class AliVEvent;
class AliRecInfoCuts;
class AliMCInfoCuts;
class AliVfriendEvent;
class AliESDVertex;
class TRootIOCtor;
#include "AliRecInfoCuts.h"
#include "AliMCInfoCuts.h"

class AliPerformanceObject : public TNamed, public AliMergeable {
public :
  AliPerformanceObject(TRootIOCtor*); 
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
  virtual void Exec(AliMCEvent* const infoMC=0, AliVEvent* const infoRC=0, AliVfriendEvent* const vfriendEvent=0, const Bool_t bUseMC=kFALSE, const Bool_t bUseVfriend=kFALSE) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list=0) = 0;

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
  virtual TCollection* GetListOfDrawableObjects(){ return 0; }

  // Selection cuts
  void SetAliRecInfoCuts(const AliRecInfoCuts* cuts) {
    if (!cuts) return;
    fCutsRC = *cuts;
  }
  void SetAliMCInfoCuts(const AliMCInfoCuts* cuts) {
    fCutsMC = *cuts;
  }

  // set and get analysisMode
  void SetAnalysisMode(const Int_t analysisMode=0) {fAnalysisMode = analysisMode;} 
  Int_t GetAnalysisMode() const {return fAnalysisMode;}

  // set and get hpt generator 
  void SetHptGenerator(const Bool_t hptGenerator=kFALSE) {fHptGenerator = hptGenerator;}
  Bool_t IsHptGenerator() const {return fHptGenerator;}

  // draw all histograms from the folder
  void PrintHisto(Bool_t logz = kTRUE, const Char_t * outFileName = "PerformanceQAHisto.ps"); 
    
  // create log axis 
  Double_t *CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax); 

  // trigger class selection
  void SetTriggerClass(const Char_t *triggerClass) { fTriggerClass = triggerClass; }
  const Char_t* GetTriggerClass() const { return fTriggerClass.IsNull()?NULL:fTriggerClass.Data(); }

  // use track vertex
  void SetUseTrackVertex(Bool_t trackVtx = kTRUE) { fUseTrackVertex = trackVtx; }
  Bool_t IsUseTrackVertex() { return fUseTrackVertex; }
  
  Bool_t IsHighMultiplicity() { return fHighMultiplicity; }  
  
  // merging of thnsparse
  Bool_t GetMergeTHnSparseObj() { return fMergeTHnSparseObj; }
  void SetMergeTHnSparseObj(Bool_t merge) {fMergeTHnSparseObj = merge; }  
  
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

  virtual void ResetOutputData() { ; }
    
protected: 

  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, TString* selString = 0);
  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, Int_t yDim, TString* selString = 0);
  void AddProjection(TObjArray* aFolderObj, TString nameSparse, THnSparse *hSparse, Int_t xDim, Int_t yDim, Int_t zDim, TString* selString = 0);

  // merge THnSparse
  Bool_t fMergeTHnSparseObj;
  
  // analysis mode
  Int_t fAnalysisMode;  // 0-TPC, 1-TPCITS, 2-Constrained, 3-TPC inner wall, 4-TPC outer wall

  Int_t fRunNumber;

  // hpt generator
  Bool_t fHptGenerator; // hpt event generator

  // trigger class
  TString fTriggerClass;

  // use track vertex
  Bool_t fUseTrackVertex; // use track vertex
  
  // PbPb mode?
  Bool_t fHighMultiplicity; // flag to switch between pp and PbPb  

  Bool_t fUseKinkDaughters; // use kink daughthers, default is yes

  Int_t  fUseCentralityBin;  // centrality bin to be used 

  Bool_t fUseTOFBunchCrossing; // use TOFBunchCrossing, default is yes
  Bool_t fUseSparse;

  // Global cuts objects
  AliRecInfoCuts fCutsRC;  // selection cuts for reconstructed tracks
  AliMCInfoCuts  fCutsMC;  // selection cuts for MC tracks

  ClassDef(AliPerformanceObject,11);
};

#endif
