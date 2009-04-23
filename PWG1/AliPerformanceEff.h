#ifndef ALIPERFORMANCEEFF_H
#define ALIPERFORMANCEEFF_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC efficiency).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TFile;
class TParticle;
class TString;
class TNamed;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliMCEvent; 
class AliESDEvent; 
class AliMCParticle; 
class AliESDtrack;
class AliESD;
class AliRecInfoCuts;
class AliMCInfoCuts;
class AliESDVertex;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceEff : public AliPerformanceObject {
public :
  AliPerformanceEff(); 
  AliPerformanceEff(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliPerformanceEff();

  // Init data members
  virtual void Init();

  // Execute analysis 
  virtual void  Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, const Bool_t bUseMC);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms 
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderEff",TString title = "Analysed Efficiency histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray *array=0);

  // Process events
  void ProcessTPC(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent);
  void ProcessTPCITS(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent);
  void ProcessConstrained(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent);

  Bool_t IsRecTPC(AliESDtrack *track);
  Bool_t IsRecTPCITS(AliESDtrack *track);
  Bool_t IsRecConstrained(AliESDtrack *track);

  Bool_t IsFindable(AliMCEvent *mcEvent, Int_t label);
  Int_t TransformToPID(TParticle *particle);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;} 
  
  // Getters
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;} 
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

  THnSparseF* GetEffHisto() const {return fEffHisto;}

private:

  // Control histograms
  THnSparseF *fEffHisto; //-> mceta:mcphi:mcpt:pid:isPrim:recStatus:findable

  // Global cuts objects
  AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceEff(const AliPerformanceEff&); // not implemented
  AliPerformanceEff& operator=(const AliPerformanceEff&); // not implemented

  ClassDef(AliPerformanceEff,1);
};

#endif
