#ifndef ALICOMPARISONEFF_H
#define ALICOMPARISONEFF_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC efficiency).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TFile;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliESD;
class AliRecInfoCuts;
class AliMCInfoCuts;
class TString;
class AliESDVertex;
class TNamed;

#include "THnSparse.h"
#include "AliComparisonObject.h"

class AliComparisonEff : public AliComparisonObject {
public :
  AliComparisonEff(); 
  AliComparisonEff(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliComparisonEff();

  // Init data members
  virtual void Init();

  // Execute analysis 
  virtual void Exec(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC);

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
  void ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);

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

  AliComparisonEff(const AliComparisonEff&); // not implemented
  AliComparisonEff& operator=(const AliComparisonEff&); // not implemented

  ClassDef(AliComparisonEff,1);
};

#endif
