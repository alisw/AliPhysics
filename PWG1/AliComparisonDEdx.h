#ifndef ALICOMPARISONDEdx_H
#define ALICOMPARISONDEdx_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC dE/dx).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliRecInfoCuts;
class AliMCInfoCuts;
class TString;
class TNamed;
class TCanvas;

#include "THnSparse.h"
#include "AliComparisonObject.h"

class AliComparisonDEdx : public AliComparisonObject {
public :
  AliComparisonDEdx(); 
  AliComparisonDEdx(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliComparisonDEdx();

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
  TFolder *CreateFolder(TString folder = "folderDEdx",TString title = "Analysed DEdx histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Process events
  void  ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC);
  void  ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC); // not implemented
  void  ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC); // not implemented

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0)   {fCutsMC = cuts;} 

  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}      
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}     

  static TH1F*     MakeResol(TH2F * his, Int_t integ, Bool_t type); 

  //
  // TPC dE/dx 
  //
  THnSparse* GetDeDxHisto() const {return fDeDxHisto;}

private:

  // TPC dE/dx 
  THnSparseF *fDeDxHisto; //-> signal:alpha:y:z:snp:tgl:ncls:pid:p
  
  // Selection cuts
  AliRecInfoCuts*  fCutsRC; // selection cuts for reconstructed tracks
  AliMCInfoCuts*   fCutsMC; // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonDEdx(const AliComparisonDEdx&); // not implemented
  AliComparisonDEdx& operator=(const AliComparisonDEdx&); // not implemented

  ClassDef(AliComparisonDEdx,1);
};

#endif
