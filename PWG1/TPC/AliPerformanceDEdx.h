#ifndef ALIPERFORMANCEDEdx_H
#define ALIPERFORMANCEDEdx_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC dE/dx).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TCanvas;
class TH1F;
class TH2F;
class TNamed;
class TString;

class AliESDEvent; 
class AliESDfriend; 
class AliMCEvent;
class AliESDtrack;
class AliStack; 
class AliRecInfoCuts;
class AliMCInfoCuts;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceDEdx : public AliPerformanceObject {
public :
  AliPerformanceDEdx(); 
  AliPerformanceDEdx(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliPerformanceDEdx();

  // Init data members
  virtual void Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend);

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
  void  ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack); // not implemented
  void  ProcessInnerTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent);
  void  ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack);      // not implemented
  void  ProcessConstrained(AliStack* const stack, AliESDtrack *const esdTrack); // not implemented

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0)   {fCutsMC = cuts;} 

  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}      
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}     

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

  AliPerformanceDEdx(const AliPerformanceDEdx&); // not implemented
  AliPerformanceDEdx& operator=(const AliPerformanceDEdx&); // not implemented

  ClassDef(AliPerformanceDEdx,1);
};

#endif
