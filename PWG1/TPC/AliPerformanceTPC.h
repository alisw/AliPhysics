#ifndef ALIPERFORMANCETPC_H
#define ALIPERFORMANCETPC_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC resolution).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;

class AliESDVertex;
class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliESDEvent; 
class AliESDfriend; 
class AliMCInfoCuts;
class AliRecInfoCuts;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceTPC : public AliPerformanceObject {
public :
  AliPerformanceTPC(); 
  AliPerformanceTPC(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  virtual ~AliPerformanceTPC();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend);
  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Process events
  void ProcessConstrained(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent);
  void ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent);
  void ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderTPC",TString title = "Analysed TPC performance histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}   
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}  
   
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}  

  // getters
  //
  THnSparse *GetTPCClustHisto() const  { return fTPCClustHisto; }
  THnSparse *GetTPCEventHisto() const  { return fTPCEventHisto; }
  THnSparse *GetTPCTrackHisto() const  { return fTPCTrackHisto; }


private:

  // TPC histogram
  THnSparseF *fTPCClustHisto; //-> gclX:gclY:TPCSide
  THnSparseF *fTPCEventHisto; //-> Xv:Yv:Zv:mult:multP:multN:vertStatus
  THnSparseF *fTPCTrackHisto; //-> nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge

  // Global cuts objects
  AliRecInfoCuts* fCutsRC;  // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;  // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms


  AliPerformanceTPC(const AliPerformanceTPC&); // not implemented
  AliPerformanceTPC& operator=(const AliPerformanceTPC&); // not implemented

  ClassDef(AliPerformanceTPC,3);
};

#endif
