#ifndef ALIPERFORMANCETPC_H
#define ALIPERFORMANCETPC_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC resolution).   
// 
// Author: J.Otwinowski 04/02/2008 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1;
class TH2;
class TH3;

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
  //AliPerformanceTPC(); 
  AliPerformanceTPC(Char_t* name="AliPerformanceTPC", Char_t* title="AliPerformanceTPC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE, Int_t run=-1, Bool_t highMult = kFALSE);

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
  
  // produce summary
  virtual TTree* CreateSummary();

  // Process events
  void ProcessConstrained(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent);
  void ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent, Bool_t vertStatus);
  void ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent *const esdEvent, Bool_t vertStatus);

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
  
  TObjArray* GetHistos() const { return fFolderObj; }
  
  static Bool_t GetMergeTHnSparse() { return fgMergeTHnSparse; }
  static void SetMergeTHnSparse(Bool_t mergeTHnSparse) {fgUseMergeTHnSparse = kTRUE; fgMergeTHnSparse = mergeTHnSparse; }
  
  void SetUseHLT(Bool_t useHLT = kTRUE) {fUseHLT = useHLT;}
  Bool_t GetUseHLT() { return fUseHLT; }


private:

  static Bool_t fgMergeTHnSparse;
  static Bool_t fgUseMergeTHnSparse;  

  // TPC histogram
  THnSparseF *fTPCClustHisto; //-> padRow:phi:TPCside
  THnSparseF *fTPCEventHisto;  //-> Xv:Yv:Zv:mult:multP:multN:vertStatus
  THnSparseF *fTPCTrackHisto;  //-> nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge:vertStatus
  TObjArray* fFolderObj; // array of analysed histograms

  // Global cuts objects
  AliRecInfoCuts* fCutsRC;  // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;  // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  Bool_t fUseHLT; // use HLT ESD

  AliPerformanceTPC(const AliPerformanceTPC&); // not implemented
  AliPerformanceTPC& operator=(const AliPerformanceTPC&); // not implemented

  ClassDef(AliPerformanceTPC,11);
};

#endif
