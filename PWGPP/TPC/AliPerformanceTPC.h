#ifndef ALIPERFORMANCETPC_H
#define ALIPERFORMANCETPC_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC resolution).   
// 
// Author: J.Otwinowski 04/02/2008 
// Changes by M.Knichel 15/10/2010
// Changes by J.Salzwedel 01/10/2014
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1;
class TH2;
class TH3;
class TH1D;
class TH2D;
class TNtuple;

class AliVTrack;
class AliMCEvent;
class AliStack;
class AliVEvent;
class AliVEvent;
class AliVfriendEvent; 
class TRootIoCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceTPC : public AliPerformanceObject {
public :
  AliPerformanceTPC(TRootIoCtor*);
  AliPerformanceTPC(const Char_t* name="AliPerformanceTPC", const Char_t* title="AliPerformanceTPC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE, Int_t run=-1, Bool_t highMult = kFALSE, Bool_t useSparse = kTRUE);

  AliPerformanceTPC(const AliPerformanceTPC&);
  AliPerformanceTPC& operator=(const AliPerformanceTPC&);

  virtual ~AliPerformanceTPC();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const infoMC, AliVEvent* const infoRC, AliVfriendEvent* const vfriendEvent, const Bool_t bUseMC=kFALSE, const Bool_t bUseVfriend=kFALSE);
  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list=0);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}
  
  // produce summary
  virtual TTree* CreateSummary();
    
  // Process events
  void ProcessConstrained(AliStack* const stack, AliVTrack *const vTrack, AliVEvent *const vEvent);
  void ProcessTPC(AliStack* const stack, AliVTrack *const vTrack, AliVEvent *const vEvent, Bool_t vertStatus);
  void ProcessTPCITS(AliStack* const stack, AliVTrack *const vTrack, AliVEvent *const vEvent, Bool_t vertStatus);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderTPC",TString title = "Analysed TPC performance histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

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
  TCollection* GetListOfDrawableObjects() {TObjArray* tmp = fFolderObj; fFolderObj = NULL; return tmp;}

  virtual void ResetOutputData();

    
private:

  static Bool_t fgMergeTHnSparse;
  static Bool_t fgUseMergeTHnSparse;  

  // TPC histogram
  THnSparseF *fTPCClustHisto; //-> padRow:phi:TPCside
  THnSparseF *fTPCEventHisto;  //-> Xv:Yv:Zv:mult:multP:multN:vertStatus
  THnSparseF *fTPCTrackHisto;  //-> nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge:vertStatus
  TObjArray* fFolderObj; // array of analysed histograms

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  Bool_t fUseHLT; // use HLT ESD
  ClassDef(AliPerformanceTPC,14);
};

#endif
