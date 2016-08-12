#ifndef ALIPERFORMANCEDEdx_H
#define ALIPERFORMANCEDEdx_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC dE/dx).   
// 
// Author: J.Otwinowski 04/02/2008 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

class TCanvas;
class TH1F;
class TH2F;
class TNamed;
class TString;

class AliVEvent;
class AliVfriendEvent;
class AliMCEvent;
class AliVTrack;
class AliStack;
class AliRecInfoCuts;
class AliMCInfoCuts;


#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceDEdx : public AliPerformanceObject {
public :
  AliPerformanceDEdx(const Char_t* name="AliPerformanceDEdx", const Char_t* title="AliPerformanceDEdx",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE, Bool_t useSparse=kTRUE);
  
  virtual ~AliPerformanceDEdx();

  // Init data members
  virtual void Init();

  // Execute analysis
    virtual void Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}
  
  static Bool_t GetMergeTHnSparse() { return fgMergeTHnSparse; }
  static void SetMergeTHnSparse(Bool_t mergeTHnSparse) {fgUseMergeTHnSparse = kTRUE; fgMergeTHnSparse = mergeTHnSparse; }

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderDEdx",TString title = "Analysed DEdx histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Process events
    void  ProcessTPC(AliStack* const stack, AliVTrack *const vTrack); // not implemented
    void  ProcessInnerTPC(AliStack* const stack, AliVTrack *const vTrack, AliVEvent *const vEvent);
    void  ProcessTPCITS(AliStack* const stack, AliVTrack *const vTrack);      // not implemented
    void  ProcessConstrained(AliStack* const stack, AliVTrack *const vTrack); // not implemented

  
  // produce summary (currently not used)
  virtual TTree* CreateSummary();

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0)   {fCutsMC = cuts;} 
  void FilldEdxHisotgram(double *vDeDxHisto);
    
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

  //
  // TPC dE/dx 
  //
  THnSparse* GetDeDxHisto() const {return fDeDxHisto;}
  TObjArray* GetHistos() const { return fFolderObj; }
  TObjArray* GetListOfDrawableObjects() {TObjArray* tmp = fFolderObj; fFolderObj = NULL; return tmp;}
    
  virtual void ResetOutputData();

private:

  static Bool_t fgMergeTHnSparse;
  static Bool_t fgUseMergeTHnSparse;
  
  // TPC dE/dx 
  THnSparseF *fDeDxHisto; //-> signal:phi:y:z:snp:tgl:ncls:p:nclsDEdx:nclsF
  TObjArray* fFolderObj; // array of analysed histograms
  
  // Selection cuts
  AliRecInfoCuts*  fCutsRC; // selection cuts for reconstructed tracks
  AliMCInfoCuts*   fCutsMC; // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  

  AliPerformanceDEdx(const AliPerformanceDEdx&); // not implemented
  AliPerformanceDEdx& operator=(const AliPerformanceDEdx&); // not implemented

  ClassDef(AliPerformanceDEdx,6);
};

#endif
