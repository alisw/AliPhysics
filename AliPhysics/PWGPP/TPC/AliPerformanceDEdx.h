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
class TRootIOCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceDEdx : public AliPerformanceObject {
public :
  AliPerformanceDEdx(TRootIOCtor*);
  AliPerformanceDEdx(const Char_t* name="AliPerformanceDEdx", const Char_t* title="AliPerformanceDEdx",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE, Bool_t useSparse=kTRUE);
  
  virtual ~AliPerformanceDEdx();

  // Init data members
  virtual void Init();

  // Execute analysis
    virtual void Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

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
  void  ProcessTPC(AliMCEvent* const mcev, AliVTrack *const vTrack); // not implemented
  void  ProcessInnerTPC(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent *const vEvent);
  void  ProcessTPCITS(AliMCEvent* const mcev, AliVTrack *const vTrack);      // not implemented
  void  ProcessConstrained(AliMCEvent* const mcev, AliVTrack *const vTrack); // not implemented

  
  // produce summary (currently not used)
  virtual TTree* CreateSummary();

  void FilldEdxHisotgram(double *vDeDxHisto);
    
  //
  // TPC dE/dx 
  //
  THnSparse* GetDeDxHisto() const {return fDeDxHisto;}
  TObjArray* GetHistos() const { return fFolderObj; }
  TCollection* GetListOfDrawableObjects();
    
  virtual void ResetOutputData();

private:

  static Bool_t fgMergeTHnSparse;
  static Bool_t fgUseMergeTHnSparse;
  
  // TPC dE/dx 
  THnSparseF *fDeDxHisto; //-> signal:phi:y:z:snp:tgl:ncls:p:nclsDEdx:nclsF
  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  
  TObjArray* fFolderObj; // array of analysed histograms
  TH1D *h_tpc_dedx_mips_0; //!
  TH1D *h_tpc_dedx_mipsele_0; //!
  TH2D *h_tpc_dedx_mips_c_0_5; //!
  TH2D *h_tpc_dedx_mips_a_0_5; //!
  TH2D *h_tpc_dedx_mips_c_0_1; //!
  TH2D *h_tpc_dedx_mips_a_0_1; //!
  

  AliPerformanceDEdx(const AliPerformanceDEdx&); // not implemented
  AliPerformanceDEdx& operator=(const AliPerformanceDEdx&); // not implemented

  ClassDef(AliPerformanceDEdx,8);
};

#endif
