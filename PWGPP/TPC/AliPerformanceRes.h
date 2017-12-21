#ifndef ALIPERFORMANCERES_H
#define ALIPERFORMANCERES_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC resolution).   
// 
// Author: J.Otwinowski 04/02/2008
// Changes by J.Salzwedel 23/10/2014
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;

class AliVTrack;
class AliMCEvent;
class AliTrackReference;
class AliVEvent; 
class AliVfriendEvent; 
class AliVfriendTrack; 
class AliMCEvent;
class AliMCParticle;
class TRootIOCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceRes : public AliPerformanceObject {
public :
  AliPerformanceRes(TRootIOCtor*);
  AliPerformanceRes(const Char_t* name="AliPerformanceRes", const Char_t* title="AliPerformanceRes",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE);
  virtual ~AliPerformanceRes();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const infoMC=0, AliVEvent* const infoRC=0, AliVfriendEvent* const vfriendEvent=0, const Bool_t bUseMC=kFALSE, const Bool_t bUseVfriend=kFALSE);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Process events
  void ProcessConstrained(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent );
  void ProcessTPC(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent);
  void ProcessTPCITS(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent);
  void ProcessInnerTPC(AliMCEvent *const mcEvent, AliVTrack *const vTrack, AliVEvent* const vEvent);
  void ProcessOuterTPC(AliMCEvent *const mcEvent, AliVTrack *const vTrack, const AliVfriendTrack *const friendTrack, AliVEvent* const vEvent);

  AliTrackReference *GetFirstTPCTrackRef(AliMCParticle *mcParticle); 
  AliTrackReference *GetLastTPCTrackRef(AliMCParticle *mcParticle); 

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderRes",TString title = "Analysed Resolution histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  TH1F*  MakeResol(TH2F * his, Int_t integ=0, Bool_t type=kFALSE, Int_t cut=0); 

  // getters
  //
  THnSparse *GetResolHisto() const  { return fResolHisto; }
  THnSparse *GetPullHisto()  const  { return fPullHisto; }
  static void SetMergeEntriesCut(Double_t entriesCut){fgkMergeEntriesCut = entriesCut;}

private:
  //
  // Control histograms
  // 5 track parameters (details in STEER/AliExternalTrackParam.h)
  //

  // resolution histogram
  THnSparseF *fResolHisto; //-> res_y:res_z:res_phi:res_lambda:res_pt:y:z:phi:eta:pt

  // pull histogram
  //THnSparseF *fPullHisto;  //-> pull_y:pull_z:pull_phi:pull_lambda:pull_1pt:y:z:eta:phi:pt
  THnSparseF *fPullHisto;  //-> pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceRes(const AliPerformanceRes&); // not implemented
  AliPerformanceRes& operator=(const AliPerformanceRes&); // not implemented
  static Double_t            fgkMergeEntriesCut;  //maximal number of entries for merging  -can be modified via setter

  ClassDef(AliPerformanceRes,3);
};

#endif
