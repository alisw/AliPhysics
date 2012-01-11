#ifndef ALIPERFORMANCERES_H
#define ALIPERFORMANCERES_H

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
class AliTrackReference;
class AliESDEvent; 
class AliESDfriend; 
class AliESDfriendTrack; 
class AliMCEvent;
class AliMCParticle;
class AliMCInfoCuts;
class AliRecInfoCuts;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceRes : public AliPerformanceObject {
public :
  AliPerformanceRes(); 
  AliPerformanceRes(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  virtual ~AliPerformanceRes();

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
  void ProcessConstrained(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent );
  void ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent);
  void ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent);
  void ProcessInnerTPC(AliMCEvent *const mcEvent, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent);
  void ProcessOuterTPC(AliMCEvent *const mcEvent, AliESDtrack *const esdTrack, AliESDfriendTrack *const friendTrack, AliESDEvent* const esdEvent);

  AliTrackReference *GetFirstTPCTrackRef(AliMCParticle *mcParticle); 
  AliTrackReference *GetLastTPCTrackRef(AliMCParticle *mcParticle); 

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderRes",TString title = "Analysed Resolution histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}   
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}  
   
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}  

  TH1F*  MakeResol(TH2F * his, Int_t integ=0, Bool_t type=kFALSE, Int_t cut=0); 

  // getters
  //
  THnSparse *GetResolHisto() const  { return fResolHisto; }
  THnSparse *GetPullHisto()  const  { return fPullHisto; }
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

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceRes(const AliPerformanceRes&); // not implemented
  AliPerformanceRes& operator=(const AliPerformanceRes&); // not implemented

  ClassDef(AliPerformanceRes,1);
};

#endif
