#ifndef ALIPERFORMANCEMATCH_H
#define ALIPERFORMANCEMATCH_H

//------------------------------------------------------------------------------
// Class keeps matching information between 
// central barrel detectors.   
// 
// Author: J.Otwinowski 17/10/2009 
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
class AliExternalTrackParam;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceMatch : public AliPerformanceObject {
public :
  AliPerformanceMatch(); 
  AliPerformanceMatch(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  virtual ~AliPerformanceMatch();

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

  // Process matching
  void ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack, AliESDfriendTrack *const friendTrack);
  void ProcessTPCTRD(AliStack* const stack, AliESDtrack *const esdTrack, AliESDfriendTrack *const friendTrack);
  void ProcessITSTPC(Int_t trackIdx, AliESDEvent* const esdEvent, AliStack* const stack, AliESDtrack *const esdTrack, AliESDfriendTrack *const friendTrack);

  // Fill histogrrams
  void FillHistograms(AliExternalTrackParam *const refParam, AliExternalTrackParam *const param, Bool_t isRec);

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
  THnSparse *GetTrackEffHisto() const  { return fTrackingEffHisto; }

private:
  //
  // Control histograms
  // 5 track parameters (details in STEER/AliExternalTrackParam.h)
  // + isRec flag to determine ITS/TRD tracking efficiency
  // w.r.t TPC

  // resolution histogram
  THnSparseF *fResolHisto; //-> res_y:res_z:res_phi:res_lambda:res_pt:y:z:phi:eta:pt:isRec

  // pull histogram
  THnSparseF *fPullHisto;  //-> pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt:isRec

  // tracking efficiency using ITS stand-alone tracks histogram
  THnSparseF *fTrackingEffHisto;  // -> has match:y:z:snp:tgl:pt:ITSclusters

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceMatch(const AliPerformanceMatch&); // not implemented
  AliPerformanceMatch& operator=(const AliPerformanceMatch&); // not implemented

  ClassDef(AliPerformanceMatch,1);
};

#endif
