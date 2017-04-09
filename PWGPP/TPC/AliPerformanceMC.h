#ifndef ALIPERFORMANCEMC_H
#define ALIPERFORMANCEMC_H

//------------------------------------------------------------------------------
// Class to keep information for MC  tracks 
// to check propagation algorithm, B-field and material budget.   
// 
// Author: J.Otwinowski 09/06/2009 
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;
class TParticle;
class TClonesArray;

class AliESDVertex;
class AliESDtrack;
class AliMCEvent;
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

class AliPerformanceMC : public AliPerformanceObject {
public :
  AliPerformanceMC(const Char_t* name="AliPerformanceMC", const Char_t* title="AliPerformanceMC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE);
  virtual ~AliPerformanceMC();

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
  void ProcessTPC(AliTrackReference* const refIn, TParticle* const particle);
  void ProcessInnerTPC(AliTrackReference* const refIn,  AliTrackReference* const refOut, TParticle* const particle);
  void ProcessOuterTPCExt(TParticle *const  part, TClonesArray * const trefs);
  AliExternalTrackParam * MakeTrack(const AliTrackReference* ref, TParticle *const part);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderMC", TString title = "Analysed MC histograms");

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
  THnSparseF *fPullHisto;  //-> pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceMC(const AliPerformanceMC&); // not implemented
  AliPerformanceMC& operator=(const AliPerformanceMC&); // not implemented

  ClassDef(AliPerformanceMC,1);
};

#endif
