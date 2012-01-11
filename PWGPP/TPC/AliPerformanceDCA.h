#ifndef ALIPERFORMANCEDCA_H
#define ALIPERFORMANCEDCA_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (DCA - Distance of Closest Approach 
// to the vertex).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class AliESDEvent; 
class AliESDfriend; 
class AliStack; 
class AliRecInfoCuts;
class AliMCInfoCuts;
class AliESDVertex;
class AliESDtrack;
class TH3;
class TH2;
class TH1;
class TString;
class TNamed;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceDCA : public AliPerformanceObject {
public :
  AliPerformanceDCA(); 
  AliPerformanceDCA(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliPerformanceDCA();

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
  TFolder *CreateFolder(TString folder = "folderDCA",TString title = "Analysed DCA histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  void ProcessConstrained(AliStack* const stack, AliESDtrack *const esdTrack);
  void ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent);
  void ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}  

  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

  // getters
  THnSparse* GetDCAHisto() const {return fDCAHisto;}

  // Make stat histograms
  TH1F* MakeStat1D(TH2 *hist, Int_t delta1, Int_t type);
  TH2F* MakeStat2D(TH3 *hist, Int_t delta0, Int_t delta1, Int_t type);

private:

  // DCA histograms
  THnSparseF *fDCAHisto; //-> dca_r:dca_z:eta:pt:phi 
 
  // Global cuts objects
  AliRecInfoCuts* fCutsRC; // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;  // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceDCA(const AliPerformanceDCA&); // not implemented
  AliPerformanceDCA& operator=(const AliPerformanceDCA&); // not implemented

  ClassDef(AliPerformanceDCA,1);
};

#endif
