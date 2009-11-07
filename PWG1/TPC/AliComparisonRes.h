#ifndef ALICOMPARISONRES_H
#define ALICOMPARISONRES_H

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
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliMCInfoCuts;
class AliRecInfoCuts;

#include "THnSparse.h"
#include "AliComparisonObject.h"

class AliComparisonRes : public AliComparisonObject {
public :
  AliComparisonRes(); 
  AliComparisonRes(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  virtual ~AliComparisonRes();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Process events
  void ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderRes",TString title = "Analysed Resolution histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}   
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}  
   
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}  

  static TH1F*     MakeResol(TH2F * his, Int_t integ, Bool_t type); 

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
  THnSparseF *fResolHisto; //-> res_y:res_z:res_phi:res_lambda:res_1pt:y:z:eta:phi:pt

  // pull histogram
  THnSparseF *fPullHisto;  //-> pull_y:pull_z:pull_phi:pull_lambda:pull_1pt:y:z:eta:phi:pt

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonRes(const AliComparisonRes&); // not implemented
  AliComparisonRes& operator=(const AliComparisonRes&); // not implemented

  ClassDef(AliComparisonRes,1);
};

#endif
