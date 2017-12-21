#ifndef ALIPERFORMANCEMC_H
#define ALIPERFORMANCEMC_H

//------------------------------------------------------------------------------
// Class to keep information for MC  tracks 
// to check propagation algorithm, B-field and material budget.   
// 
// Author: J.Otwinowski 09/06/200
// Changes by J.Salzwedel 5/11/2014
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;
class TParticle;
class TClonesArray;

class AliMCEvent;
class AliTrackReference;
class AliVEvent; 
class AliVfriendEvent; 
class AliMCEvent;
class AliMCParticle;
class AliExternalTrackParam;
class TRootIOCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceMC : public AliPerformanceObject {
public :
  AliPerformanceMC(TRootIOCtor*);
  AliPerformanceMC(const Char_t* name="AliPerformanceMC", const Char_t* title="AliPerformanceMC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE);
  virtual ~AliPerformanceMC();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vfriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

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

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceMC(const AliPerformanceMC&); // not implemented
  AliPerformanceMC& operator=(const AliPerformanceMC&); // not implemented

  ClassDef(AliPerformanceMC,2);
};

#endif
