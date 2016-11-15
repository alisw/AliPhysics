#ifndef ALIPERFORMANCEMATCH_H
#define ALIPERFORMANCEMATCH_H

//------------------------------------------------------------------------------
// Class keeps matching information between 
// central barrel detectors.   
// 
// Author: J.Otwinowski   17/10/2009  
// Changes by M.Knichel   22/10/2010
// Changes by J.Salzwedel 14/10/2014
//------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;

class AliVTrack;
class AliMCEvent;
class AliStack;
class AliTrackReference;
class AliVEvent; 
class AliVfriendEvent; 
class AliVfriendTrack; 
class AliMCEvent;
class AliMCParticle;
class AliExternalTrackParam;
class TRootIoCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceMatch : public AliPerformanceObject {
public :
  AliPerformanceMatch(TRootIoCtor*);
  AliPerformanceMatch(const Char_t* name="AliPerformanceMatch", const Char_t* title="AliPerformanceMatch",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE, Bool_t useSparse=kTRUE);
  virtual ~AliPerformanceMatch();

  // Init data members
  virtual void  Init();

  // Execute analysis
  virtual void  Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent,AliVfriendEvent *const vfriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Project Histograms store in AnalysisFolder
  virtual void Analyse();
  
  // Analyse Projected Histograms to create Efficiency and AddToFolder
  //  virtual void AnalyseFinal();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Process matching
  void ProcessTPCITS(AliStack* const stack, AliVEvent *const vEvent, AliVTrack *const vTrack);
  void ProcessTPCTRD(AliStack* const stack, AliVTrack *const vTrack, AliVfriendTrack *const friendTrack);
  void ProcessITSTPC(Int_t trackIdx, AliVEvent* const vEvent, AliStack* const stack, AliVTrack *const vTrack);
  void ProcessTPCConstrain(AliStack* const stack, AliVEvent *const vEvent, AliVTrack *const vTrack); // - 01.11.2011

  // Fill histogrrams
  void FillHistograms(AliVTrack *const refParam, AliVTrack *const param, Bool_t isRec);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderRes",TString title = "Analysed Resolution histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  TH1F*  MakeResol(TH2F * his, Int_t integ=0, Bool_t type=kFALSE, Int_t cut=0); 

  // getters
  //
  THnSparse *GetResolHisto() const  { return fResolHisto; }
  THnSparse *GetPullHisto()  const  { return fPullHisto; }
  THnSparse *GetTrackEffHisto() const  { return fTrackingEffHisto; }
  THnSparse *GetTPCConstrain() const { return fTPCConstrain; }

  TObjArray* GetHistos() const { return fFolderObj; }
  
  static Bool_t GetMergeTHnSparse() { return fgMergeTHnSparse; }
  static void SetMergeTHnSparse(Bool_t mergeTHnSparse) {fgUseMergeTHnSparse = kTRUE; fgMergeTHnSparse = mergeTHnSparse; }
  
  void SetUseHLT(Bool_t useHLT = kTRUE) {fUseHLT = useHLT;}
  Bool_t GetUseHLT() { return fUseHLT; }  

  TObjArray* GetListOfDrawableObjects() {TObjArray* tmp = fFolderObj; fFolderObj = NULL; return tmp;}
    
  virtual void ResetOutputData();

private:

  static Bool_t fgMergeTHnSparse;
  static Bool_t fgUseMergeTHnSparse;

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
  THnSparseF *fTrackingEffHisto;  //-> has match:y:z:snp:tgl:phi:pt:ITSclusters

  // TPC Inner constrained to global tracks - 01.11.2011
  THnSparseF *fTPCConstrain;  //-> pull_phi:phi,pt,eta

  
  TObjArray* fFolderObj; // array of analysed histograms  

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  
  Bool_t fUseHLT; // use HLT ESD

  AliPerformanceMatch(const AliPerformanceMatch&); // not implemented
  AliPerformanceMatch& operator=(const AliPerformanceMatch&); // not implemented

  ClassDef(AliPerformanceMatch,6);
};

#endif
