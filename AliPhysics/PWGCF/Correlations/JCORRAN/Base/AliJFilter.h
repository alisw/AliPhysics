// $Id: AliJFilter.h,v 1.5 2012/04/19 15:19:52 jkral Exp $

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIJFILTER_H
#define ALIJFILTER_H

#include "TNamed.h"
#include "AliJRunHeader.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSE.h"
#include <iostream>

#include <AliJConst.h>
#include <TVectorT.h>

//==============================================================

#ifndef AliJMaxDimBuffer
#define AliJMaxDimBuffer
const int kMaxDimBuffer = 300;//max length of a line read to a buffe
#endif

class AliJEventHeader;
class AliJRunHeader;
class AliJTrack;
class AliAnalysisTaskSE;

class TH1D;
class TH2D;
class TNtuple;
class TList;
class TTree;
class TFormula;
class TRefArray;
class TArrayI;

class AliMCEvent; 
class AliAODEvent; 
class AliAODTrack; 
class AliVCluster;
class AliVEvent;

class AliMCEvent;
class AliAnalysisFilter;

class AliPIDResponse;
class AliPIDResponse;
class AliPIDCombined;
class AliAnalysisUtils;
class AliJRunTable;

using namespace std;

class AliJFilter : public TNamed  {

 public:
  AliJFilter();
  AliJFilter(const char *name,  AliAnalysisTaskSE *task);
  AliJFilter(const AliJFilter& ap);   
  AliJFilter& operator = (const AliJFilter& ap);
  virtual ~AliJFilter();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");


  bool GetEventSuccess() const { return fEventSuccess; }
  void SetTrackThreshold(double t) { fTrackThreshold = t; }
  void SetMyTask( AliAnalysisTaskSE *t ) { fMyTask = t; }
  TClonesArray *GetTrackList() const { return fTrackList; }
  TClonesArray *GetMCTrackList() const { return fMCTrackList; }
  TClonesArray *GetHeaderList() const { return fHeaderList; }
  TList        *GetRunInfoList() const { return fRunInfoList; }

  TClonesArray **GetTrackListP()  { return &fTrackList; }
  TClonesArray **GetMCTrackListP()  { return &fMCTrackList; }
  TClonesArray **GetHeaderListP() { return &fHeaderList; }
  TList **GetRunInfoListP()  { return &fRunInfoList; }

  Bool_t      GetStoreEventPlaneSource(){ return fAliJRunHeader->GetStoreEventPlaneSource(); }
  AliAODEvent * AODEvent(){ return dynamic_cast<AliAODEvent*>(Event());}
  AliVEvent   * Event(){ return fMyTask->InputEvent(); }
  AliMCEvent  * MCEvent(){ return IsMC()?fMyTask->MCEvent():NULL; }

  Bool_t       IsMC(){ return fAliJRunHeader->IsMC(); }
  Bool_t       FromAOD(){ return fAliJRunHeader->FromAOD(); }
  Bool_t IsGoodEvent(AliAODEvent *event);
 
  AliJRunHeader* GetAliJRunHeader() const { return fAliJRunHeader; }
  void    SetAliJRunHeader( AliJRunHeader* header ){ fAliJRunHeader=header; } 
 private:

  Int_t        DebugLevel(){ return fMyTask->DebugLevel(); }
  inline void   DEBUG(int level, int type, TString msg1, TString msg2=""){
    if(DebugLevel()>level) std::cout<<type<<"\t"<<msg1<<" : "<<msg2<<std::endl;
  }

  AliJEventHeader* ReadCommonHeader(AliAODEvent *event);
  AliJEventHeader* ReadCommonHeader(AliVEvent *event);
  // methods to read data from AOD
  Bool_t ReadAODTracks(const AliAODEvent* aod);
  void ReadAODHeader(AliAODEvent* aod);
  void ReadFilter();
  void ReadMCTracksFromAOD();
	void RemapMCLabels();

  UInt_t ConvertTriggerMask();//Converts alice trigger mask to JCorran trigger mask
  //functions used for event selction:
  bool AcceptAODTrack(AliAODTrack* aodTrack);
  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char* GetOADBPath() const { return fOADBPath.Data(); }

  // method to fill jcorran
  void PrintOut() const;
  
  // UTILS
  void AddList(const char* aname, const char* cname, TClonesArray **obj, int nlist);

  // d a t a     m e m b e r s
  TVectorT<double>  fIsRealOrMC; // flags if the input are real (0) ESDs or MonteCarlo ESDs (1)
  TString fActiveTriggers[kRangeTriggerTableAlice]; // alice table mapping trigger bit to trigger name
  TString fTriggerTableJCorran[kRangeTriggerTableJCorran]; // JCorran trigger table TBit 0 =MinBias
  Bool_t fStoreEventPlaneSource; // store event plane
  TString fOADBPath; // oadb path
  Double_t fTrackThreshold; // for event tropping
  Bool_t fEventSuccess; //! if filter was successful with current event

  TArrayI *fMcMap; //! mc index map

  // jcorran output objects
  TClonesArray *    fTrackList;   //! list of charged track objects
  TClonesArray *    fMCTrackList; //! list of charged track objects
  TClonesArray *    fHeaderList;  //! event details
  TList *        fRunInfoList; //! run details
  AliPIDResponse  *fPIDResponse; //! PID response object
  AliPIDCombined  *fPIDCombined; //! PID response object

  AliJRunHeader*      fAliJRunHeader; //!  run details (mg field, trigger mask,etc...)
  AliAnalysisUtils  *fAnaUtils; //! analysis utils ALICE
  AliAnalysisTaskSE *fMyTask; //! task pointer

  TF1 *pfOutlierLowCut, *pfOutlierHighCut; //! event selection 
  Bool_t fFirstEvent; //!
  AliJRunTable *fRunTable; //!


  ClassDef(AliJFilter, 1); 
};
#endif // AliJFilter_H
