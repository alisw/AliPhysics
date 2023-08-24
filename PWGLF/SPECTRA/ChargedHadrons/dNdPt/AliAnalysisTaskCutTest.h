#ifndef AliAnalysisTaskCutTest_cxx
#define AliAnalysisTaskCutTest_cxx

class THnSparse;
class AliESDEvent;
class AliESDtrackCuts;
class AlidNdPtAcceptanceCuts;
class AliEventCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCutTest : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCutTest(const char *name = "AliAnalysisTaskCutTest");
  virtual ~AliAnalysisTaskCutTest() {}
  
  void SetTrackCuts(AliESDtrackCuts* const cuts)                { fEsdTrackCuts = cuts; }
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  void SetAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts)    { fAccCuts = cuts; }
  AlidNdPtAcceptanceCuts* GetAcceptanceCuts() const             { return fAccCuts; }  
  void SetEventCuts(AliEventCuts* const cuts)              { fEventCuts = cuts; }
  AliEventCuts* GetEventCuts() const                       { return fEventCuts; }
  void SetUseMCInfo(Bool_t useMC = kTRUE)                       { fUseMCInfo = useMC; }
  Bool_t IsUseMCInfo() const                                    { return fUseMCInfo; }
  void SetCentralityInterval(Double_t cmin, Double_t cmax)      { fCentralityMin = cmin; fCentralityMax = cmax; fUseCentrality = kTRUE; }
  void SetDeadZoneWidth(Double_t width)                         { fDeadZoneWidth = width; }
  void SetMaxZ(Double_t maxZ)                                   { fMaxZ = maxZ; }
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Double_t zv, AlidNdPtHelper::TrackObject trackObj, Int_t multMB, Int_t multRecTrack, Int_t multTrueMC);
  
 private:
  AliESDtrackCuts *fEsdTrackCuts;               // esd track cuts 
  AlidNdPtAcceptanceCuts* fAccCuts;             // acceptance cuts for multiplicity
  AliEventCuts *fEventCuts;                // event cuts
  AliESDEvent *fESD;    //! ESD event
  AliMCEvent  *fMC;      //! MC event
  Bool_t      fUseCentrality;
  Double_t    fCentralityMin;
  Double_t    fCentralityMax;
  Double_t    fDeadZoneWidth; 
  Double_t    fMaxZ;
  THnSparseD  *fTrackHistRefitTPC; 
  THnSparseD  *fTrackHistChi2TPC; 
  THnSparseD  *fTrackHistRowsTPC; 
  THnSparseD  *fTrackHistFindableTPC; 
  THnSparseD  *fTrackHistSharedTPC; 
  THnSparseD  *fTrackHistHitsITS; 
  THnSparseD  *fTrackHistRefitITS; 
  THnSparseD  *fTrackHistChi2ITS; 
  THnSparseD  *fTrackHistHitsSPD; 
  THnSparseD  *fTrackHistDCAz; 
  THnSparseD  *fTrackHistDCAxy; 
  THnSparseD  *fTrackHistChi2TPCITS;   
  THnSparseD  *fTrackHistGeoLengthTPC;
  THnSparseD  *fTrackHistGeoNcrTPC;
  THnSparseD  *fTrackHistGeoNclTPC;
  THnSparseD  *fTrackHistNclTPC;
  THnSparseD *fEventHist;    //-> zV:multMB:multRec:multTrue  
  THnSparseD *fTrackHist;    //-> zV:pt:eta:phi:charge:sigma(1/pt):multMB:mc_delta_pt
  THnSparseD *fTrackHist2;   //-> zv:pt:eta:phi:charge
  TH2D       *fPtHist;      //-> 1pt:sigma1pt
  TH2D       *fPtHist2;      //-> 1pt:sigma1pt (zoomed)
  THnSparseD *fMCTrackHist;  //-> mczV:mcpt:mcta:mcphi:charge:sigma(1/pt):multMB:mc_delta_pt
  THnSparseD *fMCTrackHist2;  //-> pt:ptmc:delta(1/pt)/sigma(1/pt)
  Bool_t      fUseMCInfo;
  TObjArray * fOutputContainer; //! output data container
    
  AliAnalysisTaskCutTest(const AliAnalysisTaskCutTest&); // not implemented
  AliAnalysisTaskCutTest& operator=(const AliAnalysisTaskCutTest&); // not implemented
  
  ClassDef(AliAnalysisTaskCutTest, 5); 
};

#endif
