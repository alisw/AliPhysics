// event plane calibration
// author: Ionut-Cristian Arsene, i.c.arsene@gsi.de
//         2012/Apr/06

//#include "AliSysInfo.h"
#ifndef ALIREDUCEDANALYSISTASKTWOPWRTRP_H
#define ALIREDUCEDANALYSISTASKTWOPWRTRP_H

#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
//#include "AliAnalysisTaskSE.h"

#include "AliReducedAnalysisTaskSE.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliCMEAnalysisCuts.h"
//#include "AliCMEVarManager.h"
#include "AliChargeOnePwrtRP.h"
#include "AliChargeTwoPwrtRP.h"
//class AliAnalysis;
//class AliReducedEvent;
//class AliCMEVarManager;
//class AliQnCorrectionsQnVector;
//class AliCMEAnalysisCuts;
//class AliReducedHelper;


//________________________________________________________________
class AliReducedAnalysisTaskTwoPwrtRP : public AliReducedAnalysisTaskSE {

public:
  AliReducedAnalysisTaskTwoPwrtRP();
  AliReducedAnalysisTaskTwoPwrtRP(const char *name, const Char_t* title);
  virtual ~AliReducedAnalysisTaskTwoPwrtRP();

  // Virtual functions, to be implemented in the inheriting analysis classes
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();

  // setters

  // getters


  enum Constants {
      maxNBins=102
  };



  void SetFriendPath(TString path) {fFriendPath=path;}

  //Bool_t IsEventSelected(AliReducedEvent* event);
  //Bool_t TriggerSelected(AliReducedEvent* event);
  void SetEventCuts(AliCMEAnalysisCuts* cuts)  {fEventCuts = cuts;}
  void SetBinEdgePt(   Int_t bin, Float_t edge)  {fPtbinning[bin]  =edge;  fNPtbins++  ;}
  void SetBinEdgeEta(  Int_t bin, Float_t edge)  {fEtabinning[bin] =edge;  fNEtabins++ ;}
  void SetBinEdgeDPt(  Int_t bin, Float_t edge)  {fDPtbinning[bin] =edge;  fNDPtbins++ ;}
  void SetBinEdgeDEta( Int_t bin, Float_t edge)  {fDEtabinning[bin]=edge;  fNDEtabins++;}
  void SetBinEdgeMPt(  Int_t bin, Float_t edge)  {fMPtbinning[bin] =edge;  fNMPtbins++ ;}
  void SetBinEdgeMEta( Int_t bin, Float_t edge)  {fMEtabinning[bin]=edge;  fNMEtabins++;}
  void SetBinEdgeCt(   Int_t bin, Float_t edge)  {fCtbinning[bin]  =edge;  fNCtbins++  ;}
  void Set1PwrtRP(AliChargeOnePwrtRP* vec,AliCMEAnalysisCuts* p1) {fFlow[fNOnePwrtRP]=vec;fTrackCutsOnePwrtRP[fNOnePwrtRP] = p1;fNOnePwrtRP++;if(fMinHarmonic>vec->GetMinHarmonic()) fMinHarmonic=vec->GetMinHarmonic(); if(fMaxHarmonic<vec->GetMaxHarmonic()) fMaxHarmonic=vec->GetMaxHarmonic();}
  void Set2PwrtRP(AliChargeTwoPwrtRP* vec,AliCMEAnalysisCuts* p1,AliCMEAnalysisCuts* p2=0x0) {fChargeCorrelation[fNTwoPwrtRP]=vec;fTrackCutsTwoPwrtRP[fNTwoPwrtRP][0] = p1;fTrackCutsTwoPwrtRP[fNTwoPwrtRP][1] = p2;fNTwoPwrtRP++;}
  void SetAxes(TAxis* Ct,TAxis* Pt,TAxis* mPt,TAxis* dPt,TAxis* Eta,TAxis* mEta,TAxis* dEta) {fAxes[0]=new TAxis(*Ct);fAxes[1]=new TAxis(*Pt);fAxes[2]=new TAxis(*mPt);fAxes[3]=new TAxis(*dPt);fAxes[4]=new TAxis(*Eta);fAxes[5]=new TAxis(*mEta);fAxes[6]=new TAxis(*dEta);}
  TString FriendPath()       const {return fFriendPath;}
  void SetTPCstandaloneTracks(Bool_t b=kTRUE)  {fTracksTPCstandalone=b;}
  void SetTPCtrackcuts(AliCMEAnalysisCuts* cuts) {fTPCtrackcuts=cuts;}
  void SetGlobaltrackcuts(AliCMEAnalysisCuts* cuts) {fGlobaltrackcuts=cuts;}
  void SetRemoveOutliers(Bool_t b=kTRUE) {fRemoveOutliers=b;}

  TAxis* GetAxis(Int_t i) const {return fAxes[i];}

  //getters
  Bool_t IsTracksTPCstandalone() const {return fTracksTPCstandalone;}
  TList* GetOutputList() {return &fListHistos; };
  Int_t GetNbinsPt(  )  {return  fNPtbins  ;}
  Int_t GetNbinsEta( )  {return  fNEtabins ;}
  Int_t GetNbinsDPt( )  {return  fNDPtbins ;}
  Int_t GetNbinsDEta()  {return  fNDEtabins;}
  Int_t GetNbinsMPt( )  {return  fNMPtbins ;}
  Int_t GetNbinsMEta()  {return  fNMEtabins;}
  Int_t GetNbinsCt(  )  {return  fNCtbins  ;}

  Int_t GetNOnePwrtRP() const {return fNOnePwrtRP;}
  Int_t GetNTwoPwrtRP() const {return fNTwoPwrtRP;}
  AliCMEAnalysisCuts* GetOnePwrtRP_TrackCuts(Int_t i) const {return fTrackCutsOnePwrtRP[i];}


protected:

  virtual Bool_t IsEventSelected(AliReducedBaseEvent* event) {kTRUE;}
  virtual Bool_t IsTrackSelected(AliReducedBaseTrack* track) {kTRUE;}
  virtual Bool_t IsPairSelected(AliReducedBaseTrack* pair)   {kTRUE;}

  void MakeBinMapping();

 private:

  Bool_t fInitialized;
  TString fFriendPath;
  TList fListHistos;                 //! Output histograms
  TList fListHistosQn;               //! Qn correlation histograms
  Int_t fOffset;
  Int_t fNfiles;
  AliCMEAnalysisCuts * fEventCuts;
  AliQnCorrectionsQnVector* fQvecTPC;
  AliQnCorrectionsQnVector* fQvecVZEROA;
  AliQnCorrectionsQnVector* fQvecVZEROC;
  AliQnCorrectionsQnVector* fQvecFMDA;
  AliQnCorrectionsQnVector* fQvecFMDC;
  AliQnCorrectionsQnVector* fQvectors[AliChargeOnePwrtRP::Nqvectors];
  Float_t fValues[2000];
  AliChargeOnePwrtRP* fFlow[10];
  AliChargeTwoPwrtRP* fChargeCorrelation[10];
  Int_t fNOnePwrtRP;
  Int_t fNTwoPwrtRP;
  AliCMEAnalysisCuts* fTrackingQuality[10];
  AliCMEAnalysisCuts* fTrackCutsOnePwrtRP[10];
  AliCMEAnalysisCuts* fTrackCutsTwoPwrtRP[10][2];
  AliCMEAnalysisCuts* fTPCtrackcuts;      // to apply our version of Mikolaj's cut
  AliCMEAnalysisCuts* fGlobaltrackcuts;   // to apply our version of Mikolaj's cut
  //Float_t fTrackSelectedMapX[10000][6];
  //Float_t fTrackSelectedMapY[10000][6];
  Float_t fCtbinning[maxNBins];
  Float_t fPtbinning[maxNBins];
  Float_t fEtabinning[maxNBins];
  Float_t fDPtbinning[maxNBins];
  Float_t fMPtbinning[maxNBins];
  Float_t fDEtabinning[maxNBins];
  Float_t fMEtabinning[maxNBins];
  typedef std::map<int, int> Map;
  Map fPtBinningMap;
  Map fMPtBinningMap;
  Map fDPtBinningMap;
  Map fEtaBinningMap;
  Map fDEtaBinningMap;
  Map fMEtaBinningMap;

  //Int_t fMapCtbinning[maxNBins*100];
  //Int_t fMapPtbinning[maxNBins*100];
  //Int_t fMapEtabinning[maxNBins*100];
  //Int_t fMapDPtbinning[maxNBins*100];
  //Int_t fMapMPtbinning[maxNBins*100];
  //Int_t fMapDEtabinning[maxNBins*100];
  //Int_t fMapMEtabinning[maxNBins*100];
  Int_t fNCtbins;
  Int_t fNPtbins;
  Int_t fNEtabins;
  Int_t fNDPtbins;
  Int_t fNDEtabins;
  Int_t fNMPtbins;
  Int_t fNMEtabins;
  TAxis* fAxes[7];
  Int_t fMinHarmonic;
  Int_t fMaxHarmonic;
  //TProfile* fPcent;
  //TProfile* fPpt;
  //TProfile* fPeta;
  //AliReducedAnalysisTaskTwoPwrtRP* fTask;
  EventPlaneCorrelation* fEPcor[AliChargeOnePwrtRP::Nqvectors];
  Bool_t fTracksTPCstandalone;
  Bool_t fRemoveOutliers;
  Int_t fNevents;
  THnF* fTrackDistribution;

  TProfile* fProfileCent[4][4][2];
  TProfile* fProfileMPt[9][4][4][2];
  THnF* fNDPhiDEta[9];


  AliReducedAnalysisTaskTwoPwrtRP(const AliReducedAnalysisTaskTwoPwrtRP &c);
  AliReducedAnalysisTaskTwoPwrtRP& operator= (const AliReducedAnalysisTaskTwoPwrtRP &c);



  ClassDef(AliReducedAnalysisTaskTwoPwrtRP, 1);
};

#endif
