#ifndef ALIANALYSISTASKITSTPCALIGNMENT_H
#define ALIANALYSISTASKITSTPCALIGNMENT_H

///////////////////////////////////////////////////////////////////////////
//  Class AliAnalysisTaskITSTPCalignment
//  runs ITS-TPC alignment with TPC vdrift calibration
//
//    Origin: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
///////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
class TList;
class TTree;
class AliESDEvent;
class AliESDfriend;
class AliMCEvent;
class AliRelAlignerKalman;
class AliRelAlignerKalmanArray;
class AliExternalTrackParam;
class TH2F;

class AliAnalysisTaskITSTPCalignment : public AliAnalysisTaskSE
{
public:

  enum { kNoESD=0,
         kNoESDfriend=1,
         kNoFriendTrack=2,
         kNoITSoutParams=3,
         kESDfriend=4,
         kFriendsSkipBit=5,
         kFriendTrack=6,
         kITSoutParams=7
       };

  AliAnalysisTaskITSTPCalignment();
  AliAnalysisTaskITSTPCalignment(const char *name);
  virtual ~AliAnalysisTaskITSTPCalignment() {}

  void SetupAlignerArray( Int_t t0, Int_t tend, Int_t slotwidth )
    { fT0=t0; fTend=tend; fSlotWidth=slotwidth; }
  void SetFillDebugTree(Bool_t m=kTRUE) {fFillDebugTree=m;}
  void SetDoQA(Bool_t d=kTRUE) {fDoQA=d;}
  void SetMinPt(Double_t m) {fMinPt=m;}
  void SetMinNclsITS(Int_t m) {fMinPointsVol1=m;}
  void SetMinNclsTPC(Int_t m) {fMinPointsVol2=m;}
  void DoQA(AliExternalTrackParam* paramsITS,
            AliExternalTrackParam* paramsTPC);
  void SetRejectOutliers(Bool_t set=kTRUE){fRejectOutliers=set;}
  void SetRejectOutliersSigma2Median(Bool_t set=kTRUE){fRejectOutliersSigma2Median=set;}
  void SetOutRejSigma(Double_t d){fOutRejSigma=d;}
  void SetOutRejSigma2Median(Double_t d){fOutRejSigma2Median=d;}
  void SetOutRejSigmaOnMerge(Double_t d){fOutRejSigmaOnMerge=d;}
  void SetUseITSoutGlobalTrack(Bool_t b){fUseITSoutGlobalTrack=b;}
  void SetUseITSoutITSSAtrack(Bool_t b){fUseITSoutITSSAtrack=b;}

  Int_t FindMatchingTracks(TObjArray& arrITS, TObjArray& arrTPC, AliESDEvent* pESD);
  void AnalyzeESDevent(AliESDEvent* event);
  
  //AnalysisTaskSE interface methods
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual Bool_t UserNotify();
  void Terminate(Option_t *);

private:
  AliRelAlignerKalmanArray* fArrayITSglobal;   //! array of aligners with ITS global
  AliRelAlignerKalmanArray* fArrayITSsa;   //! array of aligners ITS standalone
  TTree* fDebugTree;                  //! tree
  AliRelAlignerKalman* fAligner;      //! aligner
  TList* fList;                       //! list with QA histograms
 
  Bool_t fFillDebugTree;              // do we write the debug tree?
  Bool_t fDoQA;                       // do QA?
  Int_t fT0;                          // t0
  Int_t fTend;                        // tend
  Int_t fSlotWidth;                   // slotwidth
  Double_t fMinPt;                    // min pt for tracks
  Int_t fMinPointsVol1;               // min clusters its
  Int_t fMinPointsVol2;               // min clusters tpc
  Bool_t fRejectOutliers;             // reject outliers?
  Double_t fOutRejSigma;              // how many outliers
  Bool_t fRejectOutliersSigma2Median; // 
  Double_t fOutRejSigma2Median;       // 
  Double_t fOutRejSigmaOnMerge;       // 
  Bool_t fUseITSoutGlobalTrack;       //
  Bool_t fUseITSoutITSSAtrack;       //

  AliAnalysisTaskITSTPCalignment(const AliAnalysisTaskITSTPCalignment&); // not implemented
  AliAnalysisTaskITSTPCalignment& operator=(const AliAnalysisTaskITSTPCalignment&); // not implemented

  ClassDef(AliAnalysisTaskITSTPCalignment, 3);
};

#endif
