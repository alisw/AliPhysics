#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOTREELPHI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOTREELPHI_H_
#include "Rtypes.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamControlSample.h"
#include "TTree.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoTreeLPhi : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskNanoTreeLPhi();
  AliAnalysisTaskNanoTreeLPhi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoTreeLPhi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };

  void SetLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fLambdaCuts = cuts;
  }
  void SetAntiLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fAntiLambdaCuts = cuts;
  }
  void SetPosKaonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosKaonCuts = trkCuts;
  }
  void SetNegKaonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegKaonCuts = trkCuts;
  }
  void SetPhiCuts(AliFemtoDreamv0Cuts *phiCuts) { fPhiCuts = phiCuts; }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config)
  {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }

  Bool_t FillLambda(AliFemtoDreamv0* TheV0);

private:
  AliAnalysisTaskNanoTreeLPhi(const AliAnalysisTaskNanoTreeLPhi &);
  AliAnalysisTaskNanoTreeLPhi &operator=(const AliAnalysisTaskNanoTreeLPhi &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  float fInvMassCutSBdown;                //
  float fInvMassCutSBup;                  //
  UInt_t fTrigger;                        //
  TList *fResults;                        //!
  TList *fResultsQA;                      //!
  TList *fQA;                             //!
  TList *fEvtList;                        //!
  TList *fLambdaList;                     //!
  TList *fLambdaMCList;                   //!
  TList *fAntiLambdaList;                 //!
  TList *fAntiLambdaMCList;               //!
  TList *fKaonPlusList;                   //!
  TList *fKaonPlusMCList;                 //!
  TList *fKaonMinusList;                  //!
  TList *fKaonMinusMCList;                //!
  TList *fPhiList;                        //!
  TList *fPhiMCList;                      //!
  AliVEvent *fInputEvent;                 //! current event
  AliFemtoDreamEvent *fEvent;             //!
  AliFemtoDreamv0 *fLambda;               //!
  AliFemtoDreamv0 *fPhiParticle;          //!
  AliFemtoDreamTrack *fTrack;             //!
  AliFemtoDreamEventCuts *fEventCuts;     //
  AliFemtoDreamTrackCuts *fPosKaonCuts;   //
  AliFemtoDreamTrackCuts *fNegKaonCuts;   //
  AliFemtoDreamv0Cuts *fPhiCuts;          //
  AliFemtoDreamv0Cuts *fLambdaCuts;       //
  AliFemtoDreamv0Cuts *fAntiLambdaCuts;   //
  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamControlSample *fSample;    //!
  AliVTrack **fGTI;                       //!
  int fTrackBufferSize;                   //

  TTree *fTree;                           //!
  Int_t fTRunNumber;
  Float_t fTVz;
  Int_t fTMult;
  Float_t fTSpher;

  static const Int_t nMaxLambda = 100;
  Int_t fNumLambda;

  Float_t v0_pT[nMaxLambda];
  Float_t v0_decvtx[nMaxLambda];
  Float_t v0_tranrad[nMaxLambda];
  Float_t v0Daugh_dcaDecvtx[nMaxLambda];
  Float_t v0_cpa[nMaxLambda];
  Float_t v0_mass[nMaxLambda];
  ///K0s rejection
  ///IM cut around nom.mass
  Float_t Daugh_eta[nMaxLambda];
  Int_t Daugh_nTpcCls[nMaxLambda];
  Float_t Daugh_dca[nMaxLambda];
  Float_t Daugh_nSigma[nMaxLambda];

  ClassDef(AliAnalysisTaskNanoTreeLPhi, 1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOTREELPHI_H_ */
