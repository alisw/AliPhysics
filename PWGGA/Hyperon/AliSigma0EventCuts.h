#ifndef AliSigma0EventCuts_H
#define AliSigma0EventCuts_H

#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliV0ReaderV1.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TObject.h"

#include "AliEventCuts.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TProfile.h"

class AliSigma0EventCuts : public TObject {
 public:
  AliSigma0EventCuts();
  AliSigma0EventCuts(const AliSigma0EventCuts &);
  AliSigma0EventCuts &operator=(const AliSigma0EventCuts &);
  virtual ~AliSigma0EventCuts();

  bool EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent);
  bool EventIsSelectedAliEventCuts(AliVEvent *fInputEvent, AliVEvent *fMCEvent);
  bool EventIsSelectedCuts(AliVEvent *fInputEvent, AliVEvent *fMCEvent);

  void SetUseAliEventCuts(bool useAliEventCuts) {
    fUseAliEventCuts = useAliEventCuts;
  }
  void SetUseConversionCuts(bool useConversionCuts) {
    fUseConversionCuts = useConversionCuts;
  }

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetZVertexCut(float z) {
    fVertex = z;
    if (fAnaUtils) fAnaUtils->SetMaxVtxZ(z);
  }
  void SetNVertexContributors(int nContrib) { fNVertexContributors = nContrib; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }
  void FillTriggerHisto(TH1F *histo);

  UInt_t GetTrigger() const { return fTrigger; }

  void InitCutHistograms();
  TList *GetCutHistograms() const { return fHistograms; }

  AliEventCuts fAliEventCuts;

 protected:
  TList *fHistograms;
  TList *fQA;
  AliV0ReaderV1
      *fV0Reader;         //! basic photon Selection Task !!! needs //! for grid
  TString fV0ReaderName;  //

  bool fUseAliEventCuts;
  bool fUseConversionCuts;

  float fVertex;
  int fNVertexContributors;
  float fV0PercentileMax;
  UInt_t fTrigger;

  TProfile *fHistCuts;                      //
  TH1F *fHistEventQA;                       //
  TH1F *fHistTrigger;                       //
  TH1F *fHistXVertex;                       //
  TH1F *fHistYVertex;                       //
  TH1F *fHistZVertex;                       //
  TH1F *fHistXVertexAfter;                  //
  TH1F *fHistYVertexAfter;                  //
  TH1F *fHistZVertexAfter;                  //
  TH1F *fHistVertexSeparation;              //
  TH1F *fHistSPDresolution;                 //
  TH2I *fHistSPDTracklets;                  //
  TH2I *fHistSPDTrackletsAfter;             //
  TProfile *fHistRunNumber;                 //
  TH1F *fHistCentralityProfile;             //
  TH1F *fHistCentralityProfileAfter;        //
  TH1F *fHistCentralityProfileCoarseAfter;  //

  AliAnalysisUtils
      *fAnaUtils;  //! Object to use analysis utils like pile-up rejection
  AliInputEventHandler *fInputHandler;

 private:
  ClassDef(AliSigma0EventCuts, 1)
};

#endif
