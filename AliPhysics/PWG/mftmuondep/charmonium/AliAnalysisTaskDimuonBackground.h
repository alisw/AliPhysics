#ifndef AliAnalysisTaskDimuonBackground_H
#define AliAnalysisTaskDimuonBackground_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliMultiInputEventHandler.h"
#include "AliMixInputEventHandler.h"
#include "TH1D.h"
#include "TList.h"
#include "AliAODTrack.h"
#include "THnSparse.h"

//====================================================================================================================================================

class  AliAnalysisTaskDimuonBackground : public AliAnalysisTaskSE {

public:

  enum {kGenerated, kReconstructed};
  enum {kSingleEvents, kMixedEvents};
 
  AliAnalysisTaskDimuonBackground();
  AliAnalysisTaskDimuonBackground(const char *name);

  virtual ~AliAnalysisTaskDimuonBackground() { if (fHistogramList) delete fHistogramList; }
  
  void SetVertexMode(Int_t vertexMode) { fVertexMode = vertexMode; }
  void SetVtxResolutionITS(Double_t sigmaX, Double_t sigmaY, Double_t sigmaZ) {
    fVtxResolutionITS[0] = sigmaX;
    fVtxResolutionITS[1] = sigmaY;
    fVtxResolutionITS[2] = sigmaZ;
  }

  void SetMinTriggerMatch(Int_t minTriggerMatch) { fMinTriggerMatch = minTriggerMatch; }

  void SetSingleMuonMinPt(Double_t minPt)   { fSingleMuonMinPt  = minPt;  }
  void SetSingleMuonMinEta(Double_t minEta) { fSingleMuonMinEta = minEta; }
  void SetSingleMuonMaxEta(Double_t maxEta) { fSingleMuonMaxEta = maxEta; }
  void SetSingleMuonMaxChi2(Double_t maxChi2) { fSingleMuonMaxChi2 = maxChi2; }

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void UserExecMix(Option_t *);
  virtual void Terminate(Option_t *);

  AliVEvent* GetMainEvent();
  AliVEvent* GetMixedEvent(Int_t buffId=0);
  
  AliMultiInputEventHandler* SetMainInputHandler(AliAnalysisManager *mgr);
  AliMixInputEventHandler* SetMixingInputHandler(AliMultiInputEventHandler *mainIH);

  Bool_t IsSingleMuonCutPassed(AliAODTrack *mu);

private:

  Double_t fMassJpsi;

  Double_t fPrimaryVertexTrue[3], fPrimaryVertexMixEvTrue[3], fPrimaryVertex[3], fPrimaryVertexMixEv[3], fVtxResolutionITS[3];
  Int_t fVertexMode;

  Int_t fMinTriggerMatch;
  Double_t fSingleMuonMinEta, fSingleMuonMaxEta, fSingleMuonMinPt, fSingleMuonMaxChi2;

  TList *fHistogramList;

  THnSparse *fHistSingleMuonsChi2VsPtVsRapidityBeforeCut;                   //!
  THnSparse *fHistSingleMuonsChi2VsPtVsRapidityAfterCut;                    //!
  THnSparse *fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[2];           //!
  THnSparse *fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances;   //!

  AliMultiInputEventHandler *fMainInputHandler;    //! tmp pointer to main input handler
  AliMixInputEventHandler   *fMixingInputHandler;  //! tmp pointer to mixing input handler
  
  ClassDef(AliAnalysisTaskDimuonBackground, 1) // example of analysis

};

//====================================================================================================================================================

#endif
