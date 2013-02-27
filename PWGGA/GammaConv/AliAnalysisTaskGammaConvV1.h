#ifndef ALIANLYSISTASKGAMMACONVV1_cxx
#define ALIANLYSISTASKGAMMACONVV1_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionMesonCuts.h"
#include "TH3.h"
#include "TH3F.h"

class AliAnalysisTaskGammaConvV1 : public AliAnalysisTaskSE {
 public:

   AliAnalysisTaskGammaConvV1();
   AliAnalysisTaskGammaConvV1(const char *name);
   virtual ~AliAnalysisTaskGammaConvV1();

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *);
   virtual void   Terminate(const Option_t*);
   void InitBack();

   void SetIsHeavyIon(Bool_t flag){fIsHeavyIon = flag;}
   void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
   void ProcessPhotonCandidates();
   void CalculatePi0Candidates();
   void CalculateBackground();
   void CalculateBackgroundRP();
   void ProcessMCParticles();
   void ProcessTruePhotonCandidates( AliAODConversionPhoton* TruePhotonCandidate);
   void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
   void RotateParticle(AliAODConversionPhoton *gamma);
   void SetConversionCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fCutArray = CutArray;
   }
   void SetMesonCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fMesonCutArray = CutArray;
   }

   // BG HandlerSettings
   void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
   void CountESDTracks();
   void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
   void UpdateEventByEventData();

 protected:
   AliV0ReaderV1 *fV0Reader;
   AliGammaConversionAODBGHandler **fBGHandler;
   AliConversionAODBGHandlerRP    **fBGHandlerRP;
   AliVEvent *fInputEvent;
   AliMCEvent *fMCEvent;
   AliStack *fMCStack;
   TList **fCutFolder;
   TList **fESDList;
   TList **fBackList;
   TList **fMotherList;
   TList **fTrueList;
   TList **fMCList;
   TList **fHeaderNameList;
   TList *fOutputContainer;
   TClonesArray *fReaderGammas;
   TList *fGoodGammas;
   TList *fCutArray;
   AliConversionCuts *fConversionCuts;
   TList *fMesonCutArray;
   AliConversionMesonCuts *fMesonCuts;
   TH1F **hESDConvGammaPt;
   TH1F **hESDConvGammaR;
   TH2F **hESDMotherInvMassPt;
   THnSparseF **sESDMotherInvMassPtZM;
   TH2F **hESDMotherBackInvMassPt;
   THnSparseF **sESDMotherBackInvMassPtZM;
   TH2F **hESDMotherInvMassEalpha;
   TH1F **hMCAllGammaPt;
   TH1F **hMCDecayGammaPi0Pt;
   TH1F **hMCDecayGammaRhoPt;
   TH1F **hMCDecayGammaEtaPt;
   TH1F **hMCDecayGammaOmegaPt;
   TH1F **hMCDecayGammaEtapPt;
   TH1F **hMCDecayGammaPhiPt;
   TH1F **hMCDecayGammaSigmaPt;
   TH1F **hMCConvGammaPt;
   TH1F **hMCConvGammaR;
   TH1F **hMCConvGammaEta;
   TH1F **hMCConvGammaRSPt;
   TH1F **hMCConvGammaRSR;
   TH1F **hMCConvGammaRSEta;
   TH1F **hMCPi0Pt;
   TH1F **hMCEtaPt;
   TH1F **hMCPi0InAccPt;
   TH1F **hMCEtaInAccPt;
   TH2F **hESDTrueMotherInvMassPt;
   TH2F **hESDTruePi0FromEtaInvMassPt;
   TH2F **hESDTruePrimaryMotherInvMassMCPt;
   TH2F **hESDTruePrimaryPi0ESDPtMCPt;
   TH2F **hESDTrueSecondaryMotherInvMassPt;
   TH2F **hESDTrueSecondaryMotherFromK0sInvMassPt;
   TH1F **hESDTrueK0sWithPi0DaughterMCPt;
   TH2F **hESDTrueSecondaryMotherFromEtaInvMassPt;
   TH1F **hESDTrueEtaWithPi0DaughterMCPt;
   TH2F **hESDTrueBckGGInvMassPt;
   TH2F **hESDTrueBckContInvMassPt;
   TH2F **hESDTrueMotherDalitzInvMassPt;
   TH1F **hESDTrueConvGammaPt;
   TH2F **hESDCombinatorialPt;
   TH1F **hESDTruePrimaryConvGammaPt;
   TH1F **hESDTruePrimaryConvGammaR;
   TH1F **hESDTruePrimaryConvGammaEta;
   TH2F **hESDTruePrimaryConvGammaESDPtMCPt;
   TH2F **hESDTruePrimaryConvGammaRSESDPtMCPt;
   TH1F **hESDTrueSecondaryConvGammaPt;
   TH1F **hESDTrueSecondaryConvGammaR;
   TH1F **hESDTrueSecondaryConvGammaFromXFromK0sPt;
   TH1F **hESDPi0Alpha;
   TH1F **hESDBackAlpha;
   TH1F **hESDTruePi0Alpha;
   TH1I **hNEvents;
   TH1I **hNGoodESDTracks;
   TH1I **hNV0Tracks;

   TRandom3 fRandom;
   Int_t fnGoodGammas;
   Double_t *fUnsmearedPx; //[fnGoodGammas]
   Double_t *fUnsmearedPy; //[fnGoodGammas]
   Double_t *fUnsmearedPz; //[fnGoodGammas]
   Double_t *fUnsmearedE;  //[fnGoodGammas]
   Int_t fnCuts;
   Int_t fiCut;
   Int_t fNumberOfESDTracks;
   Bool_t fMoveParticleAccordingToVertex;
   Bool_t fIsHeavyIon;
   Bool_t fDoMesonAnalysis;
   Bool_t fIsFromBGEvent;

   ClassDef(AliAnalysisTaskGammaConvV1, 2);
};

#endif
