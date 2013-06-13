#ifndef ALIANLYSISTASKGAMMACONVV1_cxx
#define ALIANLYSISTASKGAMMACONVV1_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisManager.h"
#include "TH3.h"
#include "TH3F.h"

class AliAnalysisTaskGammaConvV1 : public AliAnalysisTaskSE {
 public:

   AliAnalysisTaskGammaConvV1();
   AliAnalysisTaskGammaConvV1(const char *name);
   virtual ~AliAnalysisTaskGammaConvV1();

   virtual void   UserCreateOutputObjects();
   virtual Bool_t Notify();
   virtual void   UserExec(Option_t *);
   virtual void   Terminate(const Option_t*);
   void InitBack();

   void SetIsHeavyIon(Int_t flag){
      if (flag == 1 || flag ==2 ){
         fIsHeavyIon = 1;    
      } else {
         fIsHeavyIon = 0;    
      }
   }

   void SetIsMC(Bool_t isMC){fIsMC=isMC;}
   void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
   void SetDoMesonQA(Bool_t flag){fDoMesonQA = flag;}
   void SetDoPhotonQA(Bool_t flag){fDoPhotonQA = flag;}
   void ProcessPhotonCandidates();
   void CalculatePi0Candidates();
   void CalculateBackground();
   void CalculateBackgroundRP();
   void ProcessMCParticles();
   void ProcessAODMCParticles();
   void RelabelAODPhotonCandidates(Bool_t mode);
   void ProcessTruePhotonCandidates( AliAODConversionPhoton* TruePhotonCandidate);
   void ProcessTruePhotonCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate);
   void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
   void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
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
   void CountTracks();
   void FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[]);
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
   TList **fMotherRapList;
   TList **fTrueList;
   TList **fTrueMotherRapList;
   TList **fMCList;
   TList **fHeaderNameList;
   TList **fTriggerNameList;
   TList *fOutputContainer;
   TClonesArray *fReaderGammas;
   TList *fGammaCandidates;
   TList *fCutArray;
   AliConversionCuts *fConversionCuts;
   TList *fMesonCutArray;
   AliConversionMesonCuts *fMesonCuts;
   TH1F **hESDConvGammaPt;
   TH1F **hESDConvGammaR;
   TH2F **hESDMotherInvMassPt;
   THnSparseF **sESDMotherInvMassPtZM;
   THnSparseF **sESDMotherInvMassPtY;
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
   TH1F **hMCPi0WOWeightPt;
   TH1F **hMCEtaPt;
   TH1F **hMCPi0InAccPt;
   TH1F **hMCEtaInAccPt;
   TH2F **hMCPi0PtY;
   TH2F **hMCEtaPtY;
   TH1F **hMCK0sPt;
   TH1F **hMCK0sWOWeightPt;
   TH2F **hMCK0sPtY;
   TH2F **hESDTrueMotherInvMassPt;
   TH2F **hESDTruePrimaryMotherInvMassPt;
   TH2F **hESDTruePrimaryPi0MCPtResolPt;
   TH2F **hESDTruePrimaryEtaMCPtResolPt;
   THnSparseF **sESDTruePrimaryMotherInvMassPtY;
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
   TH1I **hNEvents;
   TH1I **hNGoodESDTracks;
   TH1I **hNGammaCandidates;
   TH1I **hNV0Tracks;
   TH1F **hEtaShift;
   
   TRandom3 fRandom;
   Int_t fnGammaCandidates;
   Double_t *fUnsmearedPx; //[fnGammaCandidates]
   Double_t *fUnsmearedPy; //[fnGammaCandidates]
   Double_t *fUnsmearedPz; //[fnGammaCandidates]
   Double_t *fUnsmearedE;  //[fnGammaCandidates]
   Int_t *fMCStackPos;     //[fnGammaCandidates]
   Int_t *fMCStackNeg;     //[fnGammaCandidates]
   Int_t fnCuts;
   Int_t fiCut;
   Int_t fNumberOfESDTracks;
   Bool_t fMoveParticleAccordingToVertex;
   Bool_t fIsHeavyIon;
   Bool_t fDoMesonAnalysis;
   Bool_t fDoMesonQA;
   Bool_t fDoPhotonQA;
   Bool_t fIsFromMBHeader;
   Bool_t fIsMC;

private:

   AliAnalysisTaskGammaConvV1(const AliAnalysisTaskGammaConvV1&); // Prevent copy-construction
   AliAnalysisTaskGammaConvV1 &operator=(const AliAnalysisTaskGammaConvV1&); // Prevent assignment


   ClassDef(AliAnalysisTaskGammaConvV1, 7);
};

#endif
