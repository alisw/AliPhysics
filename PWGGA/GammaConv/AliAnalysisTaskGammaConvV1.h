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
#include "TProfile2D.h"
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
      fIsHeavyIon = flag;    
      
   }

   void SetIsMC(Bool_t isMC){fIsMC=isMC;}
   void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
   void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
   void SetDoPhotonQA(Int_t flag){fDoPhotonQA = flag;}
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
   void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
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
   void FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[]);
   void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
   void UpdateEventByEventData();
   void SetLogBinningXTH2(TH2* histoRebin);
   
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
   TList **fPhotonDCAList;
   TList **fMesonDCAList;        
   TList **fTrueList;
   TList **fMCList;
   TList **fHeaderNameList;
   TList *fOutputContainer;
   TClonesArray *fReaderGammas;
   TList *fGammaCandidates;
   TList *fCutArray;
   AliConversionCuts *fConversionCuts;
   TList *fMesonCutArray;
   AliConversionMesonCuts *fMesonCuts;
   TH1F **hESDConvGammaPt;
   TH1F **hESDConvGammaR;
   TH1F **hESDConvGammaEta;
   TTree **tESDConvGammaPtDcazCat;
   Float_t fPtGamma;
   Float_t fDCAzPhoton;
   Float_t fRConvPhoton;
   Float_t fEtaPhoton;
   UChar_t iCatPhoton;
   UChar_t iPhotonMCInfo; // 0: garbage,
                         // 1: background
                         // 2: secondary photon not from eta or k0s,
                         // 3: secondary photon from eta, 
                         // 4: secondary photon from k0s, 
                         // 5: dalitz
                         // 6: primary gamma
   TH2F **hESDMotherInvMassPt;
   THnSparseF **sESDMotherInvMassPtZM;
   TH2F **hESDMotherBackInvMassPt;
   THnSparseF **sESDMotherBackInvMassPtZM;
   TH2F **hESDMotherInvMassEalpha;
   TH2F **hESDMotherPi0PtY;
   TH2F **hESDMotherEtaPtY;
   TH2F **hESDMotherPi0PtAlpha;
   TH2F **hESDMotherEtaPtAlpha;
   TH2F **hESDMotherPi0PtOpenAngle;
   TH2F **hESDMotherEtaPtOpenAngle;
   TH1I **hMCHeaders;
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
   TH1F **hMCPi0Pt;
   TH1F **hMCPi0WOWeightPt;
   TH1F **hMCEtaPt;
   TH1F **hMCEtaWOWeightPt;
   TH1F **hMCPi0InAccPt;
   TH1F **hMCEtaInAccPt;
   TH2F **hMCPi0PtY;
   TH2F **hMCEtaPtY;
   TH1F **hMCK0sPt;
   TH1F **hMCK0sWOWeightPt;
   TH2F **hMCK0sPtY;
   TH2F **hESDTrueMotherInvMassPt;
   TH2F **hESDTruePrimaryMotherInvMassPt;
   TH2F **hESDTruePrimaryMotherW0WeightingInvMassPt;
   TProfile2D **pESDTruePrimaryMotherWeightsInvMassPt;
   TH2F **hESDTruePrimaryPi0MCPtResolPt;
   TH2F **hESDTruePrimaryEtaMCPtResolPt;
   TH2F **hESDTrueSecondaryMotherInvMassPt;
   TH2F **hESDTrueSecondaryMotherFromK0sInvMassPt;
   TH1F **hESDTrueK0sWithPi0DaughterMCPt;
   TH2F **hESDTrueSecondaryMotherFromEtaInvMassPt;
   TH1F **hESDTrueEtaWithPi0DaughterMCPt;
   TH2F **hESDTrueBckGGInvMassPt;
   TH2F **hESDTrueBckContInvMassPt;
   TH2F **hESDTruePi0PtY;
   TH2F **hESDTrueEtaPtY;
   TH2F **hESDTruePi0PtAlpha;
   TH2F **hESDTrueEtaPtAlpha;
   TH2F **hESDTruePi0PtOpenAngle;
   TH2F **hESDTrueEtaPtOpenAngle;
   TH2F **hESDTrueMotherDalitzInvMassPt;
   TH1F **hESDTrueConvGammaPt;
   TH1F **hESDTrueConvGammaEta;
   TH2F **hESDCombinatorialPt;
   TH1F **hESDTruePrimaryConvGammaPt;
   TH2F **hESDTruePrimaryConvGammaESDPtMCPt;
   TH1F **hESDTrueSecondaryConvGammaPt;
   TH1F **hESDTrueSecondaryConvGammaFromXFromK0sPt;
   TH1F **hESDTrueSecondaryConvGammaFromXFromLambdaPt;
   TH2F **hESDTrueDalitzPsiPairDeltaPhi;
   TH2F **hESDTrueGammaPsiPairDeltaPhi;
   TH1I **hNEvents;
   TH1I **hNGoodESDTracks;
   TH1I **hNGammaCandidates;
   TH1I **hNV0Tracks;
   TProfile **hEtaShift;
   TTree **tESDMesonsInvMassPtDcazMinDcazMaxFlag;
   Float_t fInvMass;
   Float_t fPt;
   Float_t fDCAzGammaMin;
   Float_t fDCAzGammaMax;
   UChar_t iFlag;
   UChar_t iMesonMCInfo; // 0: garbage,
                         // 1: background
                         // 2: secondary meson not from eta or k0s,
                         // 3: secondary meson from eta, 
                         // 4: secondary meson from k0s, 
                         // 5: dalitz
                         // 6: primary meson gamma-gamma-channel
   Double_t fEventPlaneAngle; // EventPlaneAngle
   TRandom3 fRandom;
   Int_t fnGammaCandidates;
   Double_t *fUnsmearedPx; //[fnGammaCandidates]
   Double_t *fUnsmearedPy; //[fnGammaCandidates]
   Double_t *fUnsmearedPz; //[fnGammaCandidates]
   Double_t *fUnsmearedE;  //[fnGammaCandidates]
   Int_t *fMCStackPos;     //[fnGammaCandidates]
   Int_t *fMCStackNeg;     //[fnGammaCandidates]
   Int_t *fESDArrayPos;    //[fnGammaCandidates]
   Int_t *fESDArrayNeg;    //[fnGammaCandidates]
   Int_t fnCuts;
   Int_t fiCut;
   Bool_t fMoveParticleAccordingToVertex;
   Int_t fIsHeavyIon;
   Bool_t fDoMesonAnalysis;
   Int_t fDoMesonQA;
   Int_t fDoPhotonQA;
   Bool_t fIsFromMBHeader;
   Bool_t fIsMC;

private:

   AliAnalysisTaskGammaConvV1(const AliAnalysisTaskGammaConvV1&); // Prevent copy-construction
   AliAnalysisTaskGammaConvV1 &operator=(const AliAnalysisTaskGammaConvV1&); // Prevent assignment


   ClassDef(AliAnalysisTaskGammaConvV1, 10);
};

#endif
