#ifndef ALIANLYSISTASKGAMMACONVV1_cxx
#define ALIANLYSISTASKGAMMACONVV1_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "TH3.h"


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
    void ProcessMCParticles();
    void ProcessTruePhotonCandidates( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);   
    void RotateParticle(AliAODConversionPhoton *gamma);
    void SetConversionCutList(Int_t nCuts, TList *CutArray){
       fnCuts = nCuts;
       fCutArray = CutArray;
    }
    
    // BG HandlerSettings
    void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
    void CountESDTracks();
    void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData();

protected:
    AliV0ReaderV1 *fV0Reader;
    AliGammaConversionAODBGHandler **fBGHandler;
    AliESDEvent *fESDEvent;
    AliMCEvent *fMCEvent;
    AliStack *fMCStack;
    TList **fCutFolder;
    TList **fESDList;
    TList **fBackList;
    TList **fTrueList;
    TList **fMCList;
    TList *fOutputContainer;
    TClonesArray *fReaderGammas;
    TList *fGoodGammas;
    TList *fCutArray;
    AliConversionCuts *fConversionCuts;
    TH1F **hESDConvGammaPt;
    TH2F **hESDMotherInvMassPt;
    THnSparseF **sESDMotherInvMassPtZM;
    TH2F **hESDMotherBackInvMassPt;
    THnSparseF **sESDMotherBackInvMassPtZM;
    TH2F **hESDMotherInvMassEalpha;
    TH1F **hMCAllGammaPt;
	 TH1F **hMCAllGammaPtAddedSig;
    TH1F **hMCDecayGammaPi0Pt;
	 TH1F **hMCDecayGammaPi0PtAddedSig;
    TH1F **hMCDecayGammaRhoPt;
	 TH1F **hMCDecayGammaRhoPtAddedSig;
    TH1F **hMCDecayGammaEtaPt;
	 TH1F **hMCDecayGammaEtaPtAddedSig;
    TH1F **hMCDecayGammaOmegaPt;
	 TH1F **hMCDecayGammaOmegaPtAddedSig;
    TH1F **hMCDecayGammaEtapPt;
	 TH1F **hMCDecayGammaEtapPtAddedSig;
    TH1F **hMCDecayGammaPhiPt;
	 TH1F **hMCDecayGammaPhiPtAddedSig;
    TH1F **hMCConvGammaPt;
	 TH1F **hMCConvGammaPtAddedSig;
    TH1F **hMCPi0Pt;
    TH1F **hMCEtaPt;
    TH1F **hMCPi0InAccPt;
    TH1F **hMCEtaInAccPt;
	 TH1F **hMCPi0PtAddedSig;
    TH1F **hMCEtaPtAddedSig;
    TH1F **hMCPi0InAccPtAddedSig;
    TH1F **hMCEtaInAccPtAddedSig;
    TH2F **hESDTrueMotherInvMassPt;
	 TH2F **hESDTrueMotherInvMassPtAddedSig;
    TH2F **hESDTruePrimaryMotherInvMassMCPt;
	 TH2F **hESDTruePrimaryMotherInvMassMCPtAddedSig;
    TH2F **hESDTruePrimaryPi0ESDPtMCPt;
	 TH2F **hESDTruePrimaryPi0ESDPtMCPtAddedSig;
    TH2F **hESDTrueSecondaryMotherInvMassPt;
	 TH2F **hESDTrueSecondaryMotherInvMassPtAddedSig;
    TH2F **hESDTrueSecondaryMotherFromK0sInvMassPt;
	 TH2F **hESDTrueSecondaryMotherFromK0sInvMassPtAddedSig;
    TH2F **hESDTrueBckGGInvMassPt;
    TH2F **hESDTrueBckContInvMassPt;
    TH2F **hESDTrueMotherDalitzInvMassPt;
	 TH2F **hESDTrueMotherDalitzInvMassPtAddedSig;
    TH1F **hESDTrueConvGammaPt;
	 TH1F **hESDTrueConvGammaPtAddedSig;
    TH1F **hESDTrueElecCombPt;
    TH1F **hESDTruePionCombPt;
    TH1F **hESDTruePrimaryConvGammaPt;
	 TH1F **hESDTruePrimaryConvGammaPtAddedSig;
    TH2F **hESDTruePrimaryConvGammaESDPtMCPt;
	 TH2F **hESDTruePrimaryConvGammaESDPtMCPtAddedSig;
    TH1F **hESDTrueSecondaryConvGammaPt;
	 TH1F **hESDTrueSecondaryConvGammaPtAddedSig;
    TH1F **hESDTrueSecondaryConvGammaFromK0sPt;
	 TH1F **hESDTrueSecondaryConvGammaFromK0sPtAddedSig;
    TH1F **hESDTrueSecondaryConvGammaFromXFromK0sPt;
	 TH1F **hESDTrueSecondaryConvGammaFromXFromK0sPtAddedSig;
    TH1I **hNEvents;
    TH1I **hNGoodESDTracks;
	 TH1I **hNV0Tracks;
        
    TRandom3 fRandom;
    Double_t *fUnsmearedPx;
    Double_t *fUnsmearedPy;
    Double_t *fUnsmearedPz;
    Double_t *fUnsmearedE;
    Int_t fnCuts;
    Int_t fiCut;
    Int_t fNumberOfESDTracks;
    Bool_t fMoveParticleAccordingToVertex;
    Bool_t fIsHeavyIon;
    Bool_t fDoMesonAnalysis;
    
    ClassDef(AliAnalysisTaskGammaConvV1, 2);
};

#endif
