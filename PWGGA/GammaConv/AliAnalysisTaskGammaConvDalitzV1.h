#ifndef ALIANALYSISTASKGAMMACONVDALITZV1_H
#define ALIANALYSISTASKGAMMACONVDALITZV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliDalitzElectronSelector.h"
#include "AliConversionMesonCuts.h"
#include "AliGammaConversionAODBGHandler.h"

class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliV0Reader;
class AliGammaConversionHistograms;
class AliTriggerAnalysis;

class AliAnalysisTaskGammaConvDalitzV1: public AliAnalysisTaskSE
{
  public:

	AliAnalysisTaskGammaConvDalitzV1();
	AliAnalysisTaskGammaConvDalitzV1( const char* name );
	virtual ~AliAnalysisTaskGammaConvDalitzV1();

	virtual void UserExec(Option_t *);
	virtual void UserCreateOutputObjects();
	virtual Bool_t Notify();
	virtual void Terminate(const Option_t *);



         
	void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
         
	void SetIsHeavyIon(Int_t flag){
		if (flag == 1 || flag ==2 ){
			fIsHeavyIon = 1;    
		} else {
			fIsHeavyIon = 0;    
		}
	}
	
	void SetIsMC(Bool_t isMC){fIsMC=isMC;}
	void SetConversionCutList(Int_t nCuts, TList *CutArray){
		fnCuts= nCuts;
		fCutGammaArray = CutArray;
	}
	void SetElectronCutList(TList *CutArray){
		fCutElectronArray = CutArray;
	}
	void SetMesonCutList(TList *CutArray){
		fCutMesonArray = CutArray;
	}
  

	private:

		void InitBack();
		void ProcessPhotonCandidates();
		void ProcessTruePhotonCandidates(AliAODConversionPhoton*);
                void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
                void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    		void ProcessElectronCandidates();
		void ProcessMCParticles();
		void CountESDTracks();
                void CalculatePi0DalitzCandidates();
                void CalculateBackground();
                void UpdateEventByEventData();
                Double_t GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const;
		Bool_t IsDalitz(TParticle *fMCMother) const;
                Bool_t IsPi0DalitzDaughter( Int_t label ) const;
		
		

    AliV0ReaderV1 *fV0Reader;
    AliDalitzElectronSelector* fElecSelector;
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
    vector<Int_t> fSelectorElectronIndex;
    vector<Int_t> fSelectorPositronIndex;
    TList *fGoodGammas;
    TList *fGoodVirtualGammas;
    TList *fGoodElectrons;
    TList *fGoodPositrons;
    TList *fCutGammaArray;
    TList *fCutElectronArray;
    TList *fCutMesonArray;
    TList **fGammasPool;
    AliConversionCuts *fConversionCuts;
    TH1F **hESDConvGammaPt;
    TH1F **hESDDalitzElectronPt;
    TH1F **hESDDalitzPositronPt;
    TH2F **hESDEposEnegPsiPairDPhi;
    TH2F **hESDEposEnegInvMassPt;
    TH2F **hESDEposEnegLikeSignBackInvMassPt;
    TH2F **hESDMotherInvMassPt;
    TH2F **hESDPi0MotherInvMassPt;
    TH2F **hESDPi0MotherDiffInvMassPt;
    TH2F **hESDPi0MotherDiffLimInvMassPt;
    THnSparseF **sESDMotherInvMassPtZM;
    TH2F **hESDMotherBackInvMassPt;
    THnSparseF **sESDMotherBackInvMassPtZM;
    TH1F **hMCAllGammaPt;
    TH1F **hMCConvGammaPt;
    TH1F **hMCConvGammaRSPt;
    TH1F **hMCAllPositronsPt;
    TH1F **hMCAllElectronsPt;
    TH1F **hMCPi0DalitzGammaPt;
    TH1F **hMCPi0DalitzElectronPt;
    TH1F **hMCPi0DalitzPositronPt;
    TH1F **hMCPi0Pt;
    TH1F **hMCPi0GGPt;
    TH1F **hMCEtaPt;
    TH1F **hMCEtaGGPt;
    TH1F **hMCPi0InAccPt;
    TH1F **hMCEtaInAccPt;
    TH1F **hMCChiCPt;
    TH1F **hMCChiCInAccPt;
    TH2F **hESDEposEnegTruePi0DalitzInvMassPt;
    TH2F **hESDEposEnegTrueEtaDalitzInvMassPt;
    TH2F **hESDEposEnegTruePhotonInvMassPt;
    TH2F **hESDEposEnegTrueJPsiInvMassPt;
    TH2F **hESDTrueMotherChiCInvMassPt;
    TH2F **hESDTrueMotherChiCDiffInvMassPt;
    TH2F **hESDTrueMotherInvMassPt;
    TH2F **hESDTrueMotherPi0GGInvMassPt;
    TH2F **hESDTruePrimaryMotherInvMassMCPt;
    TH2F **hESDTruePrimaryPi0DalitzESDPtMCPt;
    TH2F **hESDTrueSecondaryMotherInvMassPt;
    TH2F **hESDTrueSecondaryMotherFromK0sInvMassPt;
    TH2F **hESDTrueBckGGInvMassPt;
    TH2F **hESDTrueBckContInvMassPt;
    TH2F **hESDTrueMotherGGInvMassPt;
    TH1F **hESDTrueConvGammaPt;
    TH1F **hESDTruePositronPt;
    TH1F **hESDTrueElectronPt;
    TH1F **hESDTrueSecConvGammaPt;
    TH1F **hESDTrueSecPositronPt;
    TH1F **hESDTrueSecElectronPt;
    TH1F **hESDTruePi0DalitzConvGammaPt;
    TH1F **hESDTruePi0DalitzPositronPt;
    TH1F **hESDTruePi0DalitzElectronPt;
    TH1F **hESDTruePi0DalitzSecConvGammaPt;
    TH1F **hESDTruePi0DalitzSecPositronPt;
    TH1F **hESDTruePi0DalitzSecElectronPt;
      //if(fDoMesonAnalysis){


    TH1I **hNEvents;
    TH1I **hNGoodESDTracks;
        
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
    Bool_t fIsFromMBHeader;
		Bool_t fIsMC;

	private:
		AliAnalysisTaskGammaConvDalitzV1( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented
		AliAnalysisTaskGammaConvDalitzV1& operator=( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented

		ClassDef( AliAnalysisTaskGammaConvDalitzV1, 3 );
};

#endif // ALIANALYSISTASKGAMMACONVDALITZV1_H

