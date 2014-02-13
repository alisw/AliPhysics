#ifndef ALIANALYSISTASKHEAVYMESONTHREEBODY_H
#define ALIANALYSISTASKHEAVYMESONTHREEBODY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliPrimaryPionSelector.h"
#include "AliConversionMesonCuts.h"
#include "AliGammaConversionAODBGHandler.h"
#include "TProfile2D.h"

class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliV0Reader;
class AliTriggerAnalysis;

class AliAnalysisTaskEtaToPiPlPiMiGamma: public AliAnalysisTaskSE
{
	public:

		AliAnalysisTaskEtaToPiPlPiMiGamma();
		AliAnalysisTaskEtaToPiPlPiMiGamma( const char* name );
		virtual ~AliAnalysisTaskEtaToPiPlPiMiGamma();

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
			fGammaCutArray = CutArray;
		}
		void SetPionCutList(TList *CutArray){
			fPionCutArray = CutArray;
		}
		void SetMesonCutList(TList *CutArray){
			fMesonCutArray = CutArray;
		}
		void SetDoMesonQA(Bool_t flag){ fDoMesonQA = flag; }
	

	private:

		void InitBack();
		void ProcessPhotonCandidates();
		void ProcessTruePhotonCandidates(AliAODConversionPhoton*);
		void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
		void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    	void ProcessPionCandidates();
		void ProcessMCParticles();
		void CalculateMesonCandidates();
        void CalculateBackground();
		void UpdateEventByEventData();
        
		Bool_t IsPiPlPiMiGammaDecay(TParticle *fMCMother) const;
		Bool_t IsEtaPiPlPiMiGammaDaughter( Int_t label ) const;
		
		

		AliV0ReaderV1 *fV0Reader;
		AliPrimaryPionSelector* fPionSelector;
		AliGammaConversionAODBGHandler **fBGHandler;
		AliESDEvent *fESDEvent;
		AliMCEvent *fMCEvent;
		AliStack *fMCStack;
		TList **fCutFolder;
		TList **fESDList;
		TList **fBackList;
		TList **fMotherList;
		TList **fTrueList;
		TList **fMCList;
		TList *fOutputContainer;
		TClonesArray *fReaderGammas;
		vector<Int_t> fSelectorNegPionIndex;
		vector<Int_t> fSelectorPosPionIndex;
		TList *fGoodGammas;
		TList *fGoodVirtualParticles;
		TList *fGammaCutArray;
		TList *fPionCutArray;
		TList *fMesonCutArray;
		AliConversionCuts *fConversionCuts;
		
		// reconstructed particles
		TH1F **fHistoConvGammaPt;
		TH1F **fHistoConvGammaEta;
		TH1F **fHistoNegPionPt;
		TH1F **fHistoPosPionPt;
		TH1F **fHistoNegPionPhi;
		TH1F **fHistoPosPionPhi;
		TH1F **fHistoNegPionEta;
		TH1F **fHistoPosPionEta;
		TH2F **fHistoNegPionClsTPC;
		TH2F **fHistoPosPionClsTPC;
		TH2F **fHistoPionDCAxy;
		TH2F **fHistoPionDCAz;
		TH2F **fHistoPionTPCdEdxNSigma;
		TH2F **fHistoPionTPCdEdx;
		TH2F **fHistoPionPionInvMassPt;
		TH2F **fHistoMotherInvMassPt;
		THnSparseF **fTHnSparseMotherInvMassPtZM;
		TH2F **fHistoMotherBackInvMassPt;
		THnSparseF **fTHnSparseMotherBackInvMassPtZM;
		
		// pure MC properties
		TH1F **fHistoMCAllGammaPt;
		TH1F **fHistoMCConvGammaPt;
		TH1F **fHistoMCAllPosPionsPt;
		TH1F **fHistoMCAllNegPionsPt;
		TH1F **fHistoMCGammaFromEtaPt;
		TH1F **fHistoMCPosPionsFromEtaPt;
		TH1F **fHistoMCNegPionsFromEtaPt;
		TH1F **fHistoMCEtaPiPlPiMiGammaPt;
		TH1F **fHistoMCEtaGGPt;
		TH1F **fHistoMCEtaDalitzPt;
		TH1F **fHistoMCEtaPiPlPiMiGammaInAccPt;
		
		// reconstructed particles MC validated
		TH2F **fHistoTrueMotherPiPlPiMiGammaInvMassPt;
		TH2F **fHistoTrueMotherGammaGammaInvMassPt;
		TH2F **fHistoTrueMotherDalitzInvMassPt;
		TH1F **fHistoTrueConvGammaPt;
		TH1F **fHistoTrueConvGammaFromEtaPt;
		TH1F **fHistoTruePosPionPt;
		TH1F **fHistoTruePosPionFromEtaPt;
		TH1F **fHistoTrueNegPionPt;
		TH1F **fHistoTrueNegPionFromEtaPt;
		TH2F **fHistoTruePionPionInvMassPt;
		TH2F **fHistoTruePionPionFromEtaInvMassPt;
		// Event properties
		TH1I **fHistoNEvents;
		TH1I **fHistoNGoodESDTracks;
		TProfile **fProfileEtaShift;
			
		TRandom3 fRandom;
		Int_t fnCuts;
		Int_t fiCut;
		Int_t fNumberOfESDTracks;
		Bool_t fMoveParticleAccordingToVertex;
		Bool_t fIsHeavyIon;
		Bool_t fDoMesonAnalysis;
		Bool_t fDoMesonQA;
		Bool_t fIsFromMBHeader;
		Bool_t fIsMC;

	private:
		AliAnalysisTaskEtaToPiPlPiMiGamma( const AliAnalysisTaskEtaToPiPlPiMiGamma& ); // Not implemented
		AliAnalysisTaskEtaToPiPlPiMiGamma& operator=( const AliAnalysisTaskEtaToPiPlPiMiGamma& ); // Not implemented

		ClassDef( AliAnalysisTaskEtaToPiPlPiMiGamma, 1 );
};

#endif // ALIANALYSISTASKHEAVYMESONTHREEBODY_H

