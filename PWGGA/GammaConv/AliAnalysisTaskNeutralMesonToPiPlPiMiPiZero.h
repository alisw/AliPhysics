#ifndef ALIANALYSISTASKNEUTRALMESONTOPIPLPIMIPIZERO_H
#define ALIANALYSISTASKNEUTRALMESONTOPIPLPIMIPIZERO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliPrimaryPionSelector.h"
#include "AliConversionMesonCuts.h"
#include "AliConvEventCuts.h"
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

class AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero: public AliAnalysisTaskSE
{
	public:

		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero();
		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero( const char* name );
		virtual ~AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero();

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
		void SetEventCutList(Int_t nCuts, TList *CutArray){
			fnCuts= nCuts;
			fEventCutArray = CutArray;
		}

		void SetPionCutList(TList *CutArray){
			fPionCutArray = CutArray;
		}
		void SetNeutralPionCutList(TList *CutArray){
			fNeutralPionMesonCutArray = CutArray;
		}
		void SetMesonCutList(TList *CutArray){
			fMesonCutArray = CutArray;
		}
		void SetDoMesonQA(Bool_t flag){ fDoMesonQA = flag; }
	

	private:

		void InitBack();
		void ProcessPhotonCandidates();
		void ProcessTruePhotonCandidates(AliAODConversionPhoton*);
		void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
		void MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    	void ProcessNeutralPionCandidatesPureConversions();	
		void ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		void ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		void ProcessPionCandidates();
		void ProcessMCParticles();
		void CalculateMesonCandidates();
        void CalculateBackground();
		void UpdateEventByEventData();
        
		Bool_t IsPiPlPiMiPiZeroDecay(TParticle *fMCMother) const;
		Bool_t IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const;
		Bool_t IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const;
		Bool_t GammaIsNeutralMesonPiPlPiMiPiZeroDaughter( Int_t label ) const;

		AliV0ReaderV1 					*fV0Reader;									//
		AliPrimaryPionSelector			*fPionSelector;								//
		AliGammaConversionAODBGHandler 	**fBGHandler;								//
		AliESDEvent 					*fESDEvent;									//
		AliMCEvent 						*fMCEvent;									//
		AliStack 						*fMCStack;									//
		TList 							**fCutFolder;								//
		TList 							**fESDList;									//
		TList 							**fBackList;								//
		TList 							**fMotherList;								//
		TList 							**fTrueList;								//
		TList 							**fMCList;									//
		TList 							*fOutputContainer;							//
		TClonesArray 					*fReaderGammas;								//
		vector<Int_t> 					fSelectorNegPionIndex;						//
		vector<Int_t> 					fSelectorPosPionIndex;						//
		TList 							*fGoodGammas;								//
		TList 							*fNeutralPionCandidates;					//
		TList 							*fGoodVirtualParticles;						//
		TList 							*fEventCutArray;							//
		TList 							*fGammaCutArray;							//
		TList 							*fPionCutArray;								//
		TList 							*fNeutralPionMesonCutArray;					//
		TList 							*fMesonCutArray;							//
		AliConvEventCuts 				*fEventCuts;								//
		AliConversionPhotonCuts 		*fConversionCuts;							//
		
		// reconstructed particles
		TH1F 							**fHistoConvGammaPt;						//
		TH1F 							**fHistoConvGammaEta;						//
		TH1F 							**fHistoNegPionPt;							//
		TH1F 							**fHistoPosPionPt;							//
		TH1F 							**fHistoNegPionPhi;							//
		TH1F 							**fHistoPosPionPhi;							//
		TH1F 							**fHistoNegPionEta;							//
		TH1F 							**fHistoPosPionEta;							//
		TH2F 							**fHistoNegPionClsTPC;						//
		TH2F 							**fHistoPosPionClsTPC;						//
		TH2F 							**fHistoPionDCAxy;							//
		TH2F 							**fHistoPionDCAz;							//
		TH2F 							**fHistoPionTPCdEdxNSigma;					//
		TH2F 							**fHistoPionTPCdEdx;						//
		TH2F 							**fHistoPionPionInvMassPt;					//
		TH2F	 						**fHistoGammaGammaInvMassPt;				//
		TH2F	 						**fHistoMotherInvMassPt;					//
		THnSparseF 						**fTHnSparseMotherInvMassPtZM;				//
		TH2F 							**fHistoMotherBackInvMassPt;				//
		THnSparseF				 		**fTHnSparseMotherBackInvMassPtZM;			//
		
		// pure MC properties
		TH1F 							**fHistoMCAllGammaPt;						//
		TH1F 							**fHistoMCConvGammaPt;						//
		TH1F 							**fHistoMCAllPosPionsPt;					//
		TH1F 							**fHistoMCAllNegPionsPt;					//
		TH1F 							**fHistoMCGammaFromNeutralMesonPt;			//
		TH1F 							**fHistoMCPosPionsFromNeutralMesonPt;		//
		TH1F 							**fHistoMCNegPionsFromNeutralMesonPt;		//
		TH1F 							**fHistoMCEtaPiPlPiMiPiZeroPt;				//
		TH1F 							**fHistoMCEtaPiPlPiMiPiZeroInAccPt;			//
		TH1F 							**fHistoMCOmegaPiPlPiMiPiZeroPt;			//
		TH1F 							**fHistoMCOmegaPiPlPiMiPiZeroInAccPt;		//

		// reconstructed particles MC validated
		TH2F 							**fHistoTrueMotherPiPlPiMiPiZeroInvMassPt;	// histos with reconstructed validated eta or omega, inv mass, pT
		TH2F 							**fHistoTrueMotherGammaGammaInvMassPt;		// histos with reconstructed validated pi0, inv mass, pT
		TH1F 							**fHistoTrueConvGammaPt;					// histos with reconstructed validated gamma, pT
		TH1F 							**fHistoTrueConvGammaFromNeutralMesonPt;	// histos with reconstructed validated gamma from eta or omega via pi0, pT
		TH1F				 			**fHistoTruePosPionPt;						// histos with reconstructed validated positive pion, pT
		TH1F 							**fHistoTruePosPionFromNeutralMesonPt;		// histos with reconstructed validated positive pion from eta or omega, pT
		TH1F 							**fHistoTrueNegPionPt;						// histos with reconstructed validated negative pion, pT
		TH1F 							**fHistoTrueNegPionFromNeutralMesonPt;		// histos with reconstructed validated negative pion from eta or omega, pT
		TH2F 							**fHistoTruePionPionInvMassPt;				// histos with reconstructed validated two pion, invariant mass, pT
		TH2F 							**fHistoTruePionPionFromNeutralMesonInvMassPt;// histos with reconstructed validated two pion from eta or omega, invariant mass, pT
		// Event properties
		TH1I 							**fHistoNEvents;							// histo for event counting
		TH1I 							**fHistoNGoodESDTracks;						// histo number of reconstructed primary tracks
		TProfile 						**fProfileEtaShift;							// profile for eta shift bookkeeping
			
		TRandom3 						fRandom;									// random number
		Int_t 							fnCuts;										// number of cuts to be run in parallel
		Int_t 							fiCut;										// current cut
		Int_t 							fNumberOfESDTracks;							// integer with number of primary tracks in this event
		Bool_t 							fMoveParticleAccordingToVertex;				// Flag to move parice to the vertex
		Int_t 							fIsHeavyIon;								// Flag for collision system 0: pp, 1: PbPb, 2: pPb
		Bool_t 							fDoMesonAnalysis;							// Flag for switching on meson analysis
		Bool_t 							fDoMesonQA;									// Flag for switching on small meson QA
		Bool_t 							fIsFromMBHeader;							// Flag for particle whether it belongs to accepted header
		Bool_t 							fIsMC;										// Flag for MC  

	private:
		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero( const AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& ); // Not implemented
		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& operator=( const AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& ); // Not implemented

		ClassDef( AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero, 1 );
};

#endif // ALIANALYSISTASKNEUTRALMESONTOPIPLPIMIPIZERO_H

