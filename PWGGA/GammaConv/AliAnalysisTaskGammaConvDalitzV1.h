#ifndef ALIANALYSISTASKGAMMACONVDALITZV1_H
#define ALIANALYSISTASKGAMMACONVDALITZV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliDalitzElectronSelector.h"
#include "AliConvEventCuts.h"
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
			fIsHeavyIon = flag;
		}
		
		void SetIsMC(Bool_t isMC){fIsMC=isMC;}
		void SetEventCutList(Int_t nCuts, TList *CutArray){
			fnCuts= nCuts;
			fCutEventArray = CutArray;
		}
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
		void SetDoChicAnalysis(Bool_t flag){ fDoChicAnalysis = flag; }
		void SetDoMesonQA(Bool_t flag){ fDoMesonQA = flag; }
		void SetProductionVertextoVGamma(Bool_t flag) { fSetProductionVertextoVGamma = flag; }
	
	private:
		void InitBack();
		void ProcessPhotonCandidates();
		void ProcessTruePhotonCandidates(AliAODConversionPhoton*);
		void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
		void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
		void ProcessElectronCandidates();
		void ProcessVirtualGammasCandidates();
		void ProcessMCParticles();
		void CountESDTracks();
		void CalculatePi0DalitzCandidates();
		void CalculateBackground();
		void UpdateEventByEventData();
		void FillElectronQAHistos(AliAODConversionPhoton *Vgamma) const;
		Double_t GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const;
		Bool_t IsDalitz(TParticle *fMCMother) const;
		Bool_t IsPi0DalitzDaughter( Int_t label ) const;
		
		AliV0ReaderV1 							*fV0Reader;
		AliDalitzElectronSelector				*fElecSelector;
		AliGammaConversionAODBGHandler 			**fBGHandler;
		AliESDEvent 							*fESDEvent;
		AliMCEvent 								*fMCEvent;
		AliStack 								*fMCStack;
		TList 									**fCutFolder;
		TList 									**fESDList;
		TList 									**fBackList;
		TList 									**fMotherList;
		TList 									**fTrueList;
		TList 									**fMCList;
		TList 									**fQAFolder;
		TList 									*fOutputContainer;
		TClonesArray 							*fReaderGammas;
		vector<Int_t> 							fSelectorElectronIndex;
		vector<Int_t> 							fSelectorPositronIndex;
		TList 									*fGoodGammas;
		TList 									*fGoodVirtualGammas;
		TList 									*fGoodElectrons;
		TList 									*fGoodPositrons;
		TList 									*fCutEventArray;
		TList 									*fCutGammaArray;
		TList 									*fCutElectronArray;
		TList	 								*fCutMesonArray;
		TList	 								**fGammasPool;
		AliConvEventCuts 						*fEventCuts;
		AliConversionPhotonCuts 				*fConversionCuts;
		TH1F 									**hESDConvGammaPt;
		TH1F 									**hESDConvGammaEta;
		TH2F 									**hESDConvGammaZR;
		TH1F 									**hESDDalitzElectronPt;
		TH1F 									**hESDDalitzPositronPt;
		TH1F 									**hESDDalitzElectronPhi;
		TH1F 									**hESDDalitzPositronPhi;
		TH1F 									**hESDDalitzElectronAfterPt;
		TH1F 									**hESDDalitzPositronAfterPt;
		TH1F 									**hESDDalitzElectronAfterEta;
		TH1F 									**hESDDalitzElectronAfterEtaPCut;
		TH1F 									**hESDDalitzPositronAfterEta;
		TH1F 									**hESDDalitzPositronAfterEtaPCut;
		TH1F 									**hESDDalitzElectronAfterPhi;
		TH1F 									**hESDDalitzPositronAfterPhi;
		TH1F 									**hESDDalitzElectronAfterNClsITS;
		TH1F 									**hESDDalitzPositronAfterNClsITS;
		TH2F 									**hESDDalitzElectronAfterNFindClsTPC;
		TH2F 									**hESDDalitzPositronAfterNFindClsTPC;
		TH2F 									**hESDDalitzElectronAfterNClsTPC;
		TH2F 									**hESDDalitzPositronAfterNClsTPC;
		TH2F 									**hESDDalitzElectronAfterNCrossedRowsTPC;
		TH2F 									**hESDDalitzPositronAfterNCrossedRowsTPC;
		TH2F 									**hESDDalitzPosEleAfterDCAxy;
		TH2F 									**hESDDalitzPosEleAfterDCAz;
		TH2F 									**hESDDalitzElectronAfterTPCdEdxVsP;
		TH2F 									**hESDDalitzPositronAfterTPCdEdxVsP;
		TH2F 									**hESDDalitzElectronAfterTPCdEdxSignalVsP;
		TH2F 									**hESDDalitzPositronAfterTPCdEdxSignalVsP;
		TH2F 									**hESDDalitzElectronAfterTPCdEdxVsEta;
		TH2F 									**hESDDalitzPositronAfterTPCdEdxVsEta;
		TH2F 									**hESDDalitzElectronAfterTPCdEdxVsPhi;
		TH2F 									**hESDDalitzPositronAfterTPCdEdxVsPhi;
		TH1F 									**hESDMotherPhi;
		TH2F 									**hESDEposEnegPsiPairDPhi;
		TH2F 									**hESDEposEnegInvMassPt;
		TH2F 									**hESDEposEnegAfterMassCutInvMassPi0Pt;
		TH2F 									**hESDEposEnegInvMassPi0Pt;
		TH2F 									**hESDEposEnegLikeSignBackInvMassPt;
		TH2F 									**hESDMotherInvMassPt;
		TH2F 									**hESDPi0MotherInvMassPt;
		TH2F 									**hESDPi0MotherDiffInvMassPt;
		TH2F 									**hESDPi0MotherDiffLimInvMassPt;
		THnSparseF 								**sESDMotherInvMassPtZM;
		TH2F 									**hESDMotherBackInvMassPt;
		THnSparseF 								**sESDMotherBackInvMassPtZM;
		TH1F 									**hMCAllGammaPt;
		TH1F 									**hMCConvGammaPt;
		TH1F 									**hMCConvGammaRSPt;
		TH1F 									**hMCAllPositronsPt;
		TH1F 									**hMCAllElectronsPt;
		TH1F 									**hMCConvGammaEta;
		TH1F 									**hMCAllPositronsEta;
		TH1F 									**hMCAllElectronsEta;
		TH1F 									**hMCPi0DalitzGammaPt;
		TH1F 									**hMCPi0DalitzElectronPt;
		TH1F 									**hMCPi0DalitzPositronPt;
		TH1F 									**hMCPi0Pt;
		TH1F 									**hMCPi0GGPt;
		TH1F 									**hMCEtaPt;
		TH1F 									**hMCEtaGGPt;
		TH1F 									**hMCPi0InAccPt;
		TH1F 									**hMCEtaInAccPt;
		TH1F 									**hMCChiCPt;
		TH1F 									**hMCChiCInAccPt;
		TH2F 									**hMCPi0EposEnegInvMassPt;
		TH2F 									**hMCEtaEposEnegInvMassPt;
		TH2F 									**hESDEposEnegTruePi0DalitzInvMassPt;
		TH1F 									**hESDEposEnegTruePrimPi0DalitzInvMass;
		TH2F 									**hESDEposEnegTruePi0DalitzPsiPairDPhi;
		TH2F 									**hESDEposEnegTrueEtaDalitzInvMassPt;
		TH1F 									**hESDEposEnegTruePrimEtaDalitzInvMass;
		TH2F 									**hESDEposEnegTrueEtaDalitzPsiPairDPhi;
		TH2F 									**hESDEposEnegTruePhotonInvMassPt;
		TH2F 									**hESDEposEnegTrueInvMassPt;
		TH2F 									**hESDEposEnegTruePhotonPsiPairDPhi;
		TH2F 									**hESDEposEnegTruePhotonPsiPairDPhiPtCut;
		TH2F 									**hESDEposEnegTrueJPsiInvMassPt;
		TH2F 									**hESDTrueMotherChiCInvMassPt;
		TH2F 									**hESDTrueMotherChiCDiffInvMassPt;
		TH2F 									**hESDTrueMotherInvMassPt;
		TH2F 									**hESDTrueMotherDalitzInvMassPt;
		TH2F 									**hESDTrueMotherPi0GGInvMassPt;
		TH2F 									**hESDTruePrimaryMotherPi0GGInvMassPt;
		TH2F 									**hESDTrueSecondaryMotherPi0GGInvMassPt;
		TH2F 									**hESDTruePrimaryMotherInvMassMCPt;
		TH2F 									**hESDTruePrimaryMotherInvMassPt;
		TH2F 									**hESDTruePrimaryMotherW0WeightingInvMassPt;
		TH2F 									**hESDTruePrimaryPi0DalitzESDPtMCPt;
		TH2F 									**hESDTrueSecondaryMotherInvMassPt;
		TH2F 									**hESDTrueSecondaryMotherFromK0sInvMassPt;
		TH2F 									**hESDTrueBckGGInvMassPt;
		TH2F 									**hESDTrueBckContInvMassPt;
		TH2F 									**hESDTrueMotherGGInvMassPt;
		TH1F 									**hESDTrueConvGammaPt;
		TH1F 									**hESDTruePositronPt;
		TH1F 									**hESDTrueElectronPt;
		TH1F 									**hESDTrueSecConvGammaPt;
		TH1F 									**hESDTrueSecPositronPt;
		TH1F 									**hESDTrueSecElectronPt;
		TH1F 									**hESDTruePi0DalitzConvGammaPt;
		TH1F 									**hESDTruePi0DalitzPositronPt;
		TH1F 									**hESDTruePi0DalitzElectronPt;
		TH1F 									**hESDTruePi0DalitzSecConvGammaPt;
		TH1F 									**hESDTruePi0DalitzSecPositronPt;
		TH1F 									**hESDTruePi0DalitzSecElectronPt;

		TH1I 									**hNEvents;
		TH1I 									**hNGoodESDTracks;
		TH2F 									**hNGoodESDTracksVsNGoodGammas;
		TH2F 									**hNGoodESDTracksVsNGoodVGammas;
		TH1I  									**hNV0Tracks;
		TProfile 								**hEtaShift;
			
		TRandom3 								fRandom;
		Double_t 								*fUnsmearedPx;
		Double_t 								*fUnsmearedPy;
		Double_t 								*fUnsmearedPz;
		Double_t 								*fUnsmearedE;
		Double_t 								*fUnsmearedVPx;
		Double_t 								*fUnsmearedVPy;
		Double_t 								*fUnsmearedVPz;
		Double_t 								*fUnsmearedVE;
		
		Int_t 									fnCuts;
		Int_t 									fiCut;
		Int_t 									fNumberOfESDTracks;
		Int_t 									fNumberOfESDTrackskBoth;
		Int_t 									fNVirtualGammas;
		Bool_t 									fMoveParticleAccordingToVertex;
		Int_t									fIsHeavyIon;
		Bool_t 									fDoMesonAnalysis;
		Bool_t 									fDoChicAnalysis;
		Bool_t 									fDoMesonQA;
		Bool_t 									fSetProductionVertextoVGamma;
		Bool_t 									fIsFromMBHeader;
		Bool_t 									fIsMC;

	private:
		AliAnalysisTaskGammaConvDalitzV1( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented
		AliAnalysisTaskGammaConvDalitzV1& operator=( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented

		ClassDef( AliAnalysisTaskGammaConvDalitzV1, 4 );
};

#endif // ALIANALYSISTASKGAMMACONVDALITZV1_H

