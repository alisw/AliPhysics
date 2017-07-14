
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
#include "AliCaloPhotonCuts.h"
#include "AliGammaConversionAODBGHandler.h"
#include "TProfile2D.h"
#include <vector>

class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
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

        void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
			
		void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
		void SetIsHeavyIon(Int_t flag){
			if (flag == 1 || flag ==2 ){
				fIsHeavyIon = 1;    
			} else {
				fIsHeavyIon = 0;    
			}
		}
		
		void SetIsMC(Bool_t isMC){fIsMC=isMC;}
		void SetEventCutList(Int_t nCuts, TList *CutArray){ 
			fnCuts= nCuts;	
			fEventCutArray = CutArray;
		}
		void SetConversionCutList(TList *CutArray){ fGammaCutArray = CutArray;}
		void SetClusterCutList(TList *CutArray){ fClusterCutArray = CutArray;}
		void SetPionCutList(TList *CutArray){ fPionCutArray = CutArray;}
		void SetNeutralPionCutList(TList *CutArray){ fNeutralPionMesonCutArray = CutArray; }
		void SetMesonCutList(TList *CutArray){ fMesonCutArray = CutArray; }
		void SetDoMesonQA(Bool_t flag){ fDoMesonQA = flag; }
		void SetNeutralPionMode(Int_t mode){fNeutralPionMode = mode; }
        void SetTolerance(Double_t tol){fTolerance=tol;}

	private:

		void InitBack();
		
		// routines for photon selection from conversions
		void ProcessConversionPhotonCandidates();
		void ProcessTrueConversionPhotonCandidates(AliAODConversionPhoton*);
		
		// routines for photon selection from clusters
		void ProcessCaloPhotonCandidates();
		void ProcessTrueCaloPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate);
		
		void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
		void MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);

		void FixPzToMatchPDGInvMassPi0(AliAODConversionMother* particle);

		// routines for neutral pion candidates from pure conversion
		void ProcessNeutralPionCandidatesPureConversions();	
		void ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		void ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		
		// routines for neutral pion candidates from pure calo
		void ProcessNeutralPionCandidatesPureCalo();
		void ProcessTrueNeutralPionCandidatesPureCalo(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

		// routines for neutral pion candidates from mixed conv + calo
		void ProcessNeutralPionCandidatesMixedConvCalo();
		void ProcessTrueNeutralPionCandidatesMixedConvCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		
		void ProcessPionCandidates();
		void ProcessMCParticles();
		void CalculateMesonCandidates();
		void CalculateBackground();
		void UpdateEventByEventData();

		Bool_t IsPiPlPiMiPiZeroDecay(TParticle *fMCMother) const;
		Bool_t IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const;
		Bool_t IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const;
		Bool_t GammaIsNeutralMesonPiPlPiMiPiZeroDaughter( Int_t label ) const;

		Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);

        Bool_t KinematicCut(AliAODConversionMother *negpion, AliAODConversionMother *pospion, AliAODConversionMother *neutpion, AliAODConversionMother *omega);


		AliV0ReaderV1 					*fV0Reader;									// V0Reader for basic conversion photon selection
        TString                         fV0ReaderName;
		AliPrimaryPionSelector			*fPionSelector;								// primary charged pion selector, basic selection of pi+,pi-
		AliGammaConversionAODBGHandler	**fBGHandlerPiPl;							// BG handler Pos Pion
		AliGammaConversionAODBGHandler	**fBGHandlerPiMi;							// BG handler Neg Pion
		AliESDEvent 					*fESDEvent;									// current event
		AliMCEvent 						*fMCEvent;									// current MC event
		TList 							**fCutFolder;								// list of output folders with main cut name 
		TList 							**fESDList;									// list with main output histograms for data
//		TList 							**fBackList;								// list with THnSparseF for BG
//		TList 							**fMotherList;								// list with THnSparseF for FG
		TList 							**fTrueList;								// list with validated reconstructed MC histograms
		TList 							**fMCList;									// list with pure MC histograms
		TList 							*fOutputContainer;							// output container
		TClonesArray 					*fReaderGammas;								// array with photon from fV0Reader
		vector<Int_t> 					fSelectorNegPionIndex;						// array with pion indices for negative pions from fPionSelector
		vector<Int_t> 					fSelectorPosPionIndex;						// array with pion indices for positive pions from fPionSelector
		TList 							*fGoodConvGammas;							// good conv gammas after selection
		TList 							*fClusterCandidates; 						//! good calo gammas after selection 
		TList 							*fNeutralPionCandidates;					// good neutral pion candidates
		TList 							*fPosPionCandidates;						// good positive pion candidates
		TList 							*fNegPionCandidates;						// good negative pion candidates
		TList 							*fGoodVirtualParticles;						// combination of pi+pi- candidates
		TList 							*fEventCutArray;							// array with event cuts
		TList 							*fGammaCutArray;							// array with Conversion Cuts
		TList 							*fClusterCutArray;							// array with Cluster Cuts
		TList 							*fPionCutArray;								// array with charged pion cuts
		TList 							*fNeutralPionMesonCutArray;					// array with neutral pion cuts
		TList 							*fMesonCutArray;							// array with neutral meson cuts
		AliConvEventCuts 				*fEventCuts;								// current event cuts
		AliConversionPhotonCuts 		*fConversionCuts;							// current conversion cuts
		AliCaloPhotonCuts				*fClusterCuts;								// current cluster cuts
		
		// reconstructed particles
		TH1F 							**fHistoConvGammaPt;						// array of histos of conversion photon, pt
		TH1F 							**fHistoConvGammaEta;						// array of histos of conversion photon, eta
		TH1F 							**fHistoClusterGammaPt;						// array of histos of Cluster photon, pt
		TH1F 							**fHistoClusterGammaEta;					// array of histos of Cluster photon, eta
		TH1F 							**fHistoNegPionPt;							// array of histos of negative pion, pt
		TH1F 							**fHistoPosPionPt;							// array of histos of positive pion, pt
		TH1F 							**fHistoNegPionPhi;							// array of histos of negative pion, phi
		TH1F 							**fHistoPosPionPhi;							// array of histos of positive pion, phi
		TH1F 							**fHistoNegPionEta;							// array of histos of negative pion, eta
		TH1F 							**fHistoPosPionEta;							// array of histos of positive pion, eta
		TH2F 							**fHistoNegPionClsTPC;						// array of histos of negative pion, findable tpc cluster, pT
		TH2F 							**fHistoPosPionClsTPC;						// array of histos of positive pion, findable tpc cluster, pT
		TH2F 							**fHistoPionDCAxy;							// array of histos of pion, dca_xy, pT
		TH2F 							**fHistoPionDCAz;							// array of histos of pion, dca_z, pT
		TH2F 							**fHistoPionTPCdEdxNSigma;					// array of histos of pion, p, TPC nSigma dEdx pion
		TH2F 							**fHistoPionTPCdEdx;						// array of histos of pion, p, TPC dEdx
		TH2F 							**fHistoPionPionInvMassPt;					// array of histos of pion pion, invMass, pT_{pi+pi-}
		TH2F	 						**fHistoGammaGammaInvMassPt;				// array of histos of gamma-gamma, invMass, pT_{gamma gamma}
		TH2F	 						**fHistoMotherInvMassPt;					// array of histos of pi+pi-pi0 same event, invMass, pT_{pi+pi-pi0}
        TH2F                            **fHistoMotherInvMassPtRejectedKinematic;   // array of histos of rejected pi+pi-pi0 same event, invMass, pT_{pi+pi-pi0}
//		THnSparseF 						**fTHnSparseMotherInvMassPtZM;				// array of THnSparseF of pi+pi-pi0 same event, invMass, pT_{pi+pi-pi0}, Z, M
		TH2F 							**fHistoMotherSameDiff1Diff2BackInvMassPt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameDiff1Diff1BackInvMassPt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameSameDiff2BackInvMassPt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameDiff1SameBackInvMassPt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}
//		THnSparseF				 		**fTHnSparseMotherBackInvMassPtZM;			// array of THnSparseF of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}, Z, M

        TH2F							**fHistoMotherInvMassSubPi0;							// invariant mass of (pi+,pi-,pi0) - invariant mass of pi0
        TH2F 							**fHistoMotherSameDiff1Diff2BackInvMassSubPi0Pt;	// array of histos of pi+pi-pi0 mixed event, invMass-invMass(pi0), pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameDiff1Diff1BackInvMassSubPi0Pt;	// array of histos of pi+pi-pi0 mixed event, invMass-invMass(pi0), pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameSameDiff2BackInvMassSubPi0Pt;		// array of histos of pi+pi-pi0 mixed event, invMass-invMass(pi0), pT_{pi+pi-pi0}
		TH2F 							**fHistoMotherSameDiff1SameBackInvMassSubPi0Pt;		// array of histos of pi+pi-pi0 mixed event, invMass-invMass(pi0), pT_{pi+pi-pi0}

		TH2F							**fHistoMotherInvMassFixedPzPi0;								// invariant mass of (pi+,pi-,pi0) - invariant mass of pi0
		TH2F 							**fHistoMotherSameDiff1Diff2BackInvMassFixedPzPi0Pt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}, the Pz of the pi0 was fixed such that its invMass matches the PDG value
		TH2F 							**fHistoMotherSameDiff1Diff1BackInvMassFixedPzPi0Pt;	// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}, the Pz of the pi0 was fixed such that its invMass matches the PDG value
		TH2F 							**fHistoMotherSameSameDiff2BackInvMassFixedPzPi0Pt;		// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}, the Pz of the pi0 was fixed such that its invMass matches the PDG value
		TH2F 							**fHistoMotherSameDiff1SameBackInvMassFixedPzPi0Pt;		// array of histos of pi+pi-pi0 mixed event, invMass, pT_{pi+pi-pi0}, the Pz of the pi0 was fixed such that its invMass matches the PDG value


        // angle distributions
        TH2F 							**fHistoAngleOmegaPiPlPiMi;                 // angle between combined Pi+ and Pi- and omega
        TH2F 							**fHistoAngleOmegaPiZero;                   // angle between Pi0 and omega
        TH2F 							**fHistoAngleOmegaPiPl;                     // angle between Pi+ and omega
        TH2F 							**fHistoAngleOmegaPiMi;                     // angle between Pi- and omega
        TH2F 							**fHistoAnglePiPlPiMi;                      // angle between Pi+ and Pi-
        TH2F 							**fHistoAnglePiZeroPiMi;                    // angle between Pi0 and Pi-
        TH2F 							**fHistoAnglePiPlPiZero;                    // angle between Pi+ and Pi0
        TH2F 							**fHistoAngleSum;                           // angle between omega and Pi0 + angle between Pi+ and Pi- + angle between Pi0 and Pi-
        TH2F                            **fHistoTrueAngleSum;

        // pure MC properties
		TH1F 							**fHistoMCAllGammaPt;						// array of histos of all produced gammas in the specified y range
		TH1F 							**fHistoMCConvGammaPt;						// array of histos of all converted gammas in the specified y range 
		TH1F 							**fHistoMCAllPosPionsPt;					// array of histos with all produced primary positive pions in the specified y range
		TH1F 							**fHistoMCAllNegPionsPt;					// array of histos with all produced primary negative pions in the specified y range
		TH1F 							**fHistoMCGammaFromNeutralMesonPt;			// array of histos of all produced gammas from omega or eta via pi+pi-pi0 in the specified y range/
		TH1F 							**fHistoMCPosPionsFromNeutralMesonPt;		// array of histos of all produced positive pions from omega or eta via pi+pi-pi0 in the specified y range/
		TH1F 							**fHistoMCNegPionsFromNeutralMesonPt;		// array of histos of all produced negative pions from omega or eta via pi+pi-pi0 in the specified y range/
		TH1F 							**fHistoMCEtaPiPlPiMiPiZeroPt;				// array of histos of produced etas via pi+pi-pi0 in the specified y range
		TH1F 							**fHistoMCEtaPiPlPiMiPiZeroInAccPt;			// array of histos of produced etas via pi+pi-pi0 in the specified y range, with decay products in respective y, eta ranges 
		TH1F 							**fHistoMCOmegaPiPlPiMiPiZeroPt;			// array of histos of produced omegas via pi+pi-pi0 in the specified y range
		TH1F 							**fHistoMCOmegaPiPlPiMiPiZeroInAccPt;		// array of histos of produced omegas via pi+pi-pi0 in the specified y range, with decay products in respective y, eta ranges 

		// reconstructed particles MC validated
		TH2F 							**fHistoTrueMotherPiPlPiMiPiZeroInvMassPt;	// histos with reconstructed validated eta or omega, inv mass, pT
		TH2F 							**fHistoTrueMotherGammaGammaInvMassPt;		// histos with reconstructed validated pi0, inv mass, pT
		TH2F 							**fHistoTrueMotherGammaGammaFromEtaInvMassPt;	// histos with reconstructed validated pi0, inv mass, pT
		TH2F 							**fHistoTrueMotherGammaGammaFromOmegaInvMassPt;	// histos with reconstructed validated pi0, inv mass, pT
		TH1F 							**fHistoTrueConvGammaPt;					// histos with reconstructed validated conv gamma, pT
		TH1F 							**fHistoTrueConvGammaFromNeutralMesonPt;	// histos with reconstructed validated conv gamma from eta or omega via pi0, pT
		TH1F 							**fHistoTrueClusterGammaPt;					// histos with reconstructed validated cluster gamma, pT
		TH1F 							**fHistoTrueClusterGammaFromNeutralMesonPt;	// histos with reconstructed validated cluster gamma from eta or omega via pi0, pT
		TH1F				 			**fHistoTruePosPionPt;						// histos with reconstructed validated positive pion, pT
		TH1F 							**fHistoTruePosPionFromNeutralMesonPt;		// histos with reconstructed validated positive pion from eta or omega, pT
		TH1F 							**fHistoTrueNegPionPt;						// histos with reconstructed validated negative pion, pT
		TH1F 							**fHistoTrueNegPionFromNeutralMesonPt;		// histos with reconstructed validated negative pion from eta or omega, pT
		TH2F 							**fHistoTruePionPionInvMassPt;				// histos with reconstructed validated two pion, invariant mass, pT
		TH2F 							**fHistoTruePionPionFromSameMotherInvMassPt;// histos with reconstructed validated two pion from same mother, invariant mass, pT
		TH2F 							**fHistoTruePionPionFromEtaInvMassPt;		// histos with reconstructed validated two pion from eta , invariant mass, pT
		TH2F 							**fHistoTruePionPionFromOmegaInvMassPt;		// histos with reconstructed validated two pion from omega, invariant mass, pT
		TH2F				 			**fHistoDoubleCountTruePi0InvMassPt;		//! array of histos with double counted pi0s, invMass, pT
		TH2F				 			**fHistoDoubleCountTrueEtaInvMassPt;		//! array of histos with double counted etas, invMass, pT
		TH2F				 			**fHistoDoubleCountTrueOmegaInvMassPt;		//! array of histos with double counted omegas, invMass, pT
		TH2F				 			**fHistoDoubleCountTrueConvGammaRPt;		//! array of histos with double counted photons, R, pT
		vector<Int_t>					fVectorDoubleCountTruePi0s;					//! vector containing labels of validated pi0
		vector<Int_t>					fVectorDoubleCountTrueEtas;					//! vector containing labels of validated eta
		vector<Int_t>					fVectorDoubleCountTrueOmegas;				//! vector containing labels of validated omega
		vector<Int_t>					fVectorDoubleCountTrueConvGammas;			//! vector containing labels of validated photons
		// Event properties
		TH1I 							**fHistoNEvents;							// histo for event counting
		TH1I 							**fHistoNGoodESDTracks;						// histo number of reconstructed primary tracks
		TProfile 						**fProfileEtaShift;							// profile for eta shift bookkeeping
		TH2F							**fHistoSPDClusterTrackletBackground;		//! array of histos with SPD tracklets vs SPD clusters for background rejection

			
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
		Int_t							fNeutralPionMode;							// Flag how neutral pion is reconstructed 0=PCM-PCM, 1=PCM-Calo, 2=Calo-Calo
        Double_t                        fTolerance;                                 // tolerance in rad for angle cuts

	private:
		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero( const AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& ); // Not implemented
		AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& operator=( const AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero& ); // Not implemented

        ClassDef(AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero, 13);
};

#endif // ALIANALYSISTASKNEUTRALMESONTOPIPLPIMIPIZERO_H

