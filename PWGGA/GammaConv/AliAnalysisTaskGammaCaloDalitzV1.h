#ifndef ALIANALYSISTASKGAMMACALODALITZV1_cxx 
#define ALIANALYSISTASKGAMMACALODALITZV1_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliDalitzElectronSelector.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliDalitzElectronCuts.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include <vector>

class AliAnalysisTaskGammaCaloDalitzV1 : public AliAnalysisTaskSE {
	public:

		AliAnalysisTaskGammaCaloDalitzV1();
		AliAnalysisTaskGammaCaloDalitzV1(const char *name);
		virtual ~AliAnalysisTaskGammaCaloDalitzV1();

		virtual void   UserCreateOutputObjects();
		virtual Bool_t Notify();
		virtual void   UserExec(Option_t *);
		virtual void   Terminate(const Option_t*);
		void InitBack();

        void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
		void SetIsHeavyIon(Int_t flag){
			fIsHeavyIon = flag;    
		}

		// base functions for selecting photon and meson candidates in reconstructed data
		void ProcessClusters();
		void ProcessPhotonCandidates();
		void ProcessElectronCandidates();
		void CalculatePi0DalitzCandidates();
		
		// MC functions
		void SetIsMC(Bool_t isMC){fIsMC=isMC;}
		void ProcessMCParticles();
		void ProcessAODMCParticles();
		void RelabelAODPhotonCandidates(Bool_t mode);
		void ProcessTruePhotonCandidates( AliAODConversionPhoton* TruePhotonCandidate);
		void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate);
		void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate);
		void ProcessTruePhotonCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate);
		void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched);
		void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched);
		void ProcessConversionPhotonsForMissingTags();
		void ProcessConversionPhotonsForMissingTagsAOD();
		
		// switches for additional analysis streams or outputs
		void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
		void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
		void SetDoPhotonQA(Int_t flag){fDoPhotonQA = flag;}
		void SetDoClusterQA(Int_t flag){fDoClusterQA = flag;}
		void SetUseTHnSparse(Bool_t flag){fDoTHnSparse = flag;}
		
	    // Setting the cut lists for the conversion photons
		void SetEventCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fEventCutArray = CutArray;
		}

		// Setting the cut lists for the conversion photons
		void SetConversionCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fGammaCutArray = CutArray;
		}

	    // Setting the cut lists for the calo photons
		void SetCaloCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fClusterCutArray = CutArray;
		}
		
		// Setting the cut lists for the meson
		void SetMesonCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fMesonCutArray = CutArray;
		}
	      
		void SetElectronCutList(Int_t nCuts, TList *CutArray){
		  
			fnCuts = nCuts;
			fElectronCutArray = CutArray;
		}
		
		// BG HandlerSettings
		void CalculateBackground();
		//void CalculateBackgroundRP();
		void RotateParticle(AliAODConversionPhoton *gamma);
		void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
		void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
		void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
		void UpdateEventByEventData();
		
		// Additional functions for convenience
		void SetLogBinningXTH2(TH2* histoRebin);
		Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);
		Bool_t CheckIfContainedInString(TString input, Int_t tobechecked);
		Bool_t CheckIfContainedInStringAndAppend(TString &input, Int_t tobechecked);
		Bool_t IsPi0DalitzDaughter( Int_t label ) const;
		Double_t GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const;
		Bool_t IsDalitz(TParticle *fMCMother) const;
		Int_t FindMotherOfPhoton(Int_t particleLabel );

		Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
	
	protected:
		AliV0ReaderV1 						*fV0Reader;							// basic photon Selection Task
        TString                             fV0ReaderName;
		AliDalitzElectronSelector				*fElecSelector;							// basic electron Selection			
		AliGammaConversionAODBGHandler 		**fBGClusHandler;					// BG handler for Cluster
		AliVEvent 							*fInputEvent;						// current event
		AliMCEvent 							*fMCEvent;							// corresponding MC event
		TList 								**fCutFolder;						// Array of lists for containers belonging to cut
		TList 								**fESDList;							// Array of lists with histograms with reconstructed properties   
		TList 								**fBackList;						// Array of lists with BG THnSparseF
		TList 								**fMotherList;						// Array of lists with Signal THnSparseF
		TList 								**fPhotonDCAList;					// Array of lists with photon dca trees
		TList 								**fTrueList;						// Array of lists with histograms with MC validated reconstructed properties
		TList 								**fMCList;							// Array of lists with histograms with pure MC information
		TList 								**fHeaderNameList;					// Array of lists with header names for MC header selection
		TList								**fClusterOutputList;     			//!Array of lists of output histograms for cluster photons
		TList 								*fOutputContainer;					// Output container
		TClonesArray 						*fReaderGammas;						// Array with conversion photons selected by V0Reader Cut
		vector<Int_t> 							fSelectorElectronIndex;                         // Array of electron candidates selected by Electron Selector
		vector<Int_t> 							fSelectorPositronIndex;                         // Array of positron candidates selected by Positron Selector 
		TList 								*fGammaCandidates;					// current list of photon candidates
		TList 								*fClusterCandidates; 				//! current list of cluster candidates
		TList                                                           *fVirtualGammaCandidates;                          //current list of virtual photon candidates
		TList 								*fEventCutArray;					// List with Event Cuts
		AliConvEventCuts 					*fEventCuts;						// EventCutObject
		TList 								*fGammaCutArray;							// List with Conversion Cuts
		TList                                                   *fElectronCutArray;                                     // List with electron cuts
		AliConversionPhotonCuts 			*fConversionCuts;					// ConversionCutObject
		TList 								*fClusterCutArray;					// List with Cluster Cuts
		AliCaloPhotonCuts 					*fCaloPhotonCuts;					// CaloPhotonCutObject
		TList 								*fMesonCutArray;					// List with Meson Cuts
		AliConversionMesonCuts 				*fMesonCuts;						// MesonCutObject
		
		//histograms for Conversions reconstructed quantities
		TH1F 								**fHistoConvGammaPt;				//! histogram conversion photon pT
		TH1F 								**fHistoConvGammaR;				//! histogram conversion photon R
		TH1F 								**fHistoConvGammaEta;				//! histogram conversion photon Eta
		TH1F								**fHistoDalitzElectronPt;			//! histogram dalitz electron candidate Pt
		TH1F								**fHistoDalitzPositronPt;			//! histogram dalitz positron candidate Pt
		TH1F								**fHistoDalitzElectronPhi;			//! histogram dalitz electron candidate Phi
		TH1F								**fHistoDalitzPositronPhi;			//! histogram dalitz positron candidate Phi
		Float_t 							fPtGamma;							//! pt of conversion for tree
		Float_t 							fDCAzPhoton;						//! dcaz of conversion for tree
		Float_t 							fRConvPhoton;						//! R of conversion for tree
		Float_t 							fEtaPhoton;							//! eta of conversion for tree
		UChar_t 							fCharCatPhoton;						//! category of conversion for tree
		UChar_t 							fCharPhotonMCInfo; 					//! MC info of conversion for tree
											// 0: garbage,
											// 1: background
											// 2: secondary photon not from eta or k0s,
											// 3: secondary photon from eta, 
											// 4: secondary photon from k0s, 
											// 5: dalitz
											// 6: primary gamma
		//histograms for mesons reconstructed quantities
		TH2F 								**fHistoMotherInvMassPt;			//! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
		TH1F								**fHistoMotherInvMassOpeningAngleGammaElectron; //! array of histogram with opening angle between gamma and electron
		TH2F 								**fHistoMotherMatchedInvMassPt;		//! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
		THnSparseF 							**fSparseMotherInvMassPtZM;			//! array of THnSparseF with signal + BG for same event photon pairs, inv Mass, pt
		TH2F 								**fHistoMotherBackInvMassPt;		//! array of histogram with BG for mixed event photon pairs, inv Mass, pt
		THnSparseF 							**fSparseMotherBackInvMassPtZM;		//! array of THnSparseF with BG for same event photon pairs, inv Mass, pt
		TH2F 								**fHistoMotherInvMassEalpha;		//! array of histograms with alpha cut of 0.1 for inv mass vs pt
		TH2F 								**fHistoMotherPi0PtY;				//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, Y
		TH2F 								**fHistoMotherEtaPtY;				//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, Y
		TH2F 								**fHistoMotherPi0PtAlpha;			//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, alpha
		TH2F 								**fHistoMotherEtaPtAlpha;			//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, alpha
		TH2F 								**fHistoMotherPi0PtOpenAngle;		//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, openAngle
		TH2F 								**fHistoMotherEtaPtOpenAngle;		//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, openAngle
		TH2F                                				**fHistoMotherPi0ConvPhotonEtaPhi;  //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17 ,eta/phi of conversion photon
		TH2F                                				**fHistoMotherEtaConvPhotonEtaPhi;  //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65 ,eta/phi of conversion photon
		TH2F								**fHistoMotherInvMassECalib;		//! array of histogram with signal + BG for same event photon pairs, inv Mass, energy of cluster
		TH2F								**fHistoMotherInvMassECalibalpha;	//! array of histogram with signal + BG for same event photon pairs, inv Mass, energy of cluster, alpha cut 0.1

		// histograms for rec photons tagged by Calo
		// histograms for rec photon clusters
		TH1F								** fHistoClusGammaPt;				//! array of histos with cluster, pt
		TH1F								** fHistoClusOverlapHeadersGammaPt;	//! array of histos with cluster, pt overlapping with other headers
										
		//histograms for pure MC quantities
		TH1I 								**fHistoMCHeaders;				//! array of histos for header names
		TH1F 								**fHistoMCAllGammaPt;				//! array of histos with all gamma, pT
		TH1F								**fHistoMCAllGammaPi0Pt;                        //! array of histos with all gamma from pi0->dalitz, pT
		TH1F 								**fHistoMCAllGammaEMCALAccPt;			//! array of histos with all gamma in EMCAL acceptance, pT
		TH1F 								**fHistoMCDecayGammaPi0Pt;			//! array of histos with decay gamma from pi0, pT
		TH1F 								**fHistoMCDecayGammaRhoPt;			//! array of histos with decay gamma from rho, pT
		TH1F 								**fHistoMCDecayGammaEtaPt;			//! array of histos with decay gamma from eta, pT
		TH1F 								**fHistoMCDecayGammaOmegaPt;			//! array of histos with decay gamma from omega, pT
		TH1F 								**fHistoMCDecayGammaEtapPt;			//! array of histos with decay gamma from eta', pT
		TH1F 								**fHistoMCDecayGammaPhiPt;			//! array of histos with decay gamma from phi, pT
		TH1F 								**fHistoMCDecayGammaSigmaPt;			//! array of histos with decay gamma from Sigma0, pT
		TH1F 								**fHistoMCConvGammaPt;				//! array of histos with converted gamma, pT
		TH1F 								**fHistoMCConvGammaR;				//! array of histos with converted gamma, R
		TH1F 								**fHistoMCConvGammaEta;				//! array of histos with converted gamma, Eta
		TH1F								**fHistoMCAllPositronsPt;                       //! array of histos with positrons, pT
		TH1F								**fHistoMCDecayPositronPi0Pt;                   //! array of histos with positrons from Pi0->Dalitz, pT
		TH1F								**fHistoMCAllElectronsPt;			//! array of histos with electrons, pT
		TH1F								**fHistoMCDecayElectronPi0Pt;                   //! arrau of histos with positrons form Pi0->Dalitz, pT
		TH1F								**fHistoMCDecayNoPrimElectronPi0DalitzR;	//! array of histos with no prim electron from pi0->Dalitz, radius
		TH1F								**fHistoMCDecayNoPrimPositronPi0DalitzR;	//! array of histos with no prim positron from pi0->Dalitz, ID
		TH1F								**fHistoMCDecayNoPrimElectronPi0DalitzID;	//! array of histos with no prim positron from pi0->Dalitz, radius
		TH1F								**fHistoMCDecayNoPrimPositronPi0DalitzID;	//! array of histos with no prim electron from pi0->Dalitz, ID
		TH1F 								**fHistoMCPi0GGPt;				//! array of histos with weighted pi0, pT
		TH1F 								**fHistoMCPi0GGWOWeightPt;			//! array of histos with unweighted pi0, pT
		TH1F 								**fHistoMCPi0Pt;				//! array of histos with weighted pi0, pT
		TH1F 								**fHistoMCPi0WOWeightPt;			//! array of histos with unweighted pi0, pT
		TH1F								**fHistoMCEtaGGPt;                              //! array of histos with weighted eta->GG, pT
		TH1F                                                            **fHistoMCEtaGGWOWeightPt;                      //! array of histos with unweighted eta->GG, pT
		TH1F 								**fHistoMCEtaPt;				//! array of histos with weighted eta, pT
		TH1F 								**fHistoMCEtaWOWeightPt;			//! array of histos with unweighted eta, pT
		TH1F 								**fHistoMCPi0InAccPt;				//! array of histos with weighted pi0 in acceptance, pT
		TH1F								**fHistoMCPi0InAccOpeningAngleGammaElectron;    //! array of histos with the opening angle between gamma and positron/electron
		TH1F 								**fHistoMCEtaInAccPt;				//! array of histos with weighted eta in acceptance, pT
		TH2F 								**fHistoMCPi0PtY;				//! array of histos with weighted pi0, pT, Y
		TH2F 								**fHistoMCEtaPtY;				//! array of histos with weighted eta, pT, Y
		TH2F 								**fHistoMCPi0PtAlpha;				//! array of histos with weighted pi0, pT, alpha
		TH2F 								**fHistoMCEtaPtAlpha;				//! array of histos with weighted eta, pT, alpha
		TH1F 								**fHistoMCK0sPt;				//! array of histos with weighted K0s, pT
		TH1F 								**fHistoMCK0sWOWeightPt;			//! array of histos with unweighted K0s, pT
		TH2F 								**fHistoMCK0sPtY;				//! array of histos with weighted K0s, pT, Y
		// MC validated reconstructed quantities mesons
		TH2F 								**fHistoTruePi0InvMassPt;						//! array of histos with validated pi0, invMass, pt
		TH2F 								**fHistoTrueEtaInvMassPt;						//! array of histos with validated eta, invMass, pt
		TH2F 								**fHistoTruePi0ShowerInvMassPt;						//! array of histos with validated pi0, invMass, pt
		TH2F 								**fHistoTrueEtaShowerInvMassPt;						//! array of histos with validated eta, invMass, pt
		TH2F 								**fHistoTruePi0NoShowerInvMassPt;						//! array of histos with validated pi0, invMass, pt
		TH2F 								**fHistoTrueEtaNoShowerInvMassPt;						//! array of histos with validated eta, invMass, pt
		TH1F								**fHistoTruePi0OpeningAngleGammaElectron;                                       //! array of histos with validated pi0, opening angle
		TH2F 								**fHistoTruePi0GGInvMassPt;					//! array of histos with validated pi0 ->GG, invMass, pt
		TH2F 								**fHistoTrueEtaGGInvMassPt;						//! array of histos with validated eta ->GG, invMass, pt
		TH2F 								**fHistoTruePi0CaloPhotonInvMassPt;				//! array of histos with validated pi0, photon leading, invMass, pt
		TH2F 								**fHistoTrueEtaCaloPhotonInvMassPt;				//! array of histos with validated eta, photon leading, invMass, pt
		TH2F 								**fHistoTruePi0CaloConvertedPhotonInvMassPt;	//! array of histos with validated pi0, converted photon leading, invMass, pt
		TH2F 								**fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt;	//! array of histos with validated pi0 matched with conv photon, converted photon leading, invMass, pt
		TH2F 								**fHistoTrueEtaCaloConvertedPhotonInvMassPt;	//! array of histos with validated eta, converted photon leading, invMass, pt
		TH2F 								**fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt;	//! array of histos with validated eta matched with conv photon, converted photon leading, invMass, pt
		TH2F 								**fHistoTruePi0CaloElectronInvMassPt;			//! array of histos with validated mothers, electron leading, invMass, pt
		TH2F 								**fHistoTrueEtaCaloElectronInvMassPt;			//! array of histos with validated mothers, electron leading, invMass, pt
		TH2F 								**fHistoTruePi0CaloMergedClusterInvMassPt;		//! array of histos with validated mothers, merged cluster invMass, pt
		TH2F 								**fHistoTrueEtaCaloMergedClusterInvMassPt;		//! array of histos with validated mothers, merged cluster invMass, pt
		TH2F 								**fHistoTrueMotherCaloEMNonLeadingInvMassPt;	//! array of histos with validated mothers, EM non leading, invMass, pt
		TH2F 								**fHistoTruePi0CaloMergedClusterPartConvInvMassPt; //! array of histos with validated mothers, merged cluster part conv, invMass, pt
		TH2F 								**fHistoTrueEtaCaloMergedClusterPartConvInvMassPt; //! array of histos with validated mothers, merged cluster part conv, invMass, pt
		TH2F 								**fHistoTruePrimaryPi0InvMassPt;				//! array of histos with validated weighted primary mothers, invMass, pt
		TH2F								**fHistoTruePrimaryPi0GGInvMassPt;				//! array of histos with validated weighted primary mothers, invMass, pt pi0->GG
		TH2F 								**fHistoTruePrimaryEtaInvMassPt;				//! array of histos with validated weighted primary mothers, invMass, pt
		TH2F								**fHistoTruePrimaryEtaGGInvMassPt;				//! array of histos with validated wwighted primary mothers, invMass, pt pi0->GG
		TH2F 								**fHistoTruePrimaryPi0W0WeightingInvMassPt;		//! array of histos with validated unweighted primary mothers, invMass, pt
		TH2F 								**fHistoTruePrimaryEtaW0WeightingInvMassPt;		//! array of histos with validated unweighted primary mothers, invMass, pt
		TProfile2D 							**fProfileTruePrimaryPi0WeightsInvMassPt;		//! array of profiles with weights for validated primary mothers, invMass, pt	
		TProfile2D 							**fProfileTruePrimaryEtaWeightsInvMassPt;		//! array of profiles with weights for validated primary mothers, invMass, pt	
		TH2F 								**fHistoTruePrimaryPi0MCPtResolPt;				//! array of histos with validated weighted primary pi0, MCpt, resol pt
		TH2F	 							**fHistoTruePrimaryEtaMCPtResolPt;				//! array of histos with validated weighted primary eta, MCpt, resol pt
        TH2F                                **fHistoTrueMotherPi0ConvPhotonEtaPhi;  		//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17 ,eta/phi of conversion photon
        TH2F                                **fHistoTrueMotherEtaConvPhotonEtaPhi;  		//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65 ,eta/phi of conversion photon
		TH2F 								**fHistoTrueSecondaryPi0InvMassPt;				//! array of histos with validated secondary mothers, invMass, pt
		TH2F								**fHistoTrueSecondaryPi0GGInvMassPt;				//! array of histos with validated secondary mothers, invMass, pt pi0->GG
		TH2F 								**fHistoTrueSecondaryPi0FromK0sInvMassPt;		//! array of histos with validated secondary mothers from K0s, invMass, pt
		TH1F 								**fHistoTrueK0sWithPi0DaughterMCPt;				//! array of histos with K0s with reconstructed pi0 as daughter, pt
		TH2F 								**fHistoTrueSecondaryPi0FromEtaInvMassPt;		//! array of histos with validated secondary mothers from eta, invMass, pt
		TH1F 								**fHistoTrueEtaWithPi0DaughterMCPt;				//! array of histos with eta with reconstructed pi0 as daughter, pt
		TH2F 								**fHistoTrueSecondaryPi0FromLambdaInvMassPt;	//! array of histos with validated secondary mothers from Lambda, invMass, pt
		TH1F 								**fHistoTrueLambdaWithPi0DaughterMCPt;			//! array of histos with lambda with reconstructed pi0 as daughter, pt
		TH2F 								**fHistoTrueBckGGInvMassPt;						//! array of histos with pure gamma gamma combinatorial BG, invMass, pt
		TH2F 								**fHistoTrueBckContInvMassPt;					//! array of histos with contamination BG, invMass, pt
		TH2F 								**fHistoTruePi0PtY;								//! array of histos with validated pi0, pt, Y
		TH2F 								**fHistoTrueEtaPtY;								//! array of histos with validated eta, pt, Y
		TH2F 								**fHistoTruePi0PtAlpha;							//! array of histos with validated pi0, pt, alpha
		TH2F 								**fHistoTrueEtaPtAlpha;							//! array of histos with validated eta, pt, alpha
		TH2F 								**fHistoTruePi0PtOpenAngle;						//! array of histos with validated pi0, pt, openAngle
		TH2F 								**fHistoTrueEtaPtOpenAngle;						//! array of histos with validated eta, pt, openAngle
		// MC validated reconstructed quantities photons
		TH1F 								**fHistoTrueConvGammaPt;						//! array of histos with validated conversion photon, pt
		TH1F 								**fHistoTrueConvPi0GammaPt;						//! array of histos with validated conversion photon from pi0, pt
		TH1F 								**fHistoTrueConvGammaEta;						//! array of histos with validated conversion photon, eta
		TH1F								**fHistoTruePositronPt;							//! array of histos with validated positron, pT
		TH1F								**fHistoTrueElectronPt;      						//! array of histos with validated electron, pT
		TH1F								**fHistoTrueSecPositronPt;						//! array of histos with validated sec electron, pT
		TH1F								**fHistoTrueSecElectronPt;                                              //! array of histos with validated sec positron, pT
		TH1F								**fHistoTruePi0DalitzPositronPt;					//! array of histos with validated positron from pi0 Dalitz, pT
		TH1F								**fHistoTruePi0DalitzElectronPt;					//! array of histos with validated electron from pi0 Dalitz, pT
		TH1F								**fHistoTruePi0DalitzSecPositronPt;					//! array of histos with validated sec positron from pi0 Dalitz, pT
		TH1F								**fHistoTruePi0DalitzSecElectronPt;					//! array of histos with validated sec electron from pi0 Dalitz, pT
		TH1F 								**fHistoTruePrimaryConvGammaPt;					//! array of histos with validated primary conversion photon, pt  
		TH2F	 							**fHistoTruePrimaryConvGammaESDPtMCPt;			//! array of histos with validated primary conversion photon, rec pt, mc pt  
		TH1F 								**fHistoTrueSecondaryConvGammaPt;				//! array of histos with validated secondary conversion photon, pt  
		TH1F 								**fHistoTrueSecondaryConvGammaFromXFromK0sPt;	//! array of histos with validated secondary conversion photon from K0s, pt  
		TH1F 								**fHistoTrueSecondaryConvGammaFromXFromLambdaPt;//! array of histos with validated secondary conversion photon from Lambda, pt  
		TH1F								**fHistoTrueClusGammaPt;						//! array of histos with validated cluster (electron or photon), pt
		TH1F								**fHistoTrueClusUnConvGammaPt;					//! array of histos with validated unconverted photon, pt
		TH1F								**fHistoTrueClusUnConvGammaMCPt;					//! array of histos with validated unconverted photon, pt
		TH1F								**fHistoTrueClusElectronPt;						//! array of histos with validated electron, pt
		TH1F								**fHistoTrueClusConvGammaPt;					//! array of histos with validated converted photon, pt
		TH1F								**fHistoTrueClusConvGammaMCPt;					//! array of histos with validated converted photon, pt
		TH1F								**fHistoTrueClusConvGammaFullyPt;				//! array of histos with validated converted photon, fully contained, pt
		TH1F								**fHistoTrueClusMergedGammaPt;					//! array of histos with validated merged photons, electrons, dalitz, pt
		TH1F								**fHistoTrueClusMergedPartConvGammaPt;			//! array of histos with validated merged partially converted photons, pt
		TH1F								**fHistoTrueClusDalitzPt;						//! array of histos with validated Dalitz decay, pt
		TH1F								**fHistoTrueClusDalitzMergedPt;					//! array of histos with validated Dalitz decay, more than one decay product in cluster, pt
		TH1F								**fHistoTrueClusPhotonFromElecMotherPt;			//! array of histos with validated photon from electron, pt
		TH1F								**fHistoTrueClusShowerPt;						//! array of histos with validated shower, pt
        TH1F                                **fHistoTrueClusSubLeadingPt;                   //! array of histos with pi0/eta/eta_prime in subleading contribution
        TH1I                                **fHistoTrueClusNParticles;                     //! array of histos with number of different particles (pi0/eta/eta_prime) contributing to cluster
		TH1F								**fHistoTrueClusEMNonLeadingPt;					//! array of histos with cluster with largest energy by hadron
		TH1F								**fHistoTrueNLabelsInClus;						//! array of histos with number of labels in cluster 
		TH1F								**fHistoTruePrimaryClusGammaPt;					//! array of histos with validated primary cluster, pt
		TH2F								**fHistoTruePrimaryClusGammaESDPtMCPt;			//! array of histos with validated primary cluster, rec Pt, MC pt
		TH1F								**fHistoTruePi0DalitzClusGammaPt;			//! array of histos with validate primary cluster from pi0->dalitz, rec Pt
		TH1F								**fHistoTruePi0DalitzAllClusGammaPt;                    //! array
		TH1F								**fHistoTruePi0DalitzClusGammaMCPt;			//! array of histos with validate primary cluster from pi0->dalitz, rec Pt
		TH2F				 				**fHistoTruePrimaryPi0PhotonPairPtconv;			//! array of histos with validated primary pi0's vs conversion photon pT
		TH1F				 				**fHistoTruePrimaryPi0DCPtconv;					//! array of histos with validated primary pi0's vs conversion photon pT, double counting
		TH1F				 				**fHistoTruePrimaryPi0MissingPtconv;			//! array of histos with validated primary pi0's vs conversion photon pT, missing
		TH2F				 				**fHistoTruePrimaryEtaPhotonPairPtconv;			//! array of histos with validated primary eta's vs conversion photon pT
		TH1F				 				**fHistoTruePrimaryEtaDCPtconv;					//! array of histos with validated primary eta's vs conversion photon pT, double counting
		TH1F				 				**fHistoTruePrimaryEtaMissingPtconv;			//! array of histos with validated primary eta's vs conversion photon pT, missing
		TH2F				 				**fHistoTrueSecondaryPi0PhotonPairPtconv;		//! array of histos with validated secondary pi0's vs conversion photon pT
		TH1F				 				**fHistoTrueSecondaryPi0DCPtconv;				//! array of histos with validated secondary pi0's vs conversion photon pT, double counting
		TH1F				 				**fHistoTrueSecondaryPi0MissingPtconv;			//! array of histos with validated secondary pi0's vs conversion photon pT, missing
		TString 							*fStringRecTruePi0s;							//! array of strings containing the stack position of the reconstructed validated pi0
		TString 							*fStringRecTrueEtas;							//! array of strings containing the stack position of the reconstructed validated eta
		TH2F				 				**fHistoDoubleCountTruePi0InvMassPt;			//! array of histos with double counted pi0s, invMass, pT
		TH2F				 				**fHistoDoubleCountTrueEtaInvMassPt;			//! array of histos with double counted etas, invMass, pT
		TH2F				 				**fHistoDoubleCountTrueConvGammaRPt;				//! array of histos with double counted photons, R, pT
		vector<Int_t>						fVectorDoubleCountTruePi0s;						//! vector containing labels of validated pi0
		vector<Int_t>						fVectorDoubleCountTrueEtas;						//! vector containing labels of validated eta
		vector<Int_t>						fVectorDoubleCountTrueConvGammas;					//! vector containing labels of validated photons
		// event histograms
		TH1I 								**fHistoNEvents;								//! array of histos with event information
		TH1I 								**fHistoNGoodESDTracks;							//! array of histos with number of good tracks (2010 Standard track cuts)
		TH1I 								**fHistoNGammaCandidates;						//! array of histos with number of gamma candidates per event
		TH2F 								**fHistoNGoodESDTracksVsNGammaCanditates;		//! array of histos with number of good tracks vs gamma candidates
		TH2F								**fHistoSPDClusterTrackletBackground;			//! array of histos with SPD tracklets vs SPD clusters for background rejection
		TH1I 								**fHistoNV0Tracks;								//! array of histos with V0 counts
		TProfile 							**fProfileEtaShift;								//! array of profiles with eta shift
		
		// additional variables
		Double_t 							fEventPlaneAngle; 					// EventPlaneAngle
		TRandom3 							fRandom;							// random 
		Int_t 								fNGammaCandidates;					// number of gamma candidates in event
		Double_t 							*fUnsmearedPx; 						//[fNGammaCandidates]
		Double_t 							*fUnsmearedPy; 						//[fNGammaCandidates]
		Double_t 							*fUnsmearedPz; 						//[fNGammaCandidates]
		Double_t 							*fUnsmearedE;  						//[fNGammaCandidates]
        Int_t 								*fMCEventPos;     					//[fNGammaCandidates]
        Int_t 								*fMCEventNeg;     					//[fNGammaCandidates]
		Int_t 								*fESDArrayPos;    					//[fNGammaCandidates]
		Int_t 								*fESDArrayNeg;   					//[fNGammaCandidates]
		Int_t 								fnCuts;								// number of cuts to be analysed in parallel
		Int_t 								fiCut;								// current cut	
		Bool_t 								fMoveParticleAccordingToVertex;		// boolean for BG calculation
		Int_t	 							fIsHeavyIon;						// switch for pp = 0, PbPb = 1, pPb = 2
		Bool_t 								fDoMesonAnalysis;					// flag for meson analysis
		Int_t 								fDoMesonQA;							// flag for meson QA
		Int_t 								fDoPhotonQA;						// flag for photon QA
		Int_t 								fDoClusterQA;						// flag for cluster QA
		Bool_t 								fIsFromMBHeader;					// flag for MC headers
		Bool_t								fIsOverlappingWithOtherHeader; 		// flag for particles in MC overlapping between headers
		Bool_t 								fIsMC;								// flag for MC information
		Bool_t 								fDoTHnSparse;                 // flag for THnSparse

		
	private:
		AliAnalysisTaskGammaCaloDalitzV1(const AliAnalysisTaskGammaCaloDalitzV1&); // Prevent copy-construction
		AliAnalysisTaskGammaCaloDalitzV1 &operator=(const AliAnalysisTaskGammaCaloDalitzV1&); // Prevent assignment

        ClassDef(AliAnalysisTaskGammaCaloDalitzV1, 5);
};

#endif
