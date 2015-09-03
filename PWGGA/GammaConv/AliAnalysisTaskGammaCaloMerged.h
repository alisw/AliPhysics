#ifndef ALIANLYSISTASKGAMMACALOMERGED_cxx
#define ALIANLYSISTASKGAMMACALOMERGED_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include <vector>
#include <map>

class AliAnalysisTaskGammaCaloMerged : public AliAnalysisTaskSE {
	public:

		AliAnalysisTaskGammaCaloMerged();
		AliAnalysisTaskGammaCaloMerged(const char *name);
		virtual ~AliAnalysisTaskGammaCaloMerged();

		virtual void   UserCreateOutputObjects();
		virtual Bool_t Notify();
		virtual void   UserExec(Option_t *);
		virtual void   Terminate(const Option_t*);

		void SetIsHeavyIon(Int_t flag){
			fIsHeavyIon = flag;    
		}

		// base functions for selecting photon and meson candidates in reconstructed data
		void ProcessClusters();
		void CalculatePi0Candidate(AliAODConversionPhoton* photon1, AliAODConversionPhoton* photon2);
		
		// MC functions
		void SetIsMC(Int_t isMC){fIsMC=isMC;}
		void ProcessMCParticles();
		void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate, Float_t m02, AliAODConversionPhoton *TrueSubClusterCandidate1,
															      AliAODConversionPhoton *TrueSubClusterCandidate2);
// // 		void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		
		// switches for additional analysis streams or outputs
		void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
		void SetDoClusterQA(Int_t flag){fDoClusterQA = flag;}
		void SetPlotHistsExtQA(Bool_t flag){fSetPlotHistsExtQA = flag;}
		
	    // Setting the cut lists for the conversion photons
		void SetEventCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fEventCutArray = CutArray;
		}

	    // Setting the cut lists for the calo photons
		void SetCaloCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fClusterCutArray = CutArray;
		}

	    // Setting the cut lists for the calo photons
		void SetCaloMergedCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fClusterMergedCutArray = CutArray;
		}
		
		// Setting the cut lists for the meson
		void SetMesonCutList(Int_t nCuts, TList *CutArray){
			fnCuts = nCuts;
			fMesonCutArray = CutArray;
		}
		
		// Additional functions for convenience
		void SetLogBinningXTH2(TH2* histoRebin);

		Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
		void FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked);
		void FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist);
		
	protected:
		AliV0ReaderV1 						*fV0Reader;							// basic photon Selection Task
		AliVEvent 							*fInputEvent;						// current event
		AliMCEvent 							*fMCEvent;							// corresponding MC event
		AliStack 							*fMCStack;							// stack belonging to MC event
		TList 								**fCutFolder;						// Array of lists for containers belonging to cut
		TList 								**fESDList;							// Array of lists with histograms with reconstructed properties   
		TList 								**fTrueList;						// Array of lists with histograms with MC validated reconstructed properties
		TList 								**fMCList;							// Array of lists with histograms with pure MC information
		TList 								**fHeaderNameList;					// Array of lists with header names for MC header selection
		TList 								*fOutputContainer;					// Output container
		Int_t 								fNClusterCandidates; 				//! current number of cluster candidates
		Int_t 								fNClusterMergedCandidates; 			//! current number of merged cluster candidates
		TList 								*fEventCutArray;					// List with Event Cuts
		AliConvEventCuts 					*fEventCuts;						// EventCutObject
		TList 								*fClusterCutArray;					// List with Cluster Cuts
		TList 								*fClusterMergedCutArray;			// List with Cluster Cuts for merged clusters
		TList 								*fMesonCutArray;					// List with Meson Cuts
		AliConversionMesonCuts 				*fMesonCuts;						// MesonCutObject
		
		//histograms for mesons reconstructed quantities
		TH2F 								**fHistoMotherInvMassPt;			//! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
		TH2F 								**fHistoMotherPi0PtY;				//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, Y
		TH2F 								**fHistoMotherEtaPtY;				//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, Y
		TH2F 								**fHistoMotherPi0PtAlpha;			//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, alpha
		TH2F 								**fHistoMotherEtaPtAlpha;			//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, alpha
		TH2F 								**fHistoMotherPi0PtOpenAngle;		//! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, openAngle
		TH2F 								**fHistoMotherEtaPtOpenAngle;		//! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, openAngle

		// histograms for rec photon clusters
		TH1F								** fHistoClusGammaPt;				//! array of histos with cluster, pt
		TH1F								** fHistoClusOverlapHeadersGammaPt;	//! array of histos with cluster, pt overlapping with other headers
		TH2F								** fHistoClusMergedPtvsM02;			//! array of histos with cluster merged, pt vs M02

		//histograms for pure MC quantities
		TH1I 								**fHistoMCHeaders;					//! array of histos for header names
		TH1F 								**fHistoMCPi0Pt;					//! array of histos with weighted pi0, pT
		TH1F 								**fHistoMCPi0WOWeightPt;			//! array of histos with unweighted pi0, pT
		TH1F 								**fHistoMCPi0WOEvtWeightPt;			//! array of histos without event weights pi0, pT
		TH1F 								**fHistoMCEtaPt;					//! array of histos with weighted eta, pT
		TH1F 								**fHistoMCEtaWOWeightPt;			//! array of histos with unweighted eta, pT
		TH1F 								**fHistoMCEtaWOEvtWeightPt;			//! array of histos without event weights eta, pT
		TH1F 								**fHistoMCPi0InAccPt;				//! array of histos with weighted pi0 in acceptance, pT
		TH1F 								**fHistoMCEtaInAccPt;				//! array of histos with weighted eta in acceptance, pT
		TH2F 								**fHistoMCPi0PtY;					//! array of histos with weighted pi0, pT, Y
		TH2F 								**fHistoMCEtaPtY;					//! array of histos with weighted eta, pT, Y
		TH2F 								**fHistoMCPi0PtAlpha;				//! array of histos with weighted pi0, pT, alpha
		TH2F 								**fHistoMCEtaPtAlpha;				//! array of histos with weighted eta, pT, alpha
		TH2F 								**fHistoMCPi0PtJetPt;				//! array of histos with weighted pi0, pT, hardest jet pt
		TH2F 								**fHistoMCEtaPtJetPt;				//! array of histos with weighted eta, pT, hardest jet pt

		// MC validated cluster histos
		TH2F								** fHistoTrueClusMergedPtvsM02;					//! 
		TH2F								** fHistoTrueClusPi0PtvsM02;					//! 
		TH2F								** fHistoTrueClusPrimPi0PtvsM02;				//! 
		TH2F								** fHistoTrueClusSecPi0PtvsM02;					//! 
		TH2F								** fHistoTrueClusSecPi0FromK0sPtvsM02;			//! 
		TH2F								** fHistoTrueClusSecPi0FromLambdaPtvsM02;		//! 
		TH2F								** fHistoTrueClusEtaPtvsM02;					//! 
		TH2F								** fHistoTrueClusMergedPartConvPtvsM02;			//! 
		TH2F								** fHistoTrueClusMergedPartConvELeadPtvsM02;	//! 
		TH2F								** fHistoTrueClusPartConvPi0PtvsM02;			//! 
		TH2F								** fHistoTrueClusPartConvPrimPi0PtvsM02;		//! 
		TH2F								** fHistoTrueClusPartConvSecPi0PtvsM02;			//! 
		TH2F								** fHistoTrueClusPartConvSecPi0FromK0sPtvsM02;	//! 
		TH2F								** fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02;	//! 
		TH2F								** fHistoTrueClusPartConvEtaPtvsM02;			//! 
		TH2F								** fHistoTrueClusBGPtvsM02;						//! 
		TH2F								** fHistoTrueClusGammaPtvsM02;					//! 
		TH2F								** fHistoTrueClusMergedInvMassvsPt;					//! 
		TH2F								** fHistoTrueClusPi0InvMassvsPt;					//! 
		TH2F								** fHistoTrueClusPrimPi0InvMassvsPt;				//! 
		TH2F								** fHistoTrueClusSecPi0InvMassvsPt;					//! 
		TH2F								** fHistoTrueClusSecPi0FromK0sInvMassvsPt;			//! 
		TH2F								** fHistoTrueClusSecPi0FromLambdaInvMassvsPt;		//! 
		TH2F								** fHistoTrueClusEtaInvMassvsPt;					//! 
		TH2F								** fHistoTrueClusMergedPartConvInvMassvsPt;			//! 
		TH2F								** fHistoTrueClusMergedPartConvELeadInvMassvsPt;	//! 
		TH2F								** fHistoTrueClusPartConvPi0InvMassvsPt;			//! 
		TH2F								** fHistoTrueClusPartConvPrimPi0InvMassvsPt;		//! 
		TH2F								** fHistoTrueClusPartConvSecPi0InvMassvsPt;			//! 
		TH2F								** fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt;	//! 
		TH2F								** fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt;	//! 
		TH2F								** fHistoTrueClusPartConvEtaInvMassvsPt;			//! 
		TH2F								** fHistoTrueClusBGInvMassvsPt;						//! 
		TH2F								** fHistoTrueClusGammaInvMassvsPt;					//! 
		TH2F								** fHistoTrueClusBGPtvsSource;						//!
		
		
		// MC validated reconstructed quantities mesons
		TH2F 								**fHistoTruePi0PtY;								//! array of histos with validated pi0, pt, Y
		TH2F 								**fHistoTrueEtaPtY;								//! array of histos with validated eta, pt, Y
		TH2F 								**fHistoTruePi0PtAlpha;							//! array of histos with validated pi0, pt, alpha
		TH2F 								**fHistoTrueEtaPtAlpha;							//! array of histos with validated eta, pt, alpha
		TH2F 								**fHistoTruePi0PtOpenAngle;						//! array of histos with validated pi0, pt, openAngle
		TH2F 								**fHistoTrueEtaPtOpenAngle;						//! array of histos with validated eta, pt, openAngle

		// MC validated reconstructed quantities photons
		TH2F				 				**fHistoDoubleCountTruePi0InvMassPt;			//! array of histos with double counted pi0s, invMass, pT
		TH2F				 				**fHistoDoubleCountTrueEtaInvMassPt;			//! array of histos with double counted etas, invMass, pT
		vector<Int_t>						fVectorDoubleCountTruePi0s;						//! vector containing labels of validated pi0
		vector<Int_t>						fVectorDoubleCountTrueEtas;						//! vector containing labels of validated eta

		// event histograms
		TH1F 								**fHistoNEvents;								//! array of histos with event information
		TH1F 								**fHistoNEventsWOWeight;						//! array of histos with event information without event weights
		TH1F 								**fHistoNGoodESDTracks;							//! array of histos with number of good tracks (2010 Standard track cuts)
		TH1F								**fHistoVertexZ;								//! array of histos with vertex z distribution for selected events
		TH1F 								**fHistoNClusterCandidates;						//! array of histos with number of cluster candidates per event
		TH1F 								**fHistoNClusterMergedCandidates;				//! array of histos with number of merged cluster candidates per event
		TH2F 								**fHistoNGoodESDTracksVsNClusterCandidates;		//! array of histos with number of good tracks vs gamma candidates
		TH2F								**fHistoSPDClusterTrackletBackground;			//! array of histos with SPD tracklets vs SPD clusters for background rejection
		TH1F 								**fHistoNV0Tracks;								//! array of histos with V0 counts
		TProfile 							**fProfileEtaShift;								//! array of profiles with eta shift
		TProfile							**fProfileJetJetXSection;						//! array of profiles with xsection for jetjet
		TH1F								**fHistoJetJetNTrials;							//! array of histos with ntrials for jetjet

		// additional variables
		TRandom3 							fRandom;							// random 
		Int_t 								fnCuts;								// number of cuts to be analysed in parallel
		Int_t 								fiCut;								// current cut	
		Int_t	 							fIsHeavyIon;						// switch for pp = 0, PbPb = 1, pPb = 2
		Int_t 								fDoMesonQA;							// flag for meson QA
		Int_t 								fDoClusterQA;						// flag for cluster QA
		Bool_t 								fIsFromMBHeader;					// flag for MC headers
		Bool_t								fIsOverlappingWithOtherHeader; 		// flag for particles in MC overlapping between headers
		Int_t 								fIsMC;								// flag for MC information
		Bool_t								fSetPlotHistsExtQA;					// flag for extended QA hists
		Double_t 							fWeightJetJetMC;					// weight for Jet-Jet MC

	private:
		AliAnalysisTaskGammaCaloMerged(const AliAnalysisTaskGammaCaloMerged&); // Prevent copy-construction
		AliAnalysisTaskGammaCaloMerged &operator=(const AliAnalysisTaskGammaCaloMerged&); // Prevent assignment

		ClassDef(AliAnalysisTaskGammaCaloMerged, 2);
};

#endif
