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
		void CalculatePi0Candidates();
		
		// MC functions
		void SetIsMC(Int_t isMC){fIsMC=isMC;}
		void ProcessMCParticles();
		void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate);
		void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
		
		// switches for additional analysis streams or outputs
		void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
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
		TList 								*fClusterCandidates; 				//! current list of cluster candidates
		TList 								*fEventCutArray;					// List with Event Cuts
		AliConvEventCuts 					*fEventCuts;						// EventCutObject
		TList 								*fClusterCutArray;					// List with Cluster Cuts
		AliCaloPhotonCuts 					*fCaloPhotonCuts;					// CaloPhotonCutObject
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

		// MC validated reconstructed quantities mesons
		TH2F 								**fHistoTruePi0InvMassPt;					//! array of histos with validated mothers, invMass, pt
		TH2F 								**fHistoTrueEtaInvMassPt;					//! array of histos with validated mothers, invMass, pt
		TH2F 								**fHistoTruePrimaryPi0InvMassPt;				//! array of histos with validated weighted primary mothers, invMass, pt
		TH2F 								**fHistoTruePrimaryEtaInvMassPt;				//! array of histos with validated weighted primary mothers, invMass, pt
		TH2F 								**fHistoTruePrimaryPi0W0WeightingInvMassPt;		//! array of histos with validated unweighted primary mothers, invMass, pt
		TH2F 								**fHistoTruePrimaryEtaW0WeightingInvMassPt;		//! array of histos with validated unweighted primary mothers, invMass, pt
		TProfile2D 							**fProfileTruePrimaryPi0WeightsInvMassPt;		//! array of profiles with weights for validated primary mothers, invMass, pt	
		TProfile2D 							**fProfileTruePrimaryEtaWeightsInvMassPt;		//! array of profiles with weights for validated primary mothers, invMass, pt	
		TH2F 								**fHistoTruePrimaryPi0MCPtResolPt;				//! array of histos with validated weighted primary pi0, MCpt, resol pt
		TH2F	 							**fHistoTruePrimaryEtaMCPtResolPt;				//! array of histos with validated weighted primary eta, MCpt, resol pt
		TH2F 								**fHistoTrueSecondaryPi0InvMassPt;				//! array of histos with validated secondary mothers, invMass, pt
		TH2F 								**fHistoTrueSecondaryPi0FromK0sInvMassPt;		//! array of histos with validated secondary mothers from K0s, invMass, pt
		TH1F 								**fHistoTrueK0sWithPi0DaughterMCPt;				//! array of histos with K0s with reconstructed pi0 as daughter, pt
		TH2F 								**fHistoTrueSecondaryPi0FromLambdaInvMassPt;	//! array of histos with validated secondary mothers from Lambda, invMass, pt
		TH1F 								**fHistoTrueLambdaWithPi0DaughterMCPt;			//! array of histos with lambda with reconstructed pi0 as daughter, pt
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
		TH2F 								**fHistoNGoodESDTracksVsNClusterCandidates;		//! array of histos with number of good tracks vs gamma candidates
		TH2F								**fHistoSPDClusterTrackletBackground;			//! array of histos with SPD tracklets vs SPD clusters for background rejection
		TH1F 								**fHistoNV0Tracks;								//! array of histos with V0 counts
		TProfile 							**fProfileEtaShift;								//! array of profiles with eta shift
				
		// additional variables
		TRandom3 							fRandom;							// random 
		Int_t 								fnCuts;								// number of cuts to be analysed in parallel
		Int_t 								fiCut;								// current cut	
		Int_t	 							fIsHeavyIon;						// switch for pp = 0, PbPb = 1, pPb = 2
		Bool_t 								fDoMesonAnalysis;					// flag for meson analysis
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

		ClassDef(AliAnalysisTaskGammaCaloMerged, 1);
};

#endif
