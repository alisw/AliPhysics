#ifndef ALICALOPHOTONCUTS_H
#define ALICALOPHOTONCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner                                    *

#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TF1.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "AliEMCALGeometry.h"
#include "AliPHOSGeometry.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class TH1F;
class TH2F;
class TF1;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;

using namespace std;

class AliCaloPhotonCuts : public AliAnalysisCuts {
		
	public: 
		enum cutIds {
			kClusterType,
			kEtaMin,
			kEtaMax,
			kPhiMin,
			kPhiMax,
			kNonLinearity1,
			kNonLinearity2,
			kDistanceToBadChannel,
			kTiming,
			kTrackMatching,
			kExoticCell,
			kMinEnergy,
			kNMinCells,
			kMinM02,
			kMaxM02,
			kMinM20,
			kMaxM20,
			kDispersion,
			kNLM,
			kNCuts
		};

		enum photonCuts {
			kPhotonIn=0,
			kDetector,
			kAcceptance,
			kClusterQuality,
			kPhotonOut
		};

		//handeling of CutString
		static const char * fgkCutNames[kNCuts];
		Bool_t 			SetCutIds(TString cutString); 
		Int_t 			fCuts[kNCuts];
		Bool_t 			SetCut(cutIds cutID, Int_t cut);
		Bool_t 			UpdateCutString();
		void 			PrintCuts();
		void 			PrintCutsWithValues();
		
		Bool_t 			InitializeCutsFromCutString(const TString analysisCutSelection);
		TString 		GetCutNumber();
		Int_t 			GetClusterType() {return fClusterType;}
		
		//Constructors
		AliCaloPhotonCuts(const char *name="ClusterCuts", const char * title="Cluster Cuts");
		AliCaloPhotonCuts(const AliCaloPhotonCuts&);
		AliCaloPhotonCuts& operator=(const AliCaloPhotonCuts&);

		//virtual destructor
		virtual 		~AliCaloPhotonCuts();                          

		virtual Bool_t 	IsSelected(TObject* /*obj*/)									{return kTRUE;}
		virtual Bool_t 	IsSelected(TList* /*list*/) 									{return kTRUE;}

		Bool_t 			ClusterIsSelected(AliVCluster* cluster, AliVEvent *event, Int_t isMC);
		Bool_t 			ClusterIsSelectedMC(TParticle *particle,AliStack *fMCStack);
		Bool_t 			ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray);
			
		//correct NonLinearity
		void			CorrectEMCalNonLinearity(AliVCluster* cluster, Int_t isMC);
		Float_t			FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6);
		Float_t			FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);

		void 			InitCutHistograms(TString name="");
		void 			SetFillCutHistograms(TString name="")							{if(!fHistograms){InitCutHistograms(name);} return;}
		TList*			GetCutHistograms()												{return fHistograms;}
		TList*			GetExtQAHistograms()											{return fHistExtQA;}
		void 			FillClusterCutIndex(Int_t photoncut)							{if(fHistCutIndex)fHistCutIndex->Fill(photoncut); return;}

		void 			SetExtendedMatchAndQA(Int_t extendedMatchAndQA)					{fExtendedMatchAndQA = extendedMatchAndQA; return;}
		void			SetExtendedQA(Int_t extendedQA)									{if(extendedQA != 1 && extendedQA != 2)fExtendedMatchAndQA = extendedQA; return;}
		void			FillHistogramsExtendedQA(AliVEvent *event);

		// Cut functions
		Bool_t 			AcceptanceCuts(AliVCluster* cluster, AliVEvent *event);
		Bool_t 			ClusterQualityCuts(AliVCluster* cluster,AliVEvent *event, Int_t isMC);

		Bool_t 			MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event);

		// Set Individual Cuts
		Bool_t 			SetClusterTypeCut(Int_t);
		Bool_t 			SetMinEtaCut(Int_t);
		Bool_t 			SetMaxEtaCut(Int_t);
		Bool_t 			SetMinPhiCut(Int_t);		
		Bool_t 			SetMaxPhiCut(Int_t);		
		Bool_t 			SetDistanceToBadChannelCut(Int_t);
		Bool_t 			SetTimingCut(Int_t);
		Bool_t 			SetTrackMatchingCut(Int_t);
		Bool_t 			SetExoticCellCut(Int_t);
		Bool_t 			SetMinEnergyCut(Int_t);
		Bool_t 			SetMinNCellsCut(Int_t);
		Bool_t 			SetMaxM02(Int_t);
		Bool_t 			SetMinM02(Int_t);
		Bool_t 			SetMaxM20(Int_t);
		Bool_t 			SetMinM20(Int_t);
		Bool_t 			SetDispersion(Int_t);
		Bool_t 			SetNLM(Int_t);
		Bool_t			SetNonLinearity1(Int_t);
		Bool_t			SetNonLinearity2(Int_t);
		
	protected:
		TList			*fHistograms;
		TList			*fHistExtQA;

		AliEMCALGeometry	*geomEMCAL;					// pointer to EMCAL geometry
		AliPHOSGeometry		*geomPHOS;					// pointer to PHOS geometry
		TObjArray*			EMCALBadChannelsMap;		// pointer to EMCAL bad channel map
		TProfile*			BadChannels;				// TProfile with bad channels
		Int_t				nMaxEMCalModules;			// max number of EMCal Modules
		Int_t				nMaxPHOSModules;			// max number of PHOS Modules
		
		//cuts
		Int_t		fClusterType;						// which cluster do we have
		Double_t 	fMinEtaCut; 						// min eta cut
		Double_t 	fMaxEtaCut; 						// max eta cut
		Bool_t 		fUseEtaCut;							// flag for switching on eta cut
		Double_t 	fMinPhiCut; 						// phi cut
		Double_t 	fMaxPhiCut; 						// phi cut
		Bool_t 		fUsePhiCut;							// flag for switching on phi cut
		Double_t 	fMinDistanceToBadChannel; 			// minimum distance to bad channel
		Bool_t 		fUseDistanceToBadChannel;			// flag for switching on distance to bad channel cut
		Double_t 	fMaxTimeDiff; 						// maximum time difference to triggered collision
		Double_t	fMinTimeDiff;						// minimum time difference to triggered collision
		Bool_t 		fUseTimeDiff;						// flag for switching on time difference cut
        Double_t 	fMaxDistTrackToClusterEta; 			// minimum distance between track and cluster in eta
        Double_t 	fMinDistTrackToClusterPhi; 			// minimum distance between track and cluster in phi
        Double_t    fMaxDistTrackToClusterPhi;          // maximum distance between track and cluster in phi
		Bool_t 		fUseDistTrackToCluster;				// flag for switching on distance between track and cluster cut
		Int_t 		fExtendedMatchAndQA;				// switching on extended matching histograms (1) + extended QA (2) / only extended QA (3) or both off (0)
		Double_t 	fExoticCell;						// exotic cell cut
		Bool_t 		fUseExoticCell;						// flag for switching on exotic cell cut
		Double_t 	fMinEnergy;							// minium energy per cluster
		Bool_t 		fUseMinEnergy;						// flag for switching on minimum energy cut
		Int_t 		fMinNCells;							// minimum number of cells 
		Bool_t 		fUseNCells;							// flag for switching on minimum N Cells cut
		Double_t 	fMaxM02;							// maximum M02
		Double_t 	fMinM02;							// minimum M02
		Bool_t 		fUseM02;							// flag for switching on M02 cut
		Double_t 	fMaxM20;							// maximum M20
		Double_t 	fMinM20;							// minimum M20
		Bool_t 		fUseM20;							// flag for switching on M20 cut
		Double_t 	fMaxDispersion;						// maximum dispersion
		Bool_t 		fUseDispersion;						// flag for switching on dispersion cut
		Int_t 		fMinNLM;							// minimum number of local maxima in cluster
		Int_t 		fMaxNLM;							// maximum number of local maxima in cluster
		Bool_t 		fUseNLM;							// flag for switching on NLM cut
		Int_t		fNonLinearity1;						// selection of nonlinearity correction, part1
		Int_t		fNonLinearity2;						// selection of nonlinearity correction, part2
		Int_t		fSwitchNonLinearity;				// selection (combined) of NonLinearity
		Bool_t		fUseNonLinearity;					// flag for switching NonLinearity correction
		
		// CutString
		TObjString* fCutString; 							// cut number used for analysis
		
		// Histograms
		TH1F* 		fHistCutIndex; 							// bookkeeping for cuts
		TH1F* 		fHistAcceptanceCuts; 					// bookkeeping for acceptance cuts
		TH1F* 		fHistClusterIdentificationCuts; 		// bookkeeping for cluster identification cuts
		
		TH2F* 		fHistClusterEtavsPhiBeforeAcc; 			// eta-phi-distribution before acceptance cuts
		TH2F* 		fHistClusterEtavsPhiAfterAcc; 			// eta-phi-distribution of all after acceptance cuts
		TH2F* 		fHistClusterEtavsPhiAfterQA; 			// eta-phi-distribution of all after cluster quality cuts
// 		TH1F* 		fHistDistanceToBadChannelBeforeAcc;		// distance to bad channel before acceptance cuts
// 		TH1F* 		fHistDistanceToBadChannelAfterAcc;		// distance to bad channel after acceptance cuts
		TH2F* 		fHistClusterTimevsEBeforeQA;			// Cluster time vs E before cluster quality cuts
		TH2F* 		fHistClusterTimevsEAfterQA;				// Cluster time vs E after cluster quality cuts
// 		TH2F* 		fHistExoticCellBeforeQA;				// Exotic cell: 1-Ecross/E cell vs Ecluster before acceptance cuts
// 		TH2F* 		fHistExoticCellAfterQA;					// Exotic cell: 1-Ecross/E cell vs Ecluster after cluster quality cuts
// 		TH1F* 		fHistNMatchedTracks;					// number of matched tracks
		TH1F* 		fHistEnergyOfClusterBeforeNL;			// enery per cluster before NonLinearity correction
		TH1F* 		fHistEnergyOfClusterAfterNL;			// enery per cluster after NonLinearity correction
		TH1F* 		fHistEnergyOfClusterBeforeQA;			// enery per cluster before acceptance cuts
		TH1F* 		fHistEnergyOfClusterAfterQA;			// enery per cluster after cluster quality cuts
		TH1F* 		fHistNCellsBeforeQA;					// number of cells per cluster before acceptance cuts
		TH1F* 		fHistNCellsAfterQA;						// number of cells per cluster after cluster quality cuts
		TH1F* 		fHistM02BeforeQA;						// M02 before acceptance cuts
		TH1F* 		fHistM02AfterQA;						// M02 after cluster quality cuts
		TH1F* 		fHistM20BeforeQA;						// M20 before acceptance cuts
		TH1F* 		fHistM20AfterQA;						// M20 after cluster quality cuts
		TH1F* 		fHistDispersionBeforeQA;				// dispersion before acceptance cuts
		TH1F* 		fHistDispersionAfterQA;					// dispersion after cluster quality cuts
// 		TH1F* 		fHistNLMBeforeQA;						// number of local maxima in cluster before acceptance cuts
// 		TH1F* 		fHistNLMAfterQA;						// number of local maxima in cluster after cluster quality cuts

		//More histograms
		TH2F*		fHistClusterEnergyvsMod;				// Cluster Energy vs Module Number
		TH2F*		fHistNCellsBigger100MeVvsMod;			// NCells with >0.1 GeV vs Module Number
		TH2F*		fHistNCellsBigger1500MeVvsMod;			// NCells with >1.5 GeV vs Module Number
		TH2F*		fHistEnergyOfModvsMod;					// Deposited Energy vs Module Number
		TH2F*		fHistClusterEnergyvsNCells;				// Cluster Energy vs NCells
		TH2F*		fHistCellEnergyvsCellID;				// Cell Energy vs CellID
		TH2F*		fHistCellTimevsCellID;					// Cell Time vs CellID
		TH2F* 		fHistClusterEM02BeforeQA;				// 2-dim plot E vs. M02
		TH2F* 		fHistClusterEM02AfterQA;				// 2-dim plot E vs. M02
		TH1F*		fHistClusterIncludedCellsBeforeQA;		// CellIDs in Cluster
		TH1F*		fHistClusterIncludedCellsAfterQA;		// CellIDs in Cluster of accepted ones
		TH1F*		fHistClusterEnergyFracCellsBeforeQA;	// Energy fraction of CellIDs in Cluster
		TH1F*		fHistClusterEnergyFracCellsAfterQA;		// Energy fraction of CellIDs in Cluster of accepted ones

        //Track matching histograms
		TH1F* 		fHistClusterRBeforeQA;					// cluster position in R=SQRT(x^2+y^2) (before QA)
		TH1F* 		fHistClusterRAfterQA;					// cluster position in R=SQRT(x^2+y^2) for matched tracks (After QA)
        TH2F* 		fHistClusterdEtadPhiBeforeQA;			// 2-dim plot dEta vs. dPhi
		TH2F* 		fHistClusterdEtadPhiAfterQA;			// 2-dim plot dEta vs. dPhi for matched tracks (after QA)
		TH1F* 		fHistDistanceTrackToClusterBeforeQA;	// distance cluster to track before acceptance cuts
        TH1F* 		fHistDistanceTrackToClusterAfterQA;		// distance cluster to track after cluster quality cuts

        //Extended track matching histograms
        TH2F* 		fHistClusterdEtadPhiPosTracksBeforeQA;	// 2-dim plot dEta vs. dPhi
        TH2F* 		fHistClusterdEtadPhiNegTracksBeforeQA;	// 2-dim plot dEta vs. dPhi
        TH2F* 		fHistClusterdEtadPhiPosTracksAfterQA;	// 2-dim plot dEta vs. dPhi
        TH2F* 		fHistClusterdEtadPhiNegTracksAfterQA;	// 2-dim plot dEta vs. dPhi
        TH2F*       fHistClusterdEtadPhiPosTracksP_000_075BeforeQA;// 2-dim plot dEta vs. dPhi, positive Tracks, P < 0.75
        TH2F*       fHistClusterdEtadPhiPosTracksP_075_125BeforeQA;// 2-dim plot dEta vs. dPhi, positive Tracks, 0.75 < P < 1.25
        TH2F*       fHistClusterdEtadPhiPosTracksP_125_999BeforeQA;// 2-dim plot dEta vs. dPhi, positive Tracks, P > 1.25
        TH2F*       fHistClusterdEtadPhiNegTracksP_000_075BeforeQA;// 2-dim plot dEta vs. dPhi, negative Tracks, P < 0.75
        TH2F*       fHistClusterdEtadPhiNegTracksP_075_125BeforeQA;// 2-dim plot dEta vs. dPhi, negative Tracks, 0.75 < P < 1.25
        TH2F*       fHistClusterdEtadPhiNegTracksP_125_999BeforeQA;// 2-dim plot dEta vs. dPhi, negative Tracks, P > 1.25
        TH2F* 		fHistClusterdEtadPtBeforeQA;			// 2-dim plot dEta vs. Pt
        TH2F* 		fHistClusterdPhidPtBeforeQA;			// 2-dim plot dEta vs. Pt
        TH2F* 		fHistClusterM20M02BeforeQA;				// 2-dim plot M20 vs. M02
        TH2F* 		fHistClusterM20M02AfterQA;				// 2-dim plot M20 vs. M20

	private:

		ClassDef(AliCaloPhotonCuts,7)
};

#endif
