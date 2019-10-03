#ifndef ALIANALYSISTASKPHOSNEUTRALMESON_H
#define ALIANALYSISTASKPHOSNEUTRALMESON_H

/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                            */

/* $Id$ */

/////---------------------------------------------------------------------------------------------------
/////     AnalyisTask for neutral meson analysis with PHOS 
/////     Runs on ESDs and AODs, Tested on pp and pPb data.
/////     Authors: Malte Hecker, Fabian Pliquett
/////     Date: 01/12/2015
/////
/////-------------------------------------------------------------------------------------------------

class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TH3D;
class TString;
class TNtuple;
class TList;
class AliAnalysisUtils;
class AliESDEvent;
class AliAODEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliESDCaloCluster;
class AliAODCaloCluster;
class AliMCEvent;
class AliMCParticle;
class AliPHOSGeometry;
class AliTriggerAnalysis;

class AliGenHijingEventHeader;
class AliGenPythiaEventHeader;
class AliGenEventHeader;


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#include "AliPHOSEsdCluster.h"
#endif

class AliAnalysisTaskPHOSNeutralMeson : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskPHOSNeutralMeson() ;
		AliAnalysisTaskPHOSNeutralMeson(const char *name);
		virtual ~AliAnalysisTaskPHOSNeutralMeson();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);
		virtual void     Terminate(Option_t *);

		void SetPHOSBadMap(Int_t mod,TH2I * h)
		{
			if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
			fPHOSBadMap[mod]=new TH2I(*h) ;
			printf("Set %s \n",fPHOSBadMap[mod]->GetName());
		}
		
		void SetApplyBadMapManually(Bool_t option)  
		{
			fApplyBadMapManually = option;
		}

		void SetClusterMinCells(Int_t nCells)
		{
			fClusterMinCells = nCells;
		}

		void SetClusterMinE(Float_t energy)
		{
			fClusterMinE = energy;
		}

		void SetClusterMinM02(Float_t dispersion)
		{
			fClusterMinM02 = dispersion;
		}	
		
		void SetDoDistToBadCellCutOnCellLevel(Bool_t option)
		{
			// Uses self written function IsGoodChannelDistanceOnCellLevel(const char * det, Int_t mod,Int_t ix, Int_t iz, Int_t dist);
			fDoDistToBadCellCutOnCellLevel = option; 
		}
		
		void SetDistToBadCellOnCellLevel(Int_t ncells)
		{	
			fDistToBadCellOnCellLevel = ncells;

		}
		
		void SetDoDistToBadCellCut(Bool_t option)
		{
			//(in cm)
			//uses virCluster->GetDistanceToBadChannel() > fDistToBadCell   
			fDoDistToBadCellCut = option; 
		}
		
		void SetDistToBadCell(Float_t distance)
		{	
			fDistToBadCell = distance;
		}

		void SetDoTimingCut(Bool_t option)
		{
			fDoTimingCut = option;
		}

		void SetTimingCutMinMax(Float_t min, Float_t max)
		{
			fTimingCutMin = min;
			fTimingCutMax = max;
		}

		void SetEtaAccMinMax(Double_t min, Double_t max)
		{
			fEtaAccMin = min; 
			fEtaAccMax = max; 
		}

		void SetPhiAccMinMax(Double_t min, Double_t max)
		{
			fPhiAccMin = min; 
			fPhiAccMax = max; 
		}

		void SetFillCellIdVsE(Bool_t option)
		{
			fFillCellIdVsE = option;
		}

		void SetFillTimingHistos(Bool_t option)
		{
			fFillTimingHistos = option;
		}

		void SetUseOnlyCINT1events(Bool_t option)
		{
			fUseOnlyCINT1Events = option;
		}

		void SetFillClusterHitmaps(Bool_t option)
		{
			fFillClusterHitmaps = option;
		}

		void SetFillFeedDownHistos(Bool_t option)
		{
			fFillFeedDownHistos = option; 
		}

		void SetDoClusterEnergyRecalibration(Bool_t option)
		{
			fDoClusterEnergyRecalibration = option; 
		}

		void SetRecalibrationOption(TString option) 
		{
			fRecalibrationOption = option;
		}

		void SetDoZvertexCut(Bool_t option)
		{
			fDoZvertexCut = option; 
		}

		void SetZvertexCut(Float_t maxzvtx)
		{
			fZvertexCut = maxzvtx; 
		}

		void SetFillhMassPtModules(Bool_t option) 
		{
			fFillHMassPtModules = option;
			if(option) fMPtHistoMode = "hMPTvsMod";
		}
		
		void SetFillhMassPtTiming(Bool_t option)
		{
			fFillHMassPtTiming = option; 
			if(option) fMPtHistoMode = "hMPTvsTim";
		}
		
		void SetFillNewAsymmClasses(Bool_t option) 
		{
			fFillNewAsymmClasses = option;
			if(option) fMPtHistoMode = "hMPTnewAsy";
		}
		
		void SetNewAsymmClasses(Float_t val1, Float_t val2, Float_t val3)
		{
			fAsymmClass1 = val1;
			fAsymmClass2 = val2;
			fAsymmClass3 = val3; 
		}

		void SetFillNTupelClusterEnergy(Bool_t option) 
		{
			fFillNTupelClusterE = option;
		}

		void SetRecalibrateModuleWise(Bool_t option) 
		{
			// The cluster energies will be recalibrated with one factor per module
			// Provide those factors using SetModuleWiseRecalibrationFactors()
			// Also the cell energies will be recalibrated with these factors if you choose
			// to fill CellID_vs_Energy histos (turn on with SetFillCellIdVsE() )
			fRecalibrateModuleWise = option;
		}

		void SetModuleWiseRecalibrationFactors(Double_t factorMod1, Double_t factorMod2, Double_t factorMod3 )
		{
			fRecalFactorMod1 = factorMod1;
			fRecalFactorMod2 = factorMod2;
			fRecalFactorMod3 = factorMod3;
		}
		
		void SetDoPeakSmearing(Bool_t option)
		{
		      fDoPeakSmearing = option; 
		}
		
		void SetSmearFactor(Double_t sfactor)
		{
		      fSmearFactor = sfactor; 
		}

		void SetMcMode(Bool_t b) 
		{ 
			fMcMode = b;
		}

		void SetSelectedCollCandidatesForOutputName(TString option)
		{
			fSelectedCollisionCandidates = option;
		}
		void SetSelectedPhysSelectionForOutputName(TString option)
		{
			fSelectedPhysSelection = option;
		}
		void SetSelectedSelectedTenderForOutputName(TString option)
		{
			fSelectedTender = option;
		}
		void SetUsedBadmapForOutputName(TString option)
		{
			fUsedBadMap = option;
		}
		
		void SetUseIsVertexSelected2013pA(Bool_t option)
		{
			fUseIsVertexSelected2013pA = option;
		}
		
		void SetMaxMultforMultBinning(Int_t value) {
			fMaxMultForMultBinning = value;
		}
		void SetNMultBins(Int_t value) {
			fNMultBins = value;  // !!! maximum value is 30, otherwise array is out of range !!!
		}
		
		void SetAnalyseAddedSignals(Bool_t option) {
			//need to set fFillMCHistos=true as well!  (if(!fFillMCHistos) fAnalyseAddedSignals = kFALSE;)		
			fAnalyseAddedSignals = option;
		}
		
		void SetFillClusterHitmapsAddedSignals(Bool_t option) {
			//need to set fFillMCHistos=true as well!  (if(!fFillMCHistos) fAnalyseAddedSignals = kFALSE;)		
			fFillClusterHitmapsAddedSignals = option;
		}
		
		void SetFillDecayPhotonInfoAddedSig(Bool_t option) {
			//need to set fFillMCHistos=true as well!  (if(!fFillMCHistos) fAnalyseAddedSignals = kFALSE;)		
			fFillDecayPhotonInfoAddedSig = option;
		}
		
		void SetfFillMCHistos(Bool_t option) {
			//needs to be set to true to analyse MC!
			fFillMCHistos = option;
		}
		
		void SetTsalisParameters(Double_t param1, Double_t param2, Double_t param3, Double_t param4, Double_t param5) {
				//Sets the parameters used in tsalis function used in ::CalculateAddedSignalWeight(Double_t pT)
				fTsalisPi0Param1 = param1;
				fTsalisPi0Param2 = param2;
				fTsalisPi0Param3 = param3;
				fTsalisPi0Param4 = param4;
				fTsalisPi0Param5 = param5;
		}
		void SetExponentialParameters(Double_t param1, Double_t param2, Double_t param3, Double_t param4) {
				//Sets the parameters used in tsalis function used in ::CalculateAddedSignalWeight(Double_t pT)
				fExpParam1 = param1;
				fExpParam2 = param2;
				fExpParam3 = param3;
				fExpParam4 = param4;
		}
		
		void SetFillNDaughtersVsPtAndR(Bool_t option) {
			fFillNDaughtersVsPtAndR = option;
		}
		
		void SetFillMPtForSingleOrMultContrClus(Bool_t option) {
			fFillMPtForSingleOrMultContrClus = option;
		}
		
		void SetAnalyseMCPi0OrEta(TString particle) { 

			//Sets either the bool fAnalyseMCPi0 or fAnalyseMCEta or both to true.
			//Used to select which added signal should be analysed, to save memory if only one is needed.
			
			if(particle == "pi0"){
				fAnalyseMCPi0 = kTRUE;
				printf("Set MC particle to analyse to %s \n",particle.Data());
			} else if(particle == "eta" ) {
				fAnalyseMCEta = kTRUE;
				printf("Set MC particle  to analyse to %s \n",particle.Data());
			} else if (particle == "both"){
				fAnalyseMCPi0 = kTRUE;
				fAnalyseMCEta = kTRUE;
				printf("Set MC particle  to analyse to %s \n",particle.Data());
			} else {
				printf("ATTENTION: MC particle to analyse could not be set to %s \n",particle.Data());
			}
		}
		
		void SetAdditionalFileNameString(TString option) {
			//String appears at the end of the filname.
			//give short string to save any info you want to appear in the filename
			// string should start with "_"
			fAdditionalFileNameString = option;
		}
		
 
		TString MakeOutputFileName () {
			// ********************* IMPORTANT: ************************************
			// The OutputfileName cannot be to long!
			// The Gridjobs crash in a ERROR_RE if it is too long
			// It has been successfully tested with a filename length of 135 symbols
			// *********************************************************************
			
		
			//Example: 
			//PHS_kMB_AliPhySl_PHSTndr_2_0.30_0.00_0_-200ns_200ns_noZvtxCt_noNonLinCor_0.994_0.982_1.007_hMPTvsTim_BM_y_11a_2_9_2013_modded.root
				
			// !!!  Please keep the list of cuts/settings below up to date !!!
			//fSelectedCollisionCandidates
			//fSelectedPhysSelection
			//SelectedPHOSTender
			//Int_t fClusterMinCells;
			//Float_t fClusterMinE;
			//Float_t fClusterMinM02; 
			//Int_t fDistToBadCellOnCellLevel;
			//Int_t fDistToBadCell;
			//Float_t fTimingCutMin;
			//Float_t fTimingCutMax;	
			//Float_t fZvertexCut;
			//TString fRecalibrationOption; 
			//Double_t fRecalFactorMod1;
			//Double_t fRecalFactorMod2;
			//Double_t fRecalFactorMod3;
			//fMPtHistoMode = "fFillNewAsymmClasses" oder "fFillHMassPtModules" oder "fFillHMassPtTiming"
			//fApplyBadMapManually
			//fUseIsVertexSelected2013pA
			//fAnalyseAddedSignals
			//fUsedBadMap
			//fAdditionalFileNameString
			
			
				// Explanation
			// Some cuts are turned of using a Bool (e.g.: fDoTimingCut) but you can still set the cut variables
			// (e.g. fTimingCutMin / fTimingCutMax)  In that case the cut is still not applied but the OutputFilneName would be wrong
			// So these cuts need some special implementation here:
	
			TString zVertexCutString, timingCutMinString, timingCutMaxString;
	
			if(fDoZvertexCut) {
				zVertexCutString = Form("%.1f",fZvertexCut);
			} else {
				zVertexCutString = "noZvtxCt";
			}
			
			if(fDoTimingCut) {
				timingCutMinString = Form("%.0fns",fTimingCutMin / 0.000000001);  //express value in nanoseconds
				timingCutMaxString = Form("%.0fns",fTimingCutMax / 0.000000001);
			} else {
				timingCutMinString = "noTCt";
				timingCutMaxString = "noTCt";
			}
			
			fOutputFileName = Form("PHS_%s_%s_%s_%d_%.2f_%.2f_%d_%.2f_%s_%s_%s_%s_%.3f_%.3f_%.3f_%s_%d_%d_%d_%s%s.root",fSelectedCollisionCandidates.Data(),fSelectedPhysSelection.Data(),fSelectedTender.Data(),fClusterMinCells,fClusterMinE, fClusterMinM02, fDistToBadCellOnCellLevel, fDistToBadCell, timingCutMinString.Data(), timingCutMaxString.Data(), zVertexCutString.Data(), fRecalibrationOption.Data(), fRecalFactorMod1, fRecalFactorMod2, fRecalFactorMod3, fMPtHistoMode.Data(),fApplyBadMapManually,fUseIsVertexSelected2013pA,fAnalyseAddedSignals,fUsedBadMap.Data(),fAdditionalFileNameString.Data());
			return fOutputFileName;
		}
		
		
	private:
		static const Int_t fgkZvtxBins = 16;
		static const Int_t fgkMultBins = 30;
		static const UInt_t fgkPoolDepth = 20;

		Int_t GetMultBin(Int_t mult);
		Int_t GetZvtxBin(Double_t vertZ);
		Bool_t IsGoodChannel(const Char_t * det, Int_t mod,Int_t ix, Int_t iz);
		Bool_t IsGoodChannelDistanceOnCellLevel(const Char_t * det, Int_t mod,Int_t ix, Int_t iz, Int_t dist);
		Double_t GetDeltaPhi(TLorentzVector p1, TLorentzVector p2);
		Double_t GetDeltaEta(TLorentzVector p1, TLorentzVector p2);
		Double_t RecalibratePHOSClusterEnergy(TString calibOption, Double_t E, Int_t run); 
		Double_t CalculateAddedSignalWeightPi0Tsalis(Double_t pt);
		Double_t CalculateAddedSignalWeightPi0Exp(Double_t pt);
		Double_t CalculateAddedSignalWeightEta(Double_t pt);
		Int_t GetNDaughtersOfMCParticle(AliMCParticle * mcP);
		Bool_t DecayedToGammaGamma(AliMCParticle * mcP, AliMCEvent* mcEvent);
		Bool_t MCParticleIsPhysicalDecay(AliMCParticle * mcP, AliMCEvent* mcEvent);
		void FillHistoWithDaughterInfoPi0(AliMCParticle * mcP, TH1F* histo, AliMCEvent* mcEvent);
		void FillDecayGammaHistos(AliMCParticle * mcP, AliMCEvent* mcEvent,  Float_t wgtPi0,
																								TH1F* fH1DecGammAddPi0Eta,
																								TH1F* fH1DecGammAddPi0Y,
																								TH1F* fH1DecGammAddPi0Phi,
																								TH1F* fH1DecGammAddPi0E,
																								TH1F* fH1DecGammInPHOSAddPi0E,
																								TH1F* fH1DecGammAddPi0Asymm,
																								TH1F* fH1DecGammAddPi0OpAngle,
																								TH1F* fH1DecGammAddPi0ConvR,
																								TH1F* fH1DecGammAddPi0ConvRate,
																								TH1F* fH1TruthPtAddedPi0GammasInPHOS);

		TH2I *fPHOSBadMap[6]; 

		// ********** START Cut Variables ********** //
		Int_t fClusterMinCells;
		Float_t fClusterMinE;
		Float_t fClusterMinM02; 
		Bool_t fDoDistToBadCellCutOnCellLevel;
		Int_t fDistToBadCellOnCellLevel;
		Bool_t fDoDistToBadCellCut;
		Float_t fDistToBadCell;
		Bool_t fDoTimingCut;
		Float_t fTimingCutMin;
		Float_t fTimingCutMax;
		Double_t fEtaAccMin;
		Double_t fEtaAccMax;
		Double_t fPhiAccMin;
		Double_t fPhiAccMax;
		Bool_t fFillCellIdVsE;
		Bool_t fFillTimingHistos;
		Bool_t fUseOnlyCINT1Events; //neccesary for lhc11a
		Bool_t fFillClusterHitmaps;
		Bool_t fFillFeedDownHistos; 
		Bool_t fDoClusterEnergyRecalibration; 
		TString fRecalibrationOption; 
		Bool_t fDoZvertexCut;
		Float_t fZvertexCut;
		Bool_t fFillHMassPtModules; 
		Bool_t fFillHMassPtTiming;
		Bool_t fFillNewAsymmClasses;
		Float_t fAsymmClass1; 
		Float_t fAsymmClass2; 
		Float_t fAsymmClass3; 
		Bool_t fFillNTupelClusterE;
		Bool_t fRecalibrateModuleWise;
		Double_t fRecalFactorMod1;
		Double_t fRecalFactorMod2;
		Double_t fRecalFactorMod3;
		Bool_t fDoPeakSmearing; 
		Double_t fSmearFactor; 
		TString fSelectedCollisionCandidates;
		TString fSelectedPhysSelection;
		TString fSelectedTender;
		TString fOutputFileName;
		TString fMPtHistoMode;
		TString fUsedBadMap;
		Bool_t fUseIsVertexSelected2013pA;
		Bool_t fApplyBadMapManually;
		Int_t	fMaxMultForMultBinning;   
		Int_t	fNMultBins;
		Bool_t fAnalyseAddedSignals;  //need to set fFillMCHistos=true as well!  (if(!fFillMCHistos) fAnalyseAddedSignals = kFALSE;)
		Bool_t fFillClusterHitmapsAddedSignals;
		Bool_t fFillDecayPhotonInfoAddedSig;
		Bool_t fAnalyseMCPi0;  //
		Bool_t fAnalyseMCEta;  // 
		Bool_t fFillMCHistos;	//needs to be set to true to analyse MC!
		Int_t fIHijingMin;	// first hijing particle index      (will remain 0)
		Int_t fIHijingMax;	// last hijing particle index
		Int_t fIPi0Min;   	// first added pi0 (particle index)
		Int_t fIPi0Max; 		// last  added pi0 (particle index)
		Int_t fIEtaMin;      // first added eta (particle index)
		Int_t fIEtaMax; 	   // last  added eta (particle index)
		Int_t fIPi0EMCMin;	// first particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
		Int_t fIPi0EMCMax;	// last  particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
		Int_t fIPi0PHOSMin;	// first particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
		Int_t fIPi0PHOSMax;	// last  particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
		Int_t fIEtaEMCMin;	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)
		Int_t fIEtaEMCMax;	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)	
		Int_t fIEtaPHOSMin;	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
		Int_t fIEtaPHOSMax;	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
		Int_t fNAddedSignalsTotal;	
		Double_t fTsalisPi0Param1;  //parameters used in this::CalculateAddedSignalWeight
		Double_t fTsalisPi0Param2;
		Double_t fTsalisPi0Param3;
		Double_t fTsalisPi0Param4;
		Double_t fTsalisPi0Param5;
		Double_t fExpParam1;
		Double_t fExpParam2;
		Double_t fExpParam3;
		Double_t fExpParam4;
		Bool_t fFillNDaughtersVsPtAndR;
		Bool_t fFillMPtForSingleOrMultContrClus;
		TString fAdditionalFileNameString;  //String appears at the end of the filname. (should start with "_")
		// ********** END Cut Variables ********** //

		Int_t			fEventCounter;         // number of analyzed events
		TList			*fOutput;       	 //! Output list
		Bool_t		fMcMode;        	// monte carlo mode
		AliVEvent	*fAnyEv;	 	//!pointer to input event

		TH1F			*fH1NEvents;     	//! # of good Events in Bin 3 (for old peakExtraction codes) 
		TH1F			*fH1NEventsNamed; 	//! # of Events in different classes! 
		TH1F			*fH1NClusters;  		//!  # of clusters/evt (EMCal and PHOS)
		TH1F			*fH1NClustersPHOS;	//!  # of PHOS clusters/evt
		TH1F			*fH1NClustersPHOSafterCuts;	//!  # of PHOS clusters/evt that are not on a bad cell (after all cuts)
		TH2F			*fH2NCellsPerClusterVsClusterEnergy;					//! #cells per cluster
		TH1F			*fH1Zvtx;  		//!  # of clusters/evt
		TH1F			*fH1Mass;        		//!  Mass spectrum same events
		TH1F			*fH1MassMixed;    		//!  Mass spectrum mixed events
		TH1F			*fH1ClusterE;        		//!  cluster energy spectrum
		TH1F			*fH1ClusterEAfterCuts;        //!  cluster energy spectrum after cuts

		// MC-Histogramms: 
		
		// All Pions (in all directions) in the different Categories: 
		TH1F			*fH1Pi0DecayModes; 	//! number pi0s vs different decay modes / reaction modes
		TH1F			*fH1Pi0DecayModesAfterCuts; 	//! number pi0s vs different decay modes / reaction modes
		TH1F			*fH1Pi0DecayModesIsPrimary; 	//! number pi0s vs different decay modes / reaction modes
		TH1F			*fH1Pi0DecayModesAddedPi0; 	//! number pi0s vs different decay modes / reaction modes for addedPi0PHOS from clusters
		TH1F			*fH1Pi0DecayModesAddedPi0PHOS; 	//! number pi0s vs different decay modes / reaction modes for addedPi0PHOS from clusters
		TH2F			*fH2NDaughtersVsPtSimPi0s; 	//! number of daughters vs pT of the simulated pi0s
		TH2F			*fH2NDaughtersVsdRSimPi0s; 	//! number of daughters vs vtx dist to evt vtx of the simulated pi0s
		TH1F			*fH1MCpionVertDistToEventVert; 	//! Distance of Pi0-vertex to event-vertex
		TH1F			*fH1MCpionVertDistToEventVertLowDR; 	//! Distance of Pi0-vertex to event-vertex
		TH1F			*fH1MCpionVertDistToEventVertIsPrimary; 	//! Distance of Pi0-vertex to event-vertex
		TH1F			*fH1Pi0TruthPt;        		//!  Pt spectrum of all Pi0 from MC! 
		TH1F			*fH1K0Pi0TruthPt;        		//!  Pt spectrum of all Pi0 from K0-decays from MC! 
		TH1F			*fH1PriPi0TruthPt;        	//!  Pt spectrum of all primary Pi0 from MC!
		TH1F		 	*fH1SecPi0TruthPt; 				//!  Pt spectrum of all secondary Pi0 exept K0 from MC!
		TH1F			*fH1K0Pi0TruthPtPhi2PiY05;        	//!  Pt spectrum of all Pi0 from K0-decays from MC! 
		TH1F			*fH1PriPi0TruthPtPhi2PiY05;        //!  Pt spectrum of all primary Pi0 from MC!
		TH1F			*fH1SecPi0TruthPtPhi2PiY05; 			//!  Pt spectrum of all secondary Pi0 exept K0 from MC!
		
		// Pions that fly in PHOS-direction in the different Categories
		TH1F			*fH1Pi0TruthPtPhos;        	//!  Pt spectrum of all Pi0 pointed at PHOS from MC! 
		TH1F			*fH1K0Pi0TruthPtPhos;        	//!  Pt spectrum of all Pi0 from K0-decays pointed at PHOS from MC! 
		TH1F			*fH1PriPi0TruthPtPhos;       	//!  Pt spectrum of all primary Pi0 pointed at PHOS from MC!
		TH1F			*fH1SecPi0TruthPtPhos; 		//!  Pt spectrum of all secondary Pi0 exept K0 pointed at PHOS from MC! 

		// All Pions in y < something
		TH1F			*fH1Pi0TruthPtPhi2PiY05;	//!  Pt spectrum from MC
		TH1F			*fH1Pi0TruthPtPhi2PiY05dR1;	//!  Pt spectrum from MC
		TH1F			*fH1PriPi0TruthPtPhi2PiY05dR1;	//!  Pt spectrum from MC
		TH1F			*fH1Pi0TruthPtPhi2piY03;	//!  Pt spectrum from MC
		TH2F			*fH2Pi0TruthPhiEta;    		//!  etaphi spectrum from MC! 

		// Pions, where both Photons land in PHOS
		TH1F			*fH1ElectronConversionR;					//!	Distance of conversionpoint to beamline
		TH1F     	*fH1Pi0TruthPtPhotonsPhos;      //!  Pt spectrum of Pions with both Photons in PHOS from MC! 
		TH1F			*fH1K0Pi0TruthPtPhotonsPhos; 	//!  Pt spectrum of Pions from K0 with both Photons in PHOS from MC!
		TH1F			*fH1PriPi0TruthPtPhotonsPhos; 	//!  Pt spectrum of primary Pi0 with both Photons in PHOS from MC!
		TH1F			*fH1SecPi0TruthPtPhotonsPhos; 	//!  Pt spectrum of secondary Pi0 exept K0 with both Photons in PHOS from MC!
		
		// Added signal histos
		TH2F			*fH2YVsPhiGenPi0;		   	//! EtaVsPhi distribution of pi0s in minBias Generator
		TH2F			*fH2YVsPhiGenEta;	   		//! EtaVsPhi distribution of etas in hijing
		TH2F			*fH2YVsPhiAddedPi0;		   //! EtaVsPhi distribution of pi0s in pi0_1
		TH2F			*fH2YVsPhiAddedEta;	   	//! EtaVsPhi distribution of etas in eta_2
		TH2F			*fH2YVsPhiAddedPi0PHOS;	   //! EtaVsPhi distribution of pi0s in pi0PHS_4
		TH2F			*fH2YVsPhiAddedEtaPHOS;	   //! EtaVsPhi distribution of etas in etaPHS_6
		TH1F			*fH1YAddedPi0;					//!
		TH1F			*fH1EtaAddedPi0;				//!
		TH1F			*fH1YAddedPi0PHOS;			//!
		TH1F			*fH1EtaAddedPi0PHOS;			//!
		TH1F			*fH1ClusterEAddedPi0;		//!
		TH1F			*fH1ClusterEAddedPi0PHOS;	//!
		
		//decayphoton info for added signals 
		TH1F			*fH1DecGammAddPi0Eta;					//!
		TH1F			*fH1DecGammAddPi0Y;						//!
		TH1F			*fH1DecGammAddPi0Phi;					//!
		TH1F			*fH1DecGammAddPi0E;						//!
		TH1F			*fH1DecGammInPHOSAddPi0E;				//!
		TH1F			*fH1DecGammAddPi0Asymm;					//!
		TH1F			*fH1DecGammAddPi0OpAngle;				//!
		TH1F			*fH1DecGammAddPi0ConvR;					//!
		TH1F			*fH1DecGammAddPi0ConvRate;				//!
		TH1F			*fH1TruthPtAddedPi0GammasInPHOS;		//!
		TH1F			*fH1DecGammAddPi0PHOSEta;				//!
		TH1F			*fH1DecGammAddPi0PHOSY;					//!
		TH1F			*fH1DecGammAddPi0PHOSPhi;				//!
		TH1F			*fH1DecGammAddPi0PHOSE;					//!
		TH1F			*fH1DecGammInPHOSAddPi0PHOSE;			//!
		TH1F			*fH1DecGammAddPi0PHOSAsymm;			//!
		TH1F			*fH1DecGammAddPi0PHOSOpAngle;			//!
		TH1F			*fH1DecGammAddPi0PHOSConvR;			//!
		TH1F			*fH1DecGammAddPi0PHOSConvRate;		//!
		TH1F			*fH1TruthPtAddedPi0PHOSGammasInPHOS; //!
		
		//truth pt for added signals
		TH1F			*fH1TruthPtGenPi0;			//! pT distribution of pi0s in minBias Generator
		TH1F			*fH1TruthPtGenEta;			//! pT distribution of etas in eta_2
		TH1F			*fH1TruthPtAddedPi0;			//! pT distribution of pi0s in pi0_1
		TH1F			*fH1TruthPtAddedPi02PiY05;			//! pT distribution of pi0s in pi0_1
		TH1F			*fH1TruthPtAddedPi0MesonInPHOS;			//! pT distribution of pi0s in pi0_1
		TH1F			*fH1TruthPtAddedPi0PHOSMesonInPHOS;			//! pT distribution of pi0s in pi0PHOS_4
		TH1F			*fH1TruthPtAddedEta;			//! pT distribution of etas in eta_2
		TH1F			*fH1TruthPtAddedPi0PHOS;	//! pT distribution of pi0s in pi0PHS_4
		TH1F			*fH1TruthPtAddedEtaPHOS;	//! pT distribution of etas in etaPHS_6
		TH2F			*fH2MPtAddedPi0;				//! Mass vs pT histo for reconstructed added Pi0 (shot in wide range)
		TH2F			*fH2MPtAddedPi0_unweighed;				//! unweighed Mass vs pT histo for reconstructed added Pi0 (shot in wide range)
		TH2F			*fH2MPtAddedPi0PHOS;			//! Mass vs pT histo for reconstructed added Pi0 (shot at PHOS directly)
		TH2F			*fH2MPtAddedPi0PHOS_unweighed;			//! unweighed Mass vs pT histo for reconstructed added Pi0 (shot at PHOS directly)
		TH2F			*fH2MPtAddedPi0PHOSSingleContr;			//! Mass vs pT histo for reconstructed added Pi0, both cluster single contributer(shot at PHOS directly)
		TH2F			*fH2MPtAddedPi0PHOSMultContr;			//! Mass vs pT histo for reconstructed added Pi0, both cluster mult contributer (shot at PHOS directly)
		TH2F			*fH2MEnergyDiffAddedPi0PHOS;			//! Mass vs max(abs(ClusterE-PartcleE)) histo for reconstructed added Pi0 (shot at PHOS directly)
		TH2F			*fH2PtRecVsPtTruthAddedPi0PHOS;					//! Reconstructed pT vs true pT of added pi0
		TH2F			*fH2MPtAddedEta;				//! Mass vs pT histo for reconstructed added Eta (shot in wide range)
		TH2F			*fH2MPtAddedEtaPHOS;			//! Mass vs pT histo for reconstructed added Eta (shot at PHOS directly)
		TH2F			*fH2ClusterEVSPhotonE;		//! ClusterE vs PhotonE for monte carlo clusters with photons as main contributer

		//
		TH2D			*fH2HitmapEtaVsPhi;		//! Hitmap
		TH1F			*fH1Chi2;       		//!  pseudorapidity spectrum
		TH1F			*fH1NTrkMatch;       		//!  number of matched tracks to a cluster
		TH1F			*fH1ClusterDisp;       		//!  cluster dispersion
		TH2F			*fH2Ellipse;       		//!  ellipse axis?

		TH3F			*fH3MPtAsymm;       		//!  3dimensional M vs pt vs asym-cut
		TH3F		 	*fH3MPtModules;		   //!  3dimensional M vs pt vs module-combination
		TH3F			*fH3MPtTiming; 			//! 3dimensional M vs pt vs timingCutCombination
		TH3F			*fH3MPtAsymmMix;       	//!  2dimensional E vs mom
		TH3F			*fH3MPtModulesMix;		//!  3dimensional M vs pt vs module-combination
		TH3F			*fH3MPtTimingMix; 			//! 3dimensional M vs pt vs timingCutCombination
		TH2F			*fH2DphiDeta;		        //!  2dimensional E vs mom
		TH2F			*fH2DphiDetaMix;       //!  2dimensional E vs mom
		TH2F			*fH2CellsM02;       //!  
		TH1F			*fH1ClusterM02; 	//!

		TH1F			*fH1DistPileUpPrimVert; //! Dist prim-pile-up vtx
		TH1F			*fH1NPrimVertContribut; //! Number of Contributors to primVtx
		TH1F			*fH1nSPDPileUpVtxs; //! Number of SPD PileUpVertices
		TH1F			*fH1NClustersVsCuts; //! Number of RejectedClusters
		TProfile		*fTProfMeanClusterEnergyVsCuts; //! meanClusterEnergyAfterCuts
		TH2F			*fH2EAfterCutsVsModNum; //! ClusterE vs ModNum
		TH2F			*fH2CellIdVsE; //! CellID vs E
		TH2F			*fH2LocalMaxCellsIdVsE; //! CellID vs E for cells that are local maximum of a cluster
		TH2F			*fH2ClusterPosCellsIdVsE;  //! CellID vs E for cells that are the Position of a cluster
		TH1F			*fH1ClusterTOFWeightedWithE;	//!
		TH2F			*fH2ClusterTOFVsE; 				//! Cluster-Time against Cluster-Energy
		TH1F			*fH1CellMCLabel; //!
		TH1F			*fH1AppliedClusterCuts; //!
		TH1F			*fH1DistanceToBadChannel;   //!  Filled with cluster->GetDist
		TH2F			*fH2ClusterPositionsMod1; 	//! Hitmap Module 1
		TH2F	  		*fH2ClusterPositionsMod2; 	//! Hitmap Module 2
		TH2F	  		*fH2ClusterPositionsMod3; 	//! Hitmap Module 3
		TH2F			*fH2ClusterPosMod1Gen1;		//!
		TH2F			*fH2ClusterPosMod3Gen1;		//!
		TH2F			*fH2ClusterPosMod1Gen2;		//!
		TH2F			*fH2ClusterPosMod3Gen2;		//!

		TNtuple		*fNTupelClusterEnergyMod1;		//! NTuple to make modulewise calibration 
		TNtuple		*fNTupelClusterEnergyMod2; 	//! NTuple to make modulewise calibration 
		TNtuple		*fNTupelClusterEnergyMod3;		//! NTuple to make modulewise calibration 

		AliPHOSGeometry	*fPHOSGeo;   //!PHOS geometry
		AliPHOSCalibData	*fPHOSCalibData; //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
		AliAnalysisUtils	*fUtils;     //utils for zvtxcut

		std::vector<TLorentzVector> fPhotons[fgkPoolDepth][fgkZvtxBins][fgkMultBins];	  //3D array of vectors of TLorentzVectors.
		std::vector<Int_t> fModNumber[fgkPoolDepth][fgkZvtxBins][fgkMultBins];
		std::vector<Bool_t> fClusterTiming[fgkPoolDepth][fgkZvtxBins][fgkMultBins];
		
		std::vector<TLorentzVector> fPhotonsAdded[4];	   	//array of vectors of TLorentzVectors. Used to store photons of 4 classes of added signals (lhc13e7)
		std::vector<Bool_t> fClusterHasSingleContr[4];	   	//array of vectors of bools.  Used to store info (whether its has only one contributer) about clusters of 4 classes of added signals (lhc13e7)
		std::vector<Float_t> fDifferenceClusterParticleE[4];	//array of vectors of floats. Used to store info (ClusterE - ParticleE) about clusters of 4 classes of added signals (lhc13e7)
		Double_t fAddedSignalWeight[4]; 								//array of doubles. Used to store weights of added signals for 4 classes (lhc13e7)
		Double_t fTruthPtAddedPi0PHOS;

		AliAnalysisTaskPHOSNeutralMeson(const AliAnalysisTaskPHOSNeutralMeson&); // not implemented
		AliAnalysisTaskPHOSNeutralMeson& operator=(const AliAnalysisTaskPHOSNeutralMeson&); // not implemented

		ClassDef(AliAnalysisTaskPHOSNeutralMeson, 1);
};
#endif
