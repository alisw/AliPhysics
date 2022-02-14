/*************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-------------------------------------------------------------------------
////     AnalyisTask for neutral meson analysis with PHOS 
////     Runs on ESDs and AODs, Tested on pp and pPb data.
////     Authors: Malte Hecker, Fabian Pliquett
////     Date: 01/12/2015
////-------------------------------------------------------------------------

#include "AliAnalysisTaskPHOSNeutralMeson.h"

#include <vector>
#include <Riostream.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TProfile.h>
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliExternalTrackParam.h"
#include "AliOADBContainer.h"
#include "AliPHOSCalibData.h"
#include "AliAODMCParticle.h"


//PHOS Additional
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliPHOSGeometry.h"
#include "AliESDCaloCluster.h"
#include "TGeoManager.h"

// ROOT includes
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TH2F.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

// STEER includes
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"
#include "AliTrackerBase.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPHOSNeutralMeson)
	//________________________________________________________________________
	AliAnalysisTaskPHOSNeutralMeson::AliAnalysisTaskPHOSNeutralMeson() : 
	AliAnalysisTaskSE(),

	// ********** START Cut Variables ********** //
	fClusterMinCells(0),
	fClusterMinE(0),
	fClusterMinM02(0), 
	fDoDistToBadCellCutOnCellLevel(0),
	fDistToBadCellOnCellLevel(0),
	fDoDistToBadCellCut(0),
	fDistToBadCell(0.0),
	fDoTimingCut(0),
	fTimingCutMin(0.0),
	fTimingCutMax(0.0),
	fEtaAccMin(0.0),
	fEtaAccMax(0.0),
	fPhiAccMin(0.0),
	fPhiAccMax(0.0),
	fFillCellIdVsE(0),
	fFillTimingHistos(0),
	fUseOnlyCINT1Events(0),
	fFillClusterHitmaps(0),
	fFillFeedDownHistos(0), 
	fDoClusterEnergyRecalibration(0), 
	fRecalibrationOption("noNonLinCor"),
	fDoZvertexCut(0),
	fZvertexCut(10.0),
	fFillHMassPtModules(0),
	fFillHMassPtTiming(0),
	fFillNewAsymmClasses(0),
	fAsymmClass1(0.0),
	fAsymmClass2(0.0),
	fAsymmClass3(0.0),
	fFillNTupelClusterE(0),
	fRecalibrateModuleWise(0),
	fRecalFactorMod1(1.0),
	fRecalFactorMod2(1.0),
	fRecalFactorMod3(1.0),
	fDoPeakSmearing(0),
	fSmearFactor(0.0001), 
	fSelectedCollisionCandidates(""),
	fSelectedPhysSelection(""),
	fSelectedTender(""),
	fOutputFileName(""),
	fMPtHistoMode("nrmlMptHst"),
	fUsedBadMap(""),
	fUseIsVertexSelected2013pA(0),
	fApplyBadMapManually(kFALSE),
	fMaxMultForMultBinning(60),   //default makes sense for pPb. Set differently in your AddTask if needed.
	fNMultBins(30),
	fAnalyseAddedSignals(kFALSE),    //fFillMCHistos needs to be set true as well!
	fFillClusterHitmapsAddedSignals(kFALSE),
	fFillDecayPhotonInfoAddedSig(kFALSE),
	fAnalyseMCPi0(kFALSE),    //fFillMCHistos and fAnalyseAddedSignals needs to be set true as well!
	fAnalyseMCEta(kFALSE),    //fFillMCHistos and fAnalyseAddedSignals needs to be set true as well!
	fFillMCHistos(kFALSE),				//needs to be set to true to analyse MC!
	fIHijingMin(-1),  
	fIHijingMax(-1),	// last hijing particle index                     
	fIPi0Min(-1),   	   // first added pi0 (particle index)
	fIPi0Max(-1), 		// last  added pi0 (particle index)
	fIEtaMin(-1),      // first added eta (particle index)
	fIEtaMax(-1), 	   // last  added eta (particle index)
	fIPi0EMCMin(-1),	// first particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0EMCMax(-1),	// last  particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0PHOSMin(-1),	// first particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0PHOSMax(-1),	// last  particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIEtaEMCMin(-1),	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)
	fIEtaEMCMax(-1),	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)	
	fIEtaPHOSMin(-1),	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
	fIEtaPHOSMax(-1),	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
	fNAddedSignalsTotal(-1),
	fTsalisPi0Param1(0.0), //parameters used in this::CalculateAddedSignalWeight
	fTsalisPi0Param2(0.0),
	fTsalisPi0Param3(0.0),
	fTsalisPi0Param4(0.0),
	fTsalisPi0Param5(0.0),
	fExpParam1(0.0),
	fExpParam2(0.0),
	fExpParam3(0.0),
	fExpParam4(0.0),
	fFillNDaughtersVsPtAndR(0),
	fFillMPtForSingleOrMultContrClus(0),
	fAdditionalFileNameString(""),
	// ********** END Cut Variables ********** //

	fEventCounter(0),
	fOutput(0),
	fMcMode(0),
	fAnyEv(0),
	fH1NEvents(0),
	fH1NEventsNamed(0),
	fH1NClusters(0), 
	fH1NClustersPHOS(0), 
	fH1NClustersPHOSafterCuts(0), 
	fH2NCellsPerClusterVsClusterEnergy(0),
	fH1Zvtx(0), 
	fH1Mass(0), 
	fH1MassMixed(0), 
	fH1ClusterE(0),  
	fH1ClusterEAfterCuts(0),
	fH1Pi0DecayModes(0),
	fH1Pi0DecayModesAfterCuts(0),
	fH1Pi0DecayModesIsPrimary(0),
	fH1Pi0DecayModesAddedPi0(0),
	fH1Pi0DecayModesAddedPi0PHOS(0),
	fH2NDaughtersVsPtSimPi0s(0),
	fH2NDaughtersVsdRSimPi0s(0),
	fH1MCpionVertDistToEventVert(0),
	fH1MCpionVertDistToEventVertLowDR(0),
	fH1MCpionVertDistToEventVertIsPrimary(0),
	fH1Pi0TruthPt(0), 
	fH1K0Pi0TruthPt(0), 
	fH1PriPi0TruthPt(0), 
	fH1SecPi0TruthPt(0), 
	fH1K0Pi0TruthPtPhi2PiY05(0), 
	fH1PriPi0TruthPtPhi2PiY05(0), 
	fH1SecPi0TruthPtPhi2PiY05(0), 
	fH1Pi0TruthPtPhos(0), 
	fH1K0Pi0TruthPtPhos(0), 
	fH1PriPi0TruthPtPhos(0),
	fH1SecPi0TruthPtPhos(0), 
	fH1Pi0TruthPtPhi2PiY05(0),
	fH1Pi0TruthPtPhi2PiY05dR1(0),
	fH1PriPi0TruthPtPhi2PiY05dR1(0),
	fH1Pi0TruthPtPhi2piY03(0),
	fH2Pi0TruthPhiEta(0), 
	fH1ElectronConversionR(0),
	fH1Pi0TruthPtPhotonsPhos(0), 
	fH1K0Pi0TruthPtPhotonsPhos(0), 
	fH1PriPi0TruthPtPhotonsPhos(0),
	fH1SecPi0TruthPtPhotonsPhos(0),  
	fH2YVsPhiGenPi0(0),
	fH2YVsPhiGenEta(0),
	fH2YVsPhiAddedPi0(0),
	fH2YVsPhiAddedEta(0),
	fH2YVsPhiAddedPi0PHOS(0),
	fH2YVsPhiAddedEtaPHOS(0),
	fH1YAddedPi0(0),
	fH1EtaAddedPi0(0),
	fH1YAddedPi0PHOS(0),
	fH1EtaAddedPi0PHOS(0),
	fH1ClusterEAddedPi0(0),
	fH1ClusterEAddedPi0PHOS(0),
	fH1DecGammAddPi0Eta(0),
	fH1DecGammAddPi0Y(0),
	fH1DecGammAddPi0Phi(0),
	fH1DecGammAddPi0E(0),
	fH1DecGammInPHOSAddPi0E(0),
	fH1DecGammAddPi0Asymm(0),
	fH1DecGammAddPi0OpAngle(0),
	fH1DecGammAddPi0ConvR(0),
	fH1DecGammAddPi0ConvRate(0),
	fH1TruthPtAddedPi0GammasInPHOS(0),
	fH1DecGammAddPi0PHOSEta(0),
	fH1DecGammAddPi0PHOSY(0),
	fH1DecGammAddPi0PHOSPhi(0),
	fH1DecGammAddPi0PHOSE(0),
	fH1DecGammInPHOSAddPi0PHOSE(0),
	fH1DecGammAddPi0PHOSAsymm(0),
	fH1DecGammAddPi0PHOSOpAngle(0),
	fH1DecGammAddPi0PHOSConvR(0),
	fH1DecGammAddPi0PHOSConvRate(0),
	fH1TruthPtAddedPi0PHOSGammasInPHOS(0),
	fH1TruthPtGenPi0(0),
	fH1TruthPtGenEta(0),
	fH1TruthPtAddedPi0(0),
	fH1TruthPtAddedPi02PiY05(0),
	fH1TruthPtAddedPi0MesonInPHOS(0),
	fH1TruthPtAddedPi0PHOSMesonInPHOS(0),
	fH1TruthPtAddedEta(0),
	fH1TruthPtAddedPi0PHOS(0),
	fH1TruthPtAddedEtaPHOS(0),
	fH2MPtAddedPi0(0),
	fH2MPtAddedPi0_unweighed(0),
	fH2MPtAddedPi0PHOS(0),
	fH2MPtAddedPi0PHOS_unweighed(0),
	fH2MPtAddedPi0PHOSSingleContr(0),
	fH2MPtAddedPi0PHOSMultContr(0),
	fH2MEnergyDiffAddedPi0PHOS(0),
	fH2PtRecVsPtTruthAddedPi0PHOS(0),
	fH2MPtAddedEta(0),
	fH2MPtAddedEtaPHOS(0),
	fH2ClusterEVSPhotonE(0),
	fH2HitmapEtaVsPhi(0),
	fH1Chi2(0),
	fH1NTrkMatch(0),
	fH1ClusterDisp(0),
	fH2Ellipse(0),
	fH3MPtAsymm(0),
	fH3MPtModules(0),
	fH3MPtTiming(0),
	fH3MPtAsymmMix(0), 
	fH3MPtModulesMix(0), 
	fH3MPtTimingMix(0),
	fH2DphiDeta(0), 
	fH2DphiDetaMix(0), 
	fH2CellsM02(0),
	fH1ClusterM02(0),
	fH1NPrimVertContribut(0),
	fH1DistPileUpPrimVert(0),
	fH1nSPDPileUpVtxs(0),
	fH1NClustersVsCuts(0),
	fTProfMeanClusterEnergyVsCuts(0),
	fH2EAfterCutsVsModNum(0),
	fH2CellIdVsE(0),
	fH2LocalMaxCellsIdVsE(0),
	fH2ClusterPosCellsIdVsE(0),
	fH1ClusterTOFWeightedWithE(0),
	fH2ClusterTOFVsE(0),
	fH1CellMCLabel(0),
	fH1AppliedClusterCuts(0),
	fH1DistanceToBadChannel(0),
	fH2ClusterPositionsMod1(0),
	fH2ClusterPositionsMod2(0),
	fH2ClusterPositionsMod3(0),
	fH2ClusterPosMod1Gen1(0),
	fH2ClusterPosMod3Gen1(0),
	fH2ClusterPosMod1Gen2(0),
	fH2ClusterPosMod3Gen2(0),
	fNTupelClusterEnergyMod1(0),
	fNTupelClusterEnergyMod2(0),
	fNTupelClusterEnergyMod3(0),
	fPHOSGeo(0x0),
	fPHOSCalibData(0x0), //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
	fUtils(0),
	fTruthPtAddedPi0PHOS(0.0)
	//fAddedSignalWeight({1.0,1.0,1.0,1.0})
{
    // Dummy constructor ALWAYS needed for I/O.

    // Initialize the PHOS geometry
    // Set bad channel map
    Char_t key[55] ;
    for(Int_t i=0; i<6; i++){
		snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
		fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
	}

	fUtils = new AliAnalysisUtils();                             
}

//________________________________________________________________________
AliAnalysisTaskPHOSNeutralMeson::AliAnalysisTaskPHOSNeutralMeson(const Char_t *name) :
    AliAnalysisTaskSE(name),

	// ********** START Cut Variables ********** //
	fClusterMinCells(0),
	fClusterMinE(0),
	fClusterMinM02(0), 
	fDoDistToBadCellCutOnCellLevel(0),
	fDistToBadCellOnCellLevel(0),
	fDoDistToBadCellCut(0),
	fDistToBadCell(0.0),
	fDoTimingCut(0),
	fTimingCutMin(0.0),
	fTimingCutMax(0.0),
	fEtaAccMin(0.0),
	fEtaAccMax(0.0),
	fPhiAccMin(0.0),
	fPhiAccMax(0.0),
	fFillCellIdVsE(0),
	fFillTimingHistos(0),
	fUseOnlyCINT1Events(0),
	fFillClusterHitmaps(0),
	fFillFeedDownHistos(0), 
	fDoClusterEnergyRecalibration(0), 
	fRecalibrationOption("noNonLinCor"),
	fDoZvertexCut(0),
	fZvertexCut(10.0),
	fFillHMassPtModules(0),
	fFillHMassPtTiming(0),
	fFillNewAsymmClasses(0),
	fAsymmClass1(0.0),
	fAsymmClass2(0.0),
	fAsymmClass3(0.0),
	fFillNTupelClusterE(0),
	fRecalibrateModuleWise(0),
	fRecalFactorMod1(1.0),
	fRecalFactorMod2(1.0),
	fRecalFactorMod3(1.0),	
	fDoPeakSmearing(0),
	fSmearFactor(0.0001), 
	fSelectedCollisionCandidates(""),
	fSelectedPhysSelection(""),
	fSelectedTender(""),
	fOutputFileName(""),
	fMPtHistoMode("nrmlMptHst"),
	fUsedBadMap(""),
	fUseIsVertexSelected2013pA(0),
	fApplyBadMapManually(false),
	fMaxMultForMultBinning(60),   //default makes sense for pPb. Set differently in your AddTask if needed.
	fNMultBins(30),
	fAnalyseAddedSignals(kFALSE),    //fFillMCHistos needs to be set true as well!
	fFillClusterHitmapsAddedSignals(kFALSE),
	fFillDecayPhotonInfoAddedSig(kFALSE),
	fAnalyseMCPi0(kFALSE),    //fFillMCHistos and fAnalyseAddedSignals needs to be set true as well!
	fAnalyseMCEta(kFALSE),    //fFillMCHistos and fAnalyseAddedSignals needs to be set true as well!
	fFillMCHistos(kFALSE),				//needs to be set to true to analyse MC!
	fIHijingMin(-1),
	fIHijingMax(-1),	// last hijing particle index
	fIPi0Min(-1),   	   // first added pi0 (particle index)
	fIPi0Max(-1), 		// last  added pi0 (particle index)
	fIEtaMin(-1),      // first added eta (particle index)
	fIEtaMax(-1), 	   // last  added eta (particle index)
	fIPi0EMCMin(-1),	// first particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0EMCMax(-1),	// last  particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0PHOSMin(-1),	// first particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIPi0PHOSMax(-1),	// last  particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
	fIEtaEMCMin(-1),	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)
	fIEtaEMCMax(-1),	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)	
	fIEtaPHOSMin(-1),	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
	fIEtaPHOSMax(-1),	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
	fNAddedSignalsTotal(-1),
	fTsalisPi0Param1(0.0),
	fTsalisPi0Param2(0.0),
	fTsalisPi0Param3(0.0),
	fTsalisPi0Param4(0.0),
	fTsalisPi0Param5(0.0),
	fExpParam1(0.0),
	fExpParam2(0.0),
	fExpParam3(0.0),
	fExpParam4(0.0),
	fFillNDaughtersVsPtAndR(0),
	fFillMPtForSingleOrMultContrClus(0),
	fAdditionalFileNameString(""),
	// ********** END Cut Variables ********** //

	fEventCounter(0),
	fOutput(0),
	fMcMode(0),
	fAnyEv(0),
	fH1NEvents(0),
	fH1NEventsNamed(0),
	fH1NClusters(0), 
	fH1NClustersPHOS(0), 
	fH1NClustersPHOSafterCuts(0), 
	fH2NCellsPerClusterVsClusterEnergy(0),
	fH1Zvtx(0), 
	fH1Mass(0), 
	fH1MassMixed(0), 
	fH1ClusterE(0),  
	fH1ClusterEAfterCuts(0),
	fH1Pi0DecayModes(0),
	fH1Pi0DecayModesAfterCuts(0),
	fH1Pi0DecayModesIsPrimary(0),
	fH1Pi0DecayModesAddedPi0(0),
	fH1Pi0DecayModesAddedPi0PHOS(0),
	fH2NDaughtersVsPtSimPi0s(0),
	fH2NDaughtersVsdRSimPi0s(0),
	fH1MCpionVertDistToEventVert(0),
	fH1MCpionVertDistToEventVertLowDR(0),
	fH1MCpionVertDistToEventVertIsPrimary(0),
	fH1Pi0TruthPt(0),
	fH1K0Pi0TruthPt(0),
	fH1PriPi0TruthPt(0),
	fH1SecPi0TruthPt(0),     
	fH1K0Pi0TruthPtPhi2PiY05(0),    
	fH1PriPi0TruthPtPhi2PiY05(0), 
	fH1SecPi0TruthPtPhi2PiY05(0), 
	fH1Pi0TruthPtPhos(0), 
	fH1K0Pi0TruthPtPhos(0), 
	fH1PriPi0TruthPtPhos(0),
	fH1SecPi0TruthPtPhos(0), 
	fH1Pi0TruthPtPhi2PiY05(0),
	fH1Pi0TruthPtPhi2PiY05dR1(0),
	fH1PriPi0TruthPtPhi2PiY05dR1(0),
	fH1Pi0TruthPtPhi2piY03(0),
	fH2Pi0TruthPhiEta(0), 
	fH1ElectronConversionR(0),
	fH1Pi0TruthPtPhotonsPhos(0),
	fH1K0Pi0TruthPtPhotonsPhos(0),
	fH1PriPi0TruthPtPhotonsPhos(0),
	fH1SecPi0TruthPtPhotonsPhos(0), 
	fH2YVsPhiGenPi0(0),
	fH2YVsPhiGenEta(0),
	fH2YVsPhiAddedPi0(0),
	fH2YVsPhiAddedEta(0),
	fH2YVsPhiAddedPi0PHOS(0),
	fH2YVsPhiAddedEtaPHOS(0),
	fH1YAddedPi0(0),
	fH1EtaAddedPi0(0),
	fH1YAddedPi0PHOS(0),
	fH1EtaAddedPi0PHOS(0),
	fH1ClusterEAddedPi0(0),
	fH1ClusterEAddedPi0PHOS(0),
	fH1DecGammAddPi0Eta(0),
	fH1DecGammAddPi0Y(0),
	fH1DecGammAddPi0Phi(0),
	fH1DecGammAddPi0E(0),
	fH1DecGammInPHOSAddPi0E(0),
	fH1DecGammAddPi0Asymm(0),
	fH1DecGammAddPi0OpAngle(0),
	fH1DecGammAddPi0ConvR(0),
	fH1DecGammAddPi0ConvRate(0),
	fH1TruthPtAddedPi0GammasInPHOS(0),
	fH1DecGammAddPi0PHOSEta(0),
	fH1DecGammAddPi0PHOSY(0),
	fH1DecGammAddPi0PHOSPhi(0),
	fH1DecGammAddPi0PHOSE(0),
	fH1DecGammInPHOSAddPi0PHOSE(0),
	fH1DecGammAddPi0PHOSAsymm(0),
	fH1DecGammAddPi0PHOSOpAngle(0),
	fH1DecGammAddPi0PHOSConvR(0),
	fH1DecGammAddPi0PHOSConvRate(0),
	fH1TruthPtAddedPi0PHOSGammasInPHOS(0),
	fH1TruthPtGenPi0(0),
	fH1TruthPtGenEta(0),
	fH1TruthPtAddedPi0(0),
	fH1TruthPtAddedPi02PiY05(0),
	fH1TruthPtAddedPi0MesonInPHOS(0),
	fH1TruthPtAddedPi0PHOSMesonInPHOS(0),
	fH1TruthPtAddedEta(0),
	fH1TruthPtAddedPi0PHOS(0),
	fH1TruthPtAddedEtaPHOS(0),
	fH2MPtAddedPi0(0),
	fH2MPtAddedPi0_unweighed(0),
	fH2MPtAddedPi0PHOS(0),
	fH2MPtAddedPi0PHOS_unweighed(0),
	fH2MPtAddedPi0PHOSSingleContr(0),
	fH2MPtAddedPi0PHOSMultContr(0),
	fH2MEnergyDiffAddedPi0PHOS(0),
	fH2PtRecVsPtTruthAddedPi0PHOS(0),
	fH2MPtAddedEta(0),
	fH2MPtAddedEtaPHOS(0),
	fH2ClusterEVSPhotonE(0),
	fH2HitmapEtaVsPhi(0),
	fH1Chi2(0),
	fH1NTrkMatch(0),
	fH1ClusterDisp(0),
	fH2Ellipse(0),
	fH3MPtAsymm(0), 
	fH3MPtModules(0), 
	fH3MPtTiming(0),
	fH3MPtAsymmMix(0), 
	fH3MPtModulesMix(0),
	fH3MPtTimingMix(0),
	fH2DphiDeta(0), 
	fH2DphiDetaMix(0), 
	fH2CellsM02(0),
	fH1ClusterM02(0),
	fH1NPrimVertContribut(0),
	fH1DistPileUpPrimVert(0),
	fH1nSPDPileUpVtxs(0),
	fH1NClustersVsCuts(0),
	fTProfMeanClusterEnergyVsCuts(0),
	fH2EAfterCutsVsModNum(0),
	fH2CellIdVsE(0),
	fH2LocalMaxCellsIdVsE(0),
	fH2ClusterPosCellsIdVsE(0),
	fH1ClusterTOFWeightedWithE(0),
	fH2ClusterTOFVsE(0), 
	fH1CellMCLabel(0),
	fH1AppliedClusterCuts(0),
	fH1DistanceToBadChannel(0),
	fH2ClusterPositionsMod1(0),
	fH2ClusterPositionsMod2(0),
	fH2ClusterPositionsMod3(0),
	fH2ClusterPosMod1Gen1(0),
	fH2ClusterPosMod3Gen1(0),
	fH2ClusterPosMod1Gen2(0),
	fH2ClusterPosMod3Gen2(0),
	fNTupelClusterEnergyMod1(0),
	fNTupelClusterEnergyMod2(0),
	fNTupelClusterEnergyMod3(0),
	fPHOSGeo(0x0),
	fPHOSCalibData(0x0), //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
	fUtils(0),
	fTruthPtAddedPi0PHOS(0.0)
	//fAddedSignalWeight({1.0,1.0,1.0,1.0})
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());     // for output list

    // Set bad channel map
    Char_t key[55];
    for(Int_t i=0; i<6; i++){
		snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
		fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
    }
   
	fUtils = new AliAnalysisUtils();                             
}

//________________________________________________________________________
AliAnalysisTaskPHOSNeutralMeson::~AliAnalysisTaskPHOSNeutralMeson()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
	delete fOutput;
    }
}

//________________________________________________________________________
void AliAnalysisTaskPHOSNeutralMeson::UserCreateOutputObjects()
{
	// Create histograms
	// Called once (on the worker node)
	fOutput = new TList();
	fOutput->SetOwner();  // IMPORTANT!
	
	if(!fFillMCHistos) fAnalyseAddedSignals = kFALSE;

	cout << "__________AliAnalysisTaskFraNeutralMeson: Input settings__________" << endl;
	cout << " fFillMCHistos:        " << fFillMCHistos   << endl;
	cout << " fAnalyseAddedSignals: " << fAnalyseAddedSignals   << endl;
	cout << " number of zvtx bins:  " << fgkZvtxBins << endl;
	cout << " number of mult bins:  " << fNMultBins << endl;
	cout << " poolDepth:            " << fgkPoolDepth << endl;
	cout << endl;


	

	Double_t TotalNBins = 0.0;

	// ********** START Create Histograms ********** //
	
	Int_t nEventsbins = 10;
	Float_t nEventslow = -0.5, nEventsup = 9.5;
	fH1NEvents = new TH1F("fH1NEvents", "# of Events in different classes", nEventsbins, nEventslow, nEventsup);
	TotalNBins+= nEventsbins;

	fH1NEventsNamed = new TH1F("fH1NEventsNamed", "# of Events in different classes", 10, 0.5, 10.5);
	fH1NEventsNamed->GetYaxis()->SetTitle("Counts");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(1,"Good_Events");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(2,"All");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(3,"no (FAST && !ALL)");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(4,"prim Vertex");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(5,"pileUpSPD");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(6,"zVertex");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(7,"");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(8,"");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(9,"");
	fH1NEventsNamed->GetXaxis()->SetBinLabel(10,"");
	TotalNBins+=10;

    Int_t nClustersbins = 501;
    Float_t nClusterslow = -0.5, nClustersup = 500.5;
    fH1NClusters = new TH1F("fH1NClusters", "# of clusters per event", nClustersbins, nClusterslow, nClustersup);
    fH1NClusters->GetXaxis()->SetTitle("number of clusters/evt");
    fH1NClusters->GetYaxis()->SetTitle("counts");
    fH1NClusters->SetMarkerStyle(kFullCircle);
    TotalNBins+=nClustersbins;
    
    fH1NClustersPHOS = new TH1F("fH1NClustersPHOS", "# of clusters per event in PHOS (before cuts)", nClustersbins, nClusterslow, nClustersup);
    fH1NClustersPHOS->GetXaxis()->SetTitle("number of clusters/evt (before cuts)");
    fH1NClustersPHOS->GetYaxis()->SetTitle("counts");
    fH1NClustersPHOS->SetMarkerStyle(kFullCircle);
    TotalNBins+=nClustersbins;  
      
    fH1NClustersPHOSafterCuts = new TH1F("fH1NClustersPHOSafterCuts", "# of clusters per event in PHOS (after cuts)", nClustersbins, nClusterslow, nClustersup);
    fH1NClustersPHOSafterCuts->GetXaxis()->SetTitle("number of clusters/evt (after all cuts)");
    fH1NClustersPHOSafterCuts->GetYaxis()->SetTitle("energy");
    fH1NClustersPHOSafterCuts->SetMarkerStyle(kFullCircle);
    TotalNBins+=nClustersbins;
    
    fH2NCellsPerClusterVsClusterEnergy= new TH2F("fH2NCellsPerClusterVsClusterEnergy", "# of cells per  cluster PHOS", 1000, 0, 1000,1000,0,40);
    fH2NCellsPerClusterVsClusterEnergy->GetXaxis()->SetTitle("E (GeV)");
    fH2NCellsPerClusterVsClusterEnergy->GetYaxis()->SetTitle("number of Cells per Cluster");
    fH2NCellsPerClusterVsClusterEnergy->SetMarkerStyle(kFullCircle);
    TotalNBins+=1000;
    
    Int_t nZvertexbins = 501;
    Float_t Zvertexlow = -50.0, Zvertexup = 50.0;
    fH1Zvtx = new TH1F("fH1Zvtx", "# of events", nZvertexbins+500, Zvertexlow-50, Zvertexup+50);
    fH1Zvtx->GetXaxis()->SetTitle("z_{vertex}");
    fH1Zvtx->GetYaxis()->SetTitle("counts");
    fH1Zvtx->SetMarkerStyle(kFullCircle);
    TotalNBins+=nZvertexbins;

    Int_t Mbins = 1000;   
    Float_t Mlow = 0.0, Mup = 1.0;   
    fH1Mass = new TH1F("fH1Mass", "Invariant Mass", Mbins, Mlow, Mup);
    fH1Mass->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    fH1Mass->GetYaxis()->SetTitle("counts");
    fH1Mass->SetMarkerStyle(kFullCircle);
    TotalNBins+=Mbins;

    fH1MassMixed = new TH1F("fH1MassMixed", "Invariant Mass (mixed events)", Mbins, Mlow, Mup);
    fH1MassMixed->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    fH1MassMixed->GetYaxis()->SetTitle("counts");
    fH1MassMixed->SetMarkerStyle(kFullCircle);
    TotalNBins+=Mbins;

    Int_t ptbins = 400;
    Float_t ptlow = 0.0, ptup = 40.0;
    Int_t Ebins = 1000;
    Float_t Elow = 0.0, Eup = 30.0;
    fH1ClusterE = new TH1F("fH1ClusterE", "Cluster Energy in Phos", Ebins, Elow, Eup);
    fH1ClusterE->GetXaxis()->SetTitle("E [GeV]");
    fH1ClusterE->GetYaxis()->SetTitle("counts");
    fH1ClusterE->SetMarkerStyle(kFullCircle);
    TotalNBins+=Ebins;
    
    fH1ClusterEAfterCuts = new TH1F("fH1ClusterEAfterCuts", "Cluster Energy in Phos afterCuts", Ebins, Elow, Eup);
    fH1ClusterEAfterCuts->GetXaxis()->SetTitle("E [GeV]");
    fH1ClusterEAfterCuts->GetYaxis()->SetTitle("counts");
    fH1ClusterEAfterCuts->SetMarkerStyle(kFullCircle);
    TotalNBins+=Ebins;
    
    if(fFillMCHistos) {
    
		 fH1Pi0DecayModes = new TH1F("fH1Pi0DecayModes","fH1Pi0DecayModes",10,0.5,10.5);
		 fH1Pi0DecayModes->GetXaxis()->SetTitle("type of decay / reaction");
		 fH1Pi0DecayModes->GetYaxis()->SetTitle("##pi^{0}");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(1,"0");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(2,"1");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
		 fH1Pi0DecayModes->GetXaxis()->SetBinLabel(7,">=4");
		 TotalNBins+=10;
		 
		 fH1Pi0DecayModesAfterCuts = new TH1F("fH1Pi0DecayModesAfterCuts","fH1Pi0DecayModesAfterCuts",10,0.5,10.5);
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetTitle("type of decay / reaction");
		 fH1Pi0DecayModesAfterCuts->GetYaxis()->SetTitle("##pi^{0}");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(1,"0");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(2,"1");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
		 fH1Pi0DecayModesAfterCuts->GetXaxis()->SetBinLabel(7,">=4");
		 TotalNBins+=10;
		 
		 fH1Pi0DecayModesIsPrimary = new TH1F("fH1Pi0DecayModesIsPrimary","fH1Pi0DecayModesIsPrimary",10,0.5,10.5);
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetTitle("type of decay / reaction");
		 fH1Pi0DecayModesIsPrimary->GetYaxis()->SetTitle("##pi^{0}");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(1,"0");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(2,"1");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
		 fH1Pi0DecayModesIsPrimary->GetXaxis()->SetBinLabel(7,">=4");
		 TotalNBins+=10;
		 
		if(fAnalyseAddedSignals) {
			fH1Pi0DecayModesAddedPi0 = new TH1F("fH1Pi0DecayModesAddedPi0","fH1Pi0DecayModesAddedPi0",10,0.5,10.5);
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetTitle("type of decay / reaction");
			fH1Pi0DecayModesAddedPi0->GetYaxis()->SetTitle("##pi^{0}");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(1,"0");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(2,"1");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
			fH1Pi0DecayModesAddedPi0->GetXaxis()->SetBinLabel(7,">=4");
			TotalNBins+=10;
			
			fH1Pi0DecayModesAddedPi0PHOS = new TH1F("fH1Pi0DecayModesAddedPi0PHOS","fH1Pi0DecayModesAddedPi0PHOS",10,0.5,10.5);
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetTitle("type of decay / reaction");
			fH1Pi0DecayModesAddedPi0PHOS->GetYaxis()->SetTitle("##pi^{0}");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(1,"0");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(2,"1");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
			fH1Pi0DecayModesAddedPi0PHOS->GetXaxis()->SetBinLabel(7,">=4");
			TotalNBins+=10;
		}
		 
		 if(fFillNDaughtersVsPtAndR) {
			 fH2NDaughtersVsPtSimPi0s = new TH2F("fH2NDaughtersVsPtSimPi0s","fH2NDaughtersVsPtSimPi0s",11,-1,10,ptbins, ptlow, ptup);
			 fH2NDaughtersVsPtSimPi0s->GetXaxis()->SetTitle("nDaughters");
			 fH2NDaughtersVsPtSimPi0s->GetYaxis()->SetTitle("pT");
			 TotalNBins+=6;
			 fH2NDaughtersVsdRSimPi0s = new TH2F("fH2NDaughtersVsdRSimPi0s","fH2NDaughtersVsdRSimPi0s",11,-1,10,20000, 0, 10);
			 fH2NDaughtersVsdRSimPi0s->GetXaxis()->SetTitle("nDaughters");
			 fH2NDaughtersVsdRSimPi0s->GetYaxis()->SetTitle("dR");
			 TotalNBins+=20000*11;
		 }
		 
		 fH1MCpionVertDistToEventVert = new TH1F("fH1MCpionVertDistToEventVert", "fH1MCpionVertDistToEventVert", 50000, 0, 5000);
		 fH1MCpionVertDistToEventVert->GetXaxis()->SetTitle("dR");
		 fH1MCpionVertDistToEventVert->GetYaxis()->SetTitle("counts");
		 fH1MCpionVertDistToEventVert->SetMarkerStyle(kFullCircle);
		 TotalNBins+=10000;
		 
		 fH1MCpionVertDistToEventVertLowDR = new TH1F("fH1MCpionVertDistToEventVertLowDR", "fH1MCpionVertDistToEventVertLowDR", 3000, 0, 3);
		 fH1MCpionVertDistToEventVertLowDR->GetXaxis()->SetTitle("dR");
		 fH1MCpionVertDistToEventVertLowDR->GetYaxis()->SetTitle("counts");
		 fH1MCpionVertDistToEventVertLowDR->SetMarkerStyle(kFullCircle);
		 TotalNBins+=10000;
		 
		 fH1MCpionVertDistToEventVertIsPrimary = new TH1F("fH1MCpionVertDistToEventVertIsPrimary", "fH1MCpionVertDistToEventVertIsPrimary", 3000, 0, 3);
		 fH1MCpionVertDistToEventVertIsPrimary->GetXaxis()->SetTitle("dR");
		 fH1MCpionVertDistToEventVertIsPrimary->GetYaxis()->SetTitle("counts");
		 fH1MCpionVertDistToEventVertIsPrimary->SetMarkerStyle(kFullCircle);
		 TotalNBins+=10000;

		 fH1Pi0TruthPt = new TH1F("fH1Pi0TruthPt", "P_{T} distribution for Truth Pi0's", ptbins, ptlow, ptup);
		 fH1Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		 fH1Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		 fH1Pi0TruthPt->SetMarkerStyle(kFullCircle);
		 TotalNBins+=ptbins;

		//if(fFillFeedDownHistos)
		//{
			fH1K0Pi0TruthPt = new TH1F("fH1K0Pi0TruthPt", "P_{T} distribution for Truth Pi0's from K^{0}_{s} decays", ptbins, ptlow, ptup);
			fH1K0Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1K0Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1K0Pi0TruthPt->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;

			fH1PriPi0TruthPt = new TH1F("fH1PriPi0TruthPt", "P_{T} distribution for Truth Primary Pi0's", ptbins, ptlow, ptup);
			fH1PriPi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1PriPi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1PriPi0TruthPt->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;

			fH1SecPi0TruthPt = new TH1F("fH1SecPi0TruthPt", "P_{T} distribution for Truth Secondary Pi0's (without K0-decays)", ptbins, ptlow, ptup);
			fH1SecPi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1SecPi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1SecPi0TruthPt->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;
		
			fH1K0Pi0TruthPtPhi2PiY05 = new TH1F("fH1K0Pi0TruthPtPhi2PiY05", "P_{T} distribution for Truth Pi0's from K^{0}_{s} decays Phi2PiY05", ptbins, ptlow, ptup);
			fH1K0Pi0TruthPtPhi2PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1K0Pi0TruthPtPhi2PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1K0Pi0TruthPtPhi2PiY05->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;

			fH1PriPi0TruthPtPhi2PiY05 = new TH1F("fH1PriPi0TruthPtPhi2PiY05", "P_{T} distribution for Truth Primary Pi0's Phi2PiY05 ", ptbins, ptlow, ptup);
			fH1PriPi0TruthPtPhi2PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1PriPi0TruthPtPhi2PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1PriPi0TruthPtPhi2PiY05->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;

			fH1SecPi0TruthPtPhi2PiY05 = new TH1F("fH1SecPi0TruthPtPhi2PiY05", "P_{T} distribution for Truth Secondary Pi0's (without K0-decays) Phi2PiY05", ptbins, ptlow, ptup);
			fH1SecPi0TruthPtPhi2PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fH1SecPi0TruthPtPhi2PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fH1SecPi0TruthPtPhi2PiY05->SetMarkerStyle(kFullCircle);
			TotalNBins+=ptbins;
		 //}

		
		fH1Pi0TruthPtPhos = new TH1F("fH1Pi0TruthPtPhos", "P_{T} distribution for Truth Pi0's (hit Phos)", ptbins, ptlow, ptup);
		fH1Pi0TruthPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1Pi0TruthPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1Pi0TruthPtPhos->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;
		
		fH1K0Pi0TruthPtPhos = new TH1F("fH1K0Pi0TruthPtPhos", "P_{T} distribution for Truth Pi0's from K^{0}_{s} decays (hit Phos)", ptbins, ptlow, ptup);
		fH1K0Pi0TruthPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1K0Pi0TruthPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1K0Pi0TruthPtPhos->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;
		
		fH1PriPi0TruthPtPhos = new TH1F("fH1PriPi0TruthPtPhos", "P_{T} distribution for Truth Primary Pi0's (hit Phos)", ptbins, ptlow, ptup);
		fH1PriPi0TruthPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1PriPi0TruthPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1PriPi0TruthPtPhos->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;
		
		fH1SecPi0TruthPtPhos = new TH1F("fH1SecPi0TruthPtPhos", "P_{T} distribution for Truth Secondary Pi0's (hit Phos, no K0-decays)", ptbins, ptlow, ptup);
		fH1SecPi0TruthPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1SecPi0TruthPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1SecPi0TruthPtPhos->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;

		fH1Pi0TruthPtPhi2PiY05 = new TH1F("fH1Pi0TruthPtPhi2PiY05", "P_{T} for Truth Pi0's [|y_{#pi^{0}}|<0.5 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
		fH1Pi0TruthPtPhi2PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1Pi0TruthPtPhi2PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1Pi0TruthPtPhi2PiY05->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;

		fH1Pi0TruthPtPhi2PiY05dR1 = new TH1F("fH1Pi0TruthPtPhi2PiY05dR1", "P_{T} for Truth Pi0's [|y_{#pi^{0}}|<0.5 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
		fH1Pi0TruthPtPhi2PiY05dR1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1Pi0TruthPtPhi2PiY05dR1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1Pi0TruthPtPhi2PiY05dR1->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;

		fH1PriPi0TruthPtPhi2PiY05dR1 = new TH1F("fH1PriPi0TruthPtPhi2PiY05dR1", "P_{T} for Truth Pi0's [|y_{#pi^{0}}|<0.5 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
		fH1PriPi0TruthPtPhi2PiY05dR1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1PriPi0TruthPtPhi2PiY05dR1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1PriPi0TruthPtPhi2PiY05dR1->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;
		
		 
		fH1Pi0TruthPtPhi2piY03 = new TH1F("fH1Pi0TruthPtPhi2piY03", "P_{T} for Truth Pi0's [|y_{#pi^{0}}|<0.3 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
		fH1Pi0TruthPtPhi2piY03->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		fH1Pi0TruthPtPhi2piY03->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		fH1Pi0TruthPtPhi2piY03->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;

		fH2Pi0TruthPhiEta = new TH2F("fH2Pi0TruthPhiEta","Pi0Truth Phi vs Eta ", 380,-0.02,6.30, 200,-10,10);
		fH2Pi0TruthPhiEta->GetXaxis()->SetTitle("#phi [rad]");
		fH2Pi0TruthPhiEta->GetYaxis()->SetTitle("#eta ");
		TotalNBins+=380*200;

	
		fH1ElectronConversionR = new TH1F("fH1ElectronConversionR", "conversion point (radius)", 600,0,600);
		fH1ElectronConversionR->GetXaxis()->SetTitle("production radius [cm]");
		fH1ElectronConversionR->SetMarkerStyle(kFullCircle);
		TotalNBins+=600; 

		//// FILLING NEEDS TO BE IMPLEMENTED
		//fH1Pi0TruthPtPhotonsPhos = new TH1F("fH1Pi0TruthPtPhotonsPhos", "P_{T} distribution for Pions with both photons (in Phos)", ptbins, ptlow, ptup);
		//fH1Pi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		//fH1Pi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		//fH1Pi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
		//TotalNBins+=ptbins;

		//// FILLING NEEDS TO BE IMPLEMENTED
		//fH1K0Pi0TruthPtPhotonsPhos = new TH1F("fH1K0Pi0TruthPtPhotonsPhos", "P_{T} distribution for Pions from K0 with both photons (in Phos)", ptbins, ptlow, ptup);
		//fH1K0Pi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		//fH1K0Pi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		//fH1K0Pi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
		//TotalNBins+=ptbins;

		//// FILLING NEEDS TO BE IMPLEMENTED
		//fH1PriPi0TruthPtPhotonsPhos = new TH1F("fH1PriPi0TruthPtPhotonsPhos", "P_{T} distribution for primary Pions with both photons (in Phos)", ptbins, ptlow, ptup);
		//fH1PriPi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		//fH1PriPi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		//fH1PriPi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
		//TotalNBins+=ptbins;

		//// FILLING NEEDS TO BE IMPLEMENTED
		//fH1SecPi0TruthPtPhotonsPhos = new TH1F("fH1SecPi0TruthPtPhotonsPhos", "P_{T} distribution for secondary Pions with both photons (in Phos)", ptbins, ptlow, ptup);
		//fH1SecPi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
		//fH1SecPi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
		//fH1SecPi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
		//TotalNBins+=ptbins;
		
		
		if(fAnalyseAddedSignals) {
			
			if(fAnalyseMCPi0) {
				fH2YVsPhiGenPi0 = new TH2F("fH2YVsPhiGenPi0","y vs phi distribution of pi0 Hijing",300, 0., 6.4,500,-10,10);
				fH2YVsPhiGenPi0->GetYaxis()->SetTitle("y");
				fH2YVsPhiGenPi0->GetXaxis()->SetTitle("phi");
			}
			if(fAnalyseMCEta) {
				fH2YVsPhiGenEta = new TH2F("fH2YVsPhiGenEta","y vs phi distribution of eta Hijing",300, 0., 6.4,500,-10,10);
				fH2YVsPhiGenEta->GetYaxis()->SetTitle("y");
				fH2YVsPhiGenEta->GetXaxis()->SetTitle("phi");
			}
			if(fAnalyseMCPi0) {
				fH2YVsPhiAddedPi0 = new TH2F("fH2YVsPhiAddedPi0","y vs phi distribution for pi0_1",300, 0., 6.4,500,-10,10);
				fH2YVsPhiAddedPi0->GetXaxis()->SetTitle("phi");
				fH2YVsPhiAddedPi0->GetYaxis()->SetTitle("y");
			}
			if(fAnalyseMCEta) {
				fH2YVsPhiAddedEta = new TH2F("fH2YVsPhiAddedEta","y vs phi distribution for eta_2",300, 0., 6.4,500,-10,10);
				fH2YVsPhiAddedEta->GetXaxis()->SetTitle("phi");
				fH2YVsPhiAddedEta->GetYaxis()->SetTitle("y");
			}
			if(fAnalyseMCPi0) {
				fH2YVsPhiAddedPi0PHOS = new TH2F("fH2YVsPhiAddedPi0PHOS","y vs phi distribution for pi0PHS_4",500, 0., 6.4,500,-10,10);
				fH2YVsPhiAddedPi0PHOS->GetXaxis()->SetTitle("phi");
				fH2YVsPhiAddedPi0PHOS->GetYaxis()->SetTitle("y");
			}
			if(fAnalyseMCEta) {
				fH2YVsPhiAddedEtaPHOS = new TH2F("fH2YVsPhiAddedEtaPHOS","y vs phi distribution for etaPHS_6",500, 0., 6.4,500,-10,10);
				fH2YVsPhiAddedEtaPHOS->GetXaxis()->SetTitle("phi");
				fH2YVsPhiAddedEtaPHOS->GetYaxis()->SetTitle("y");
			}
			
			if(fAnalyseMCPi0) { 
				fH1YAddedPi0 = new TH1F("fH1YAddedPi0", "rapidity distribution of pi0_1", 2600,-1.3,1.3);
				fH1YAddedPi0->GetXaxis()->SetTitle("y");
				fH1YAddedPi0->GetYaxis()->SetTitle("counts");
				fH1YAddedPi0->SetMarkerStyle(kFullCircle);
				
				fH1EtaAddedPi0 = new TH1F("fH1EtaAddedPi0", "pseudorapidity distribution of pi0_1", 2600,-1.3,1.3);
				fH1EtaAddedPi0->GetXaxis()->SetTitle("#eta");
				fH1EtaAddedPi0->GetYaxis()->SetTitle("counts");
				fH1EtaAddedPi0->SetMarkerStyle(kFullCircle);
				
				fH1YAddedPi0PHOS = new TH1F("fH1YAddedPi0PHOS", "rapidity distribution of pi0PHOS_4", 2600,-1.3,1.3);
				fH1YAddedPi0PHOS->GetXaxis()->SetTitle("y");
				fH1YAddedPi0PHOS->GetYaxis()->SetTitle("counts");
				fH1YAddedPi0PHOS->SetMarkerStyle(kFullCircle);
				
				fH1EtaAddedPi0PHOS = new TH1F("fH1EtaAddedPi0PHOS", "pseudorapidity distribution of pi0PHOS_4", 2600,-1.3,1.3);
				fH1EtaAddedPi0PHOS->GetXaxis()->SetTitle("#eta");
				fH1EtaAddedPi0PHOS->GetYaxis()->SetTitle("counts");
				fH1EtaAddedPi0PHOS->SetMarkerStyle(kFullCircle);
				
				fH1ClusterEAddedPi0 = new TH1F("fH1ClusterEAddedPi0", "Cluster Energy in Phos Added Pi0 ", Ebins, Elow, Eup);
				fH1ClusterEAddedPi0->GetXaxis()->SetTitle("E [GeV]");
				fH1ClusterEAddedPi0->GetYaxis()->SetTitle("counts");
				fH1ClusterEAddedPi0->SetMarkerStyle(kFullCircle);
				fH1ClusterEAddedPi0->Sumw2();
				TotalNBins+=Ebins;
	
				fH1ClusterEAddedPi0PHOS = new TH1F("fH1ClusterEAddedPi0PHOS", "Cluster Energy in Phos Added Pi0PHOS", Ebins, Elow, Eup);
				fH1ClusterEAddedPi0PHOS->GetXaxis()->SetTitle("E [GeV]");
				fH1ClusterEAddedPi0PHOS->GetYaxis()->SetTitle("counts");
				fH1ClusterEAddedPi0PHOS->SetMarkerStyle(kFullCircle);
				fH1ClusterEAddedPi0PHOS->Sumw2();
				TotalNBins+=Ebins;
				
				if(fFillDecayPhotonInfoAddedSig) {
					
					
					fH1DecGammAddPi0Eta = new TH1F("fH1DecGammAddPi0Eta", "pseudorapidity distribution of gammas of pi0_1", 1000,-0.5,0.5);
					fH1DecGammAddPi0Eta->GetXaxis()->SetTitle("#eta");
					fH1DecGammAddPi0Eta->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0Y = new TH1F("fH1DecGammAddPi0Y", "rapidity distribution of gammas of pi0_1", 1000,-0.5,0.5);
					fH1DecGammAddPi0Y->GetXaxis()->SetTitle("y");
					fH1DecGammAddPi0Y->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0Phi = new TH1F("fH1DecGammAddPi0Phi", "rapidity distribution of gammas of pi0_1", 720,0,360);
					fH1DecGammAddPi0Phi->GetXaxis()->SetTitle("#phi");
					fH1DecGammAddPi0Phi->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0E = new TH1F("fH1DecGammAddPi0E", "Energy of gammas of pi0_1 ", Ebins, Elow, Eup);
					fH1DecGammAddPi0E->GetXaxis()->SetTitle("E [GeV]");
					fH1DecGammAddPi0E->GetYaxis()->SetTitle("counts");
					fH1DecGammAddPi0E->SetMarkerStyle(kFullCircle);
					fH1DecGammAddPi0E->Sumw2();
					
					fH1DecGammInPHOSAddPi0E = new TH1F("fH1DecGammInPHOSAddPi0E", "Energy of gammas pointing to PHOS of pi0_1 ", Ebins, Elow, Eup);
					fH1DecGammInPHOSAddPi0E->GetXaxis()->SetTitle("E [GeV]");
					fH1DecGammInPHOSAddPi0E->GetYaxis()->SetTitle("counts");
					fH1DecGammInPHOSAddPi0E->SetMarkerStyle(kFullCircle);
					fH1DecGammInPHOSAddPi0E->Sumw2();
					
					fH1DecGammAddPi0Asymm = new TH1F("fH1DecGammAddPi0Asymm", "Asymmetry of gammas of pi0_1 ", 300, 0, 1);
					fH1DecGammAddPi0Asymm->GetXaxis()->SetTitle("abs(E1-E2)/(E1+E1)");
					fH1DecGammAddPi0Asymm->GetYaxis()->SetTitle("counts");
					fH1DecGammAddPi0Asymm->SetMarkerStyle(kFullCircle);
					
					fH1DecGammAddPi0OpAngle = new TH1F("fH1DecGammAddPi0OpAngle", "Opening angle of gammas of pi0_1 ", 400,0.,1.1*TMath::Pi());
					fH1DecGammAddPi0OpAngle->GetXaxis()->SetTitle("opening angle");
					fH1DecGammAddPi0OpAngle->GetYaxis()->SetTitle("counts");
					
					fH1DecGammAddPi0ConvR = new TH1F("fH1DecGammAddPi0ConvR", "conversion point (radius) of gammas of pi0_1", 1400,0,700);
					fH1DecGammAddPi0ConvR->GetXaxis()->SetTitle("production radius [cm]");
					fH1DecGammAddPi0ConvR->SetMarkerStyle(kFullCircle);					
					
					fH1DecGammAddPi0ConvRate = new TH1F("fH1DecGammAddPi0ConvRate", "gamma from pi0_1 converted?", 2,-0.5,1.5);
					fH1DecGammAddPi0ConvRate->GetXaxis()->SetTitle("production radius [cm]");
					fH1DecGammAddPi0ConvRate->GetXaxis()->SetBinLabel(1,"yes");
					fH1DecGammAddPi0ConvRate->GetXaxis()->SetBinLabel(2,"no");

					fH1TruthPtAddedPi0GammasInPHOS = new TH1F("fH1TruthPtAddedPi0GammasInPHOS", "P_{T} distribution for pi0_1 gammas in PHOS", ptbins, ptlow, ptup);
					fH1TruthPtAddedPi0GammasInPHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
					fH1TruthPtAddedPi0GammasInPHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
					fH1TruthPtAddedPi0GammasInPHOS->SetMarkerStyle(kFullCircle);
					fH1TruthPtAddedPi0GammasInPHOS->Sumw2();
					
					
			
					fH1DecGammAddPi0PHOSEta = new TH1F("fH1DecGammAddPi0PHOSEta", "pseudorapidity distribution of gammas of pi0PHOS_4", 1000,-0.5,0.5);
					fH1DecGammAddPi0PHOSEta->GetXaxis()->SetTitle("#eta");
					fH1DecGammAddPi0PHOSEta->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0PHOSY = new TH1F("fH1DecGammAddPi0PHOSY", "rapidity distribution of gammas of pi0PHOS_4", 1000,-0.5,0.5);
					fH1DecGammAddPi0PHOSY->GetXaxis()->SetTitle("y");
					fH1DecGammAddPi0PHOSY->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0PHOSPhi = new TH1F("fH1DecGammAddPi0PHOSPhi", "rapidity distribution of gammas of pi0PHOS_4", 720,0,360);
					fH1DecGammAddPi0PHOSPhi->GetXaxis()->SetTitle("#phi");
					fH1DecGammAddPi0PHOSPhi->GetYaxis()->SetTitle("Counts");
					
					fH1DecGammAddPi0PHOSE = new TH1F("fH1DecGammAddPi0PHOSE", "Energy of gammas of pi0PHOS_4 ", Ebins, Elow, Eup);
					fH1DecGammAddPi0PHOSE->GetXaxis()->SetTitle("E [GeV]");
					fH1DecGammAddPi0PHOSE->GetYaxis()->SetTitle("counts");
					fH1DecGammAddPi0PHOSE->SetMarkerStyle(kFullCircle);
					fH1DecGammAddPi0PHOSE->Sumw2();
					
					fH1DecGammInPHOSAddPi0PHOSE = new TH1F("fH1DecGammInPHOSAddPi0PHOSE", "Energy of gammas pointing to PHOS of pi0PHOS_4 ", Ebins, Elow, Eup);
					fH1DecGammInPHOSAddPi0PHOSE->GetXaxis()->SetTitle("E [GeV]");
					fH1DecGammInPHOSAddPi0PHOSE->GetYaxis()->SetTitle("counts");
					fH1DecGammInPHOSAddPi0PHOSE->SetMarkerStyle(kFullCircle);
					fH1DecGammInPHOSAddPi0PHOSE->Sumw2();
					
					fH1DecGammAddPi0PHOSAsymm = new TH1F("fH1DecGammAddPi0PHOSAsymm", "Asymmetry of gammas of pi0PHOS_4 ", 300, 0, 1);
					fH1DecGammAddPi0PHOSAsymm->GetXaxis()->SetTitle("abs(E1-E2)/(E1+E1)");
					fH1DecGammAddPi0PHOSAsymm->GetYaxis()->SetTitle("counts");
					fH1DecGammAddPi0PHOSAsymm->SetMarkerStyle(kFullCircle);
					
					fH1DecGammAddPi0PHOSOpAngle = new TH1F("fH1DecGammAddPi0PHOSOpAngle", "Opening angle of gammas of pi0PHOS_4 ", 400,0.,1.1*TMath::Pi());
					fH1DecGammAddPi0PHOSOpAngle->GetXaxis()->SetTitle("opening angle");
					fH1DecGammAddPi0PHOSOpAngle->GetYaxis()->SetTitle("counts");
					
					fH1DecGammAddPi0PHOSConvR = new TH1F("fH1DecGammAddPi0PHOSConvR", "conversion point (radius) of gammas of pi0PHOS_4", 1400,0,700);
					fH1DecGammAddPi0PHOSConvR->GetXaxis()->SetTitle("production radius [cm]");
					fH1DecGammAddPi0PHOSConvR->SetMarkerStyle(kFullCircle);					
					
					fH1DecGammAddPi0PHOSConvRate = new TH1F("fH1DecGammAddPi0PHOSConvRate", "gamma from pi0PHOS_4 converted?", 2,-0.5,1.5);
					fH1DecGammAddPi0PHOSConvRate->GetXaxis()->SetTitle("");
					fH1DecGammAddPi0PHOSConvRate->GetXaxis()->SetBinLabel(1,"yes");
					fH1DecGammAddPi0PHOSConvRate->GetXaxis()->SetBinLabel(2,"no");

					fH1TruthPtAddedPi0PHOSGammasInPHOS = new TH1F("fH1TruthPtAddedPi0PHOSGammasInPHOS", "P_{T} distribution for pi0PHOS_4 gammas in PHOS", ptbins, ptlow, ptup);
					fH1TruthPtAddedPi0PHOSGammasInPHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
					fH1TruthPtAddedPi0PHOSGammasInPHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
					fH1TruthPtAddedPi0PHOSGammasInPHOS->SetMarkerStyle(kFullCircle);
					fH1TruthPtAddedPi0PHOSGammasInPHOS->Sumw2();
			
			
				}
				
			}
			

			if(fAnalyseMCPi0) {
				fH1TruthPtGenPi0 = new TH1F("fH1TruthPtGenPi0", "P_{T} distribution of pi0 hijing", ptbins, ptlow, ptup);
				fH1TruthPtGenPi0->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtGenPi0->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtGenPi0->SetMarkerStyle(kFullCircle);
			}
			if(fAnalyseMCEta) {
				fH1TruthPtGenEta = new TH1F("fH1TruthPtGenEta", "P_{T} distribution of eta hijing", ptbins, ptlow, ptup);
				fH1TruthPtGenEta->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtGenEta->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtGenEta->SetMarkerStyle(kFullCircle);
			}
			if(fAnalyseMCPi0) {
				fH1TruthPtAddedPi0 = new TH1F("fH1TruthPtAddedPi0", "P_{T} distribution for pi0_1", ptbins, ptlow, ptup);
				fH1TruthPtAddedPi0->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedPi0->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedPi0->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedPi0->Sumw2();
				
				fH1TruthPtAddedPi02PiY05 = new TH1F("fH1TruthPtAddedPi02PiY05", "P_{T} distribution for pi0_1 phi2pieta05", ptbins, ptlow, ptup);
				fH1TruthPtAddedPi02PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedPi02PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedPi02PiY05->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedPi02PiY05->Sumw2();
				
				
				fH1TruthPtAddedPi0MesonInPHOS = new TH1F("fH1TruthPtAddedPi0MesonInPHOS", "P_{T} distribution for pi0_1 went into PHOS acc", ptbins, ptlow, ptup);
				fH1TruthPtAddedPi0MesonInPHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedPi0MesonInPHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedPi0MesonInPHOS->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedPi0MesonInPHOS->Sumw2();
				
				
				fH1TruthPtAddedPi0PHOSMesonInPHOS = new TH1F("fH1TruthPtAddedPi0PHOSMesonInPHOS", "P_{T} distribution for pi0PHOS_4 went into PHOS acc", ptbins, ptlow, ptup);
				fH1TruthPtAddedPi0PHOSMesonInPHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedPi0PHOSMesonInPHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedPi0PHOSMesonInPHOS->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedPi0PHOSMesonInPHOS->Sumw2();
			}
			if(fAnalyseMCEta) {
				fH1TruthPtAddedEta = new TH1F("fH1TruthPtAddedEta", "P_{T} distribution for pi0_1", ptbins, ptlow, ptup);
				fH1TruthPtAddedEta->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedEta->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedEta->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedEta->Sumw2();
			}
			if(fAnalyseMCPi0) {
				fH1TruthPtAddedPi0PHOS = new TH1F("fH1TruthPtAddedPi0PHOS", "P_{T} distribution for pi0PHS_4", ptbins, ptlow, ptup);
				fH1TruthPtAddedPi0PHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedPi0PHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedPi0PHOS->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedPi0PHOS->Sumw2();
			}
			if(fAnalyseMCEta) {
				fH1TruthPtAddedEtaPHOS = new TH1F("fH1TruthPtAddedEtaPHOS", "P_{T} distribution for etaPHS_6", ptbins, ptlow, ptup);
				fH1TruthPtAddedEtaPHOS->GetXaxis()->SetTitle("P_{T} (GeV/c)");
				fH1TruthPtAddedEtaPHOS->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
				fH1TruthPtAddedEtaPHOS->SetMarkerStyle(kFullCircle);
				fH1TruthPtAddedEtaPHOS->Sumw2();
			}
			if(fAnalyseMCPi0) {
				fH2MPtAddedPi0	= new TH2F("fH2MPtAddedPi0", "m_{inv} vs P_{T} for reconstructed added Pi0 (shot in wide range)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedPi0->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedPi0->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedPi0->Sumw2();
			}
			
			if(fAnalyseMCPi0) {
				fH2MPtAddedPi0_unweighed	= new TH2F("fH2MPtAddedPi0_unweighed", "m_{inv} vs P_{T} for reconstructed added Pi0 (shot in wide range)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedPi0_unweighed->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedPi0_unweighed->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedPi0_unweighed->Sumw2();
				
				fH2MPtAddedPi0PHOS	= new TH2F("fH2MPtAddedPi0PHOS", "m_{inv} vs P_{T} for reconstructed added Pi0 (shot at PHOS directly)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedPi0PHOS->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedPi0PHOS->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedPi0PHOS->Sumw2();
				
				fH2MPtAddedPi0PHOS_unweighed	= new TH2F("fH2MPtAddedPi0PHOS_unweighed", "m_{inv} vs P_{T} for reconstructed added Pi0 (shot at PHOS directly)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedPi0PHOS_unweighed->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedPi0PHOS_unweighed->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedPi0PHOS_unweighed->Sumw2();
				
				if(fFillMPtForSingleOrMultContrClus) {
					fH2MPtAddedPi0PHOSSingleContr	= new TH2F("fH2MPtAddedPi0PHOSSingleContr", "m_{inv} vs P_{T} for rec. added Pi0 both cluster single contributer (shot at PHOS)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
					fH2MPtAddedPi0PHOSSingleContr->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
					fH2MPtAddedPi0PHOSSingleContr->GetYaxis()->SetTitle("p_{T} [GeV/c]");
					fH2MPtAddedPi0PHOSSingleContr->Sumw2();
					
					fH2MPtAddedPi0PHOSMultContr	= new TH2F("fH2MPtAddedPi0PHOSMultContr", "m_{inv} vs P_{T} for rec. added Pi0, both cluster mult contributer (shot at PHOS)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
					fH2MPtAddedPi0PHOSMultContr->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
					fH2MPtAddedPi0PHOSMultContr->GetYaxis()->SetTitle("p_{T} [GeV/c]");
					fH2MPtAddedPi0PHOSMultContr->Sumw2();
				}
				
				fH2MEnergyDiffAddedPi0PHOS	= new TH2F("fH2MEnergyDiffAddedPi0PHOS", "m_{inv} vs P_{T} for reconstructed added Pi0 (shot at PHOS directly)",Mbins, Mlow, Mup, 600,-30,30);
				fH2MEnergyDiffAddedPi0PHOS->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MEnergyDiffAddedPi0PHOS->GetYaxis()->SetTitle("+/- max(abs(E^{clu}_{1}-E^{#gamma}_{1}),abs(E^{clu}_{2}-E^{#gamma}_{2}))");
				fH2MEnergyDiffAddedPi0PHOS->Sumw2();
			
				fH2PtRecVsPtTruthAddedPi0PHOS = new TH2F("fH2PtRecVsPtTruthAddedPi0PHOS","fH2PtRecVsPtTruthAddedPi0PHOS",ptbins, ptlow, ptup, ptbins, ptlow, ptup);
				fH2PtRecVsPtTruthAddedPi0PHOS->GetXaxis()->SetTitle("true p_{T} [GeV]");
				fH2PtRecVsPtTruthAddedPi0PHOS->GetYaxis()->SetTitle("rec p_{T} [GeV]");
				fH2PtRecVsPtTruthAddedPi0PHOS->Sumw2();
			}
			
			if(fAnalyseMCEta) {
				fH2MPtAddedEta	= new TH2F("fH2MPtAddedEta", "m_{inv} vs P_{T} for reconstructed added Eta (shot in wide range)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedEta->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedEta->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedEta->Sumw2();
				
				fH2MPtAddedEtaPHOS	= new TH2F("fH2MPtAddedEtaPHOS", "m_{inv} vs P_{T} for reconstructed added Eta (shot at PHOS directly)",Mbins, Mlow, Mup, ptbins, ptlow, ptup);
				fH2MPtAddedEtaPHOS->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
				fH2MPtAddedEtaPHOS->GetYaxis()->SetTitle("p_{T} [GeV/c]");
				fH2MPtAddedEtaPHOS->Sumw2();
			}
			
		}
		
		fH2ClusterEVSPhotonE = new TH2F("fH2ClusterEVSPhotonE","fH2ClusterEVSPhotonE",Ebins, Elow, Eup,Ebins, Elow, Eup);
		fH2ClusterEVSPhotonE->GetXaxis()->SetTitle("E_{Photon} [GeV]");
		fH2ClusterEVSPhotonE->GetYaxis()->SetTitle("E_{Cluster} [GeV]");
		
		
	} //fFillMCHistos


	const Int_t etabins = 1000;
	const Float_t etalow = -1.6;
	const Float_t etaup = 1.6;

	const Int_t phibins = 1200;
	const Float_t philow = -0.0;
	const Float_t phiup = -5.0;
	
	fH2HitmapEtaVsPhi = new TH2D("fH2HitmapEtaVsPhi","#eta vs phi distribution for reconstructed",etabins, etalow, etaup,phibins,phiup,philow);
	fH2HitmapEtaVsPhi->GetXaxis()->SetTitle("eta");
	fH2HitmapEtaVsPhi->GetYaxis()->SetTitle("phi");
	TotalNBins+=400*600;

    Int_t chi2bins = 100;
    Float_t chi2low = -2, chi2up = 2;
    fH1Chi2 = new TH1F("fH1Chi2","#chi^{2} distribution for reconstructed",chi2bins, chi2low, chi2up);
    fH1Chi2->GetXaxis()->SetTitle("#chi^{2}");
    fH1Chi2->GetYaxis()->SetTitle("counts");
    TotalNBins+=chi2bins;

    fH1NTrkMatch = new TH1F("fH1NTrkMatch","number of matched tracks",14, -1.5, 5.5);
    fH1NTrkMatch->GetXaxis()->SetTitle("nTracksMatched");
    fH1NTrkMatch->GetYaxis()->SetTitle("counts");
    TotalNBins+=14;

    fH1ClusterDisp = new TH1F("fH1ClusterDisp","Dispersion of CaloCluster",1000, -1, 3);
    fH1ClusterDisp->GetXaxis()->SetTitle("cluster->GetClusterDisp()");
    fH1ClusterDisp->GetYaxis()->SetTitle("counts");
    TotalNBins+=1000;

    fH2Ellipse = new TH2F("fH2Ellipse","Ellipse axis M20 vs M02",500, -0.01, 1, 500, -0.01, 1);
    fH2Ellipse->GetXaxis()->SetTitle("cluster->GetM20()");
    fH2Ellipse->GetYaxis()->SetTitle("cluster->GetM02()");
    fH2Ellipse->GetZaxis()->SetTitle("counts");
    TotalNBins+=500*500;

   if(fFillHMassPtModules) {
		fH3MPtModules = new TH3F("fH3MPtModules","mass vs p_{T} vs Module Combinations",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 6,0.5,6.5);
		fH3MPtModules->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtModules->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtModules->GetZaxis()->SetTitle("Module Combinations: 11,12,13,22,23,33");
		TotalNBins+=Mbins*ptbins*6.0;

		fH3MPtModulesMix = new TH3F("fH3MPtModulesMix","mass vs p_{T} vs Module Combinations (mixed events)",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 6,0.5,6.5);
		fH3MPtModulesMix->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtModulesMix->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtModulesMix->GetZaxis()->SetTitle("Module Combinations: 11,12,13,22,23,33");
		TotalNBins+=Mbins*ptbins*6.0;
		
   } 
   
   else if(fFillHMassPtTiming) {
		fH3MPtTiming = new TH3F("fH3MPtTiming", "mass vs p_{T} vs Timing-Cut Combinations", Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5); 
		fH3MPtTiming->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtTiming->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtTiming->GetZaxis()->SetTitle("Timing-Cut: Both Good, One Good, Both Bad");
		TotalNBins+=Mbins*ptbins*3.0;
		
		fH3MPtTimingMix = new TH3F("fH3MPtTimingMix", "mass vs p_{T} vs Timing-Cut Combinations", Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5); 
		fH3MPtTimingMix->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtTimingMix->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtTimingMix->GetZaxis()->SetTitle("Timing-Cut: Both Good, One Good, Both Bad");
		TotalNBins+=Mbins*ptbins*3.0;	
   }
   // *** default ***
   else {  
		fH3MPtAsymm = new TH3F("fH3MPtAsymm","mass vs p_{T} vs Asymm cut",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
		fH3MPtAsymm->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtAsymm->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtAsymm->GetZaxis()->SetTitle("Asymmetry Cut  (edges: 0.0, 0.1, 0.7, 1.0)");
		TotalNBins+=Mbins*ptbins*3.0;

		fH3MPtAsymmMix = new TH3F("fH3MPtAsymmMix","mass vs p_{T} vs Asymm cut (mixed events)",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
		fH3MPtAsymmMix->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
		fH3MPtAsymmMix->GetYaxis()->SetTitle("p_{T} [GeV/c]");
		fH3MPtAsymmMix->GetZaxis()->SetTitle("Asymmetry Cut  (edges: 0.0, 0.1, 0.7, 1.0)");
		TotalNBins+=Mbins*ptbins*3.0;

		if(fFillNewAsymmClasses) {
			Int_t asymmbins = 3;					// It's not possible to make more than 3 
			Float_t asymmClasses[4]; 			
			asymmClasses[0] = 0; 
			asymmClasses[1] = fAsymmClass1; 
			asymmClasses[2] = fAsymmClass2;
			asymmClasses[3] = fAsymmClass3;

			fH3MPtAsymm->GetZaxis()->Set(asymmbins, asymmClasses);
			fH3MPtAsymm->GetZaxis()->SetTitle("Asymmetry Cut");
			fH3MPtAsymmMix->GetZaxis()->Set(asymmbins, asymmClasses);
			fH3MPtAsymmMix->GetZaxis()->SetTitle("Asymmetry Cut");
		}
	}   
    
    fH2DphiDeta = new TH2F("fH2DphiDeta","#Delta#phi vs #Delta#eta", 349,-1.5,5, 400,-2.0,2.0);
    fH2DphiDeta->GetXaxis()->SetTitle("#Delta#phi");
    fH2DphiDeta->GetYaxis()->SetTitle("#Delta#eta");
    TotalNBins+=349*400;

    fH2DphiDetaMix = new TH2F("fH2DphiDetaMix","#Delta#phi vs #Delta#eta (mixed events)", 349,-1.5,5, 400,-2.0,2.0);
    fH2DphiDetaMix->GetXaxis()->SetTitle("#Delta#phi");
    fH2DphiDetaMix->GetYaxis()->SetTitle("#Delta#eta");
    TotalNBins+=349*400;

    fH2CellsM02 = new TH2F("fH2CellsM02", "nCells vs M02", 204,-1.5,100.5, 500,-1,1.5);
    fH2CellsM02->GetXaxis()->SetTitle("nCells");
    fH2CellsM02->GetYaxis()->SetTitle("M02");
    fH2CellsM02->GetZaxis()->SetTitle("counts");
    TotalNBins+=204*500;
    
    fH1ClusterM02 = new TH1F("fH1ClusterM02", "fH1ClusterM02 vs M02", 2000,-1,4.5);
    fH1ClusterM02->GetXaxis()->SetTitle("M02");
    fH1ClusterM02->GetYaxis()->SetTitle("counts");
    TotalNBins+=500;

    fH1NPrimVertContribut = new TH1F("fH1NPrimVertContribut","Number of Contributors to primary Vertex",130,0,130);
    fH1NPrimVertContribut->GetXaxis()->SetTitle("nCells");
    fH1NPrimVertContribut->GetYaxis()->SetTitle("M");
    TotalNBins+=130;
    
    fH1DistPileUpPrimVert = new TH1F("fH1DistPileUpPrimVert", "prim_vtz - pileUp_vtz (ESD only)", 12000,-30,30);
    fH1DistPileUpPrimVert->GetXaxis()->SetTitle("Z_{0} - Z_{pileup}");
    fH1DistPileUpPrimVert->GetYaxis()->SetTitle("N events");
    TotalNBins+=12000;

    fH1nSPDPileUpVtxs = new TH1F("fH1nSPDPileUpVtxs","Number of PileUp Vertices from SPD (ESD only)",10,0.5,10.5);
    fH1nSPDPileUpVtxs->GetXaxis()->SetTitle("nPileUpVtxices");
    fH1nSPDPileUpVtxs->GetYaxis()->SetTitle("Counts");
    TotalNBins+=10;

    fH1NClustersVsCuts = new TH1F("fH1NClustersVsCuts","Number of Clusters after Cuts",10,0.5,10.5);
    fH1NClustersVsCuts->GetYaxis()->SetTitle("#it{N}_{cl}");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(1,"All PHOS Clus");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(2,"in mod 1-3");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(3,Form("nCells>%i",fClusterMinCells-1));
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(4,Form("E>%.2f",fClusterMinE));
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(5,Form("M02>%.2f",fClusterMinM02));
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(6,"ClusPos good");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(7,"dist_to_BadCell");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(8,"good Timing");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(9,"--");
    fH1NClustersVsCuts->GetXaxis()->SetBinLabel(10,"--");
    TotalNBins+=10;
        
    fTProfMeanClusterEnergyVsCuts  = new TProfile("fTProfMeanClusterEnergyVsCuts","tprofile meanClusterEnergyAfterCuts",10,0.5,10,0,50);
    fTProfMeanClusterEnergyVsCuts->GetYaxis()->SetTitle("#LT#it{E}_{cl}#GT");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(1,"All PHOS Clus");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(2,"in mod 1-3");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(3,Form("nCells>%i",fClusterMinCells-1));
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(4,Form("E>%.2f",fClusterMinE));
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(5,Form("M02>%.2f",fClusterMinM02));
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(6,"ClusPos good");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(7,"dist_to_BadCell");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(8,"good Timing");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(9,"--");
    fTProfMeanClusterEnergyVsCuts->GetXaxis()->SetBinLabel(10,"--");
    TotalNBins+=10;
        
    fH2EAfterCutsVsModNum = new TH2F("fH2EAfterCutsVsModNum","clusterE_vs_ModNumber",Ebins, Elow, Eup, 3,1,4);
    fH2EAfterCutsVsModNum->GetXaxis()->SetTitle("clusterE");
    fH2EAfterCutsVsModNum->GetYaxis()->SetTitle("ModNumber");  
    fH2EAfterCutsVsModNum->GetZaxis()->SetTitle("counts");
    TotalNBins+=(Ebins)*3;

    Int_t nEcl      = 4000;
    Double_t EclMin = 0;
    Double_t EclMax = 20; 
    
    if(fFillCellIdVsE) {
		fH2CellIdVsE = new TH2F("cellID_vs_E", "cellID_vs_E", 11000,0.5,11000.5, nEcl,EclMin,EclMax);
		fH2CellIdVsE->GetXaxis()->SetTitle("CellID");
		fH2CellIdVsE->GetYaxis()->SetTitle("Cell E");
		fH2CellIdVsE->GetZaxis()->SetTitle("counts");
		TotalNBins+=11000*nEcl;
		
		fH2LocalMaxCellsIdVsE = new TH2F("fH2LocalMaxCellsIdVsE", "fH2LocalMaxCellsIdVsE", 11000,0.5,11000.5, nEcl,EclMin,EclMax);
		fH2LocalMaxCellsIdVsE->GetXaxis()->SetTitle("CellID");
		fH2LocalMaxCellsIdVsE->GetYaxis()->SetTitle("Cell Energy");
		fH2LocalMaxCellsIdVsE->GetZaxis()->SetTitle("counts");
		
		fH2ClusterPosCellsIdVsE = new TH2F("fH2ClusterPosCellsIdVsE", "fH2ClusterPosCellsIdVsE", 11000,0.5,11000.5, nEcl,EclMin,EclMax);
		fH2ClusterPosCellsIdVsE->GetXaxis()->SetTitle("CellID");
		fH2ClusterPosCellsIdVsE->GetYaxis()->SetTitle("Cell Energy");
		fH2ClusterPosCellsIdVsE->GetZaxis()->SetTitle("counts");
	}

	if(fFillTimingHistos) {
		fH1ClusterTOFWeightedWithE = new TH1F("fH1ClusterTOFWeightedWithE", "fH1ClusterTOFWeightedWithE", 5000, -6e-6 , 6e-6);
		fH1ClusterTOFWeightedWithE->GetXaxis()->SetTitle("Cluster Time");
		fH1ClusterTOFWeightedWithE->GetYaxis()->SetTitle("counts (weighted with Cluster E)");
		TotalNBins+=12000;

		fH2ClusterTOFVsE = new TH2F("fH2ClusterTOFVsE", "Cluster Time against Cluster Energy", 5000, -6e-6, 6e-6, nEcl, EclMin, EclMax);
		fH2ClusterTOFVsE->GetXaxis()->SetTitle("Cluster-Time");
		fH2ClusterTOFVsE->GetYaxis()->SetTitle("Cluster-Energy");   
	}

    fH1CellMCLabel = new TH1F("fH1CellMCLabel","fH1CellMCLabel",10, -5.5, 4.5);
    fH1CellMCLabel->GetXaxis()->SetTitle("CellMCLabel");
    fH1CellMCLabel->GetYaxis()->SetTitle("Counts");  
    TotalNBins+=10;
    
    fH1AppliedClusterCuts = new TH1F("fH1AppliedClusterCuts","appliedClusterCuts",11,0.5,11.5);
    fH1AppliedClusterCuts->GetXaxis()->SetTitle("");
    fH1AppliedClusterCuts->GetYaxis()->SetTitle("Counts"); 
	
	Char_t saythis[50]; Int_t iter = 1; 
    sprintf(saythis,"nMinCells=%d",fClusterMinCells); 
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"minE=%.3f GeV",fClusterMinE);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fDoDistToBadCellCutOnCellLevel=%d",fDoDistToBadCellCutOnCellLevel);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fDistToBadCellOnCellLevel=%d Cells",fDistToBadCellOnCellLevel);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fDoDistToBadCellCut=%d",fDoDistToBadCellCut);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fDistToBadCell=%.2f cm",fDistToBadCell);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"do timing cut? %d",fDoTimingCut);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"timing_min=%.4fe-6s",fTimingCutMin / 0.000001);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"timing_max=%.4fe-6s",fTimingCutMax/ 0.000001);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fDoZvertexCut? %d",fDoZvertexCut);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    sprintf(saythis,"fZvertexCut-value= %.1f",fZvertexCut);
    fH1AppliedClusterCuts->GetXaxis()->SetBinLabel(iter++,saythis);
    TotalNBins+=11;
    
    fH1DistanceToBadChannel = new TH1F("fH1DistanceToBadChannel","fH1DistanceToBadChannel",1000,0,50);
    fH1DistanceToBadChannel->GetXaxis()->SetTitle("Distance (cm)");
    fH1DistanceToBadChannel->GetYaxis()->SetTitle("Counts"); 
    TotalNBins+=1000;

	// hitmaps 
	if(fFillClusterHitmaps)
	{
		fH2ClusterPositionsMod1 = new TH2F("fH2ClusterPositionsMod1", "fH2ClusterPositionsMod1", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPositionsMod1->SetTitle("cluster positions module 1");
		fH2ClusterPositionsMod1->GetXaxis()->SetTitle("cell X");
		fH2ClusterPositionsMod1->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 

		fH2ClusterPositionsMod2 = new TH2F("fH2ClusterPositionsMod2", "fH2ClusterPositionsMod2", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPositionsMod2->SetTitle("cluster positions module 2");
		fH2ClusterPositionsMod2->GetXaxis()->SetTitle("cell X");
		fH2ClusterPositionsMod2->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 

		fH2ClusterPositionsMod3 = new TH2F("fH2ClusterPositionsMod3", "fH2ClusterPositionsMod3", 64, 0., 64., 56, 0., 56.);
		fH2ClusterPositionsMod3->SetTitle("cluster positions module 3");
		fH2ClusterPositionsMod3->GetXaxis()->SetTitle("cell X");
		fH2ClusterPositionsMod3->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 
	}
	
	if(fFillClusterHitmapsAddedSignals)
	{
		fH2ClusterPosMod1Gen1 = new TH2F("fH2ClusterPosMod1Gen1", "fH2ClusterPosMod1Gen1", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPosMod1Gen1->SetTitle("cluster positions module 1 Generator 1");
		fH2ClusterPosMod1Gen1->GetXaxis()->SetTitle("cell X");
		fH2ClusterPosMod1Gen1->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 
		
		fH2ClusterPosMod3Gen1 = new TH2F("fH2ClusterPosMod3Gen1", "fH2ClusterPosMod3Gen1", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPosMod3Gen1->SetTitle("cluster positions module 1 Generator 1");
		fH2ClusterPosMod3Gen1->GetXaxis()->SetTitle("cell X");
		fH2ClusterPosMod3Gen1->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 
		
		fH2ClusterPosMod1Gen2 = new TH2F("fH2ClusterPosMod1Gen2", "fH2ClusterPosMod1Gen2", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPosMod1Gen2->SetTitle("cluster positions module 1 Generator 1");
		fH2ClusterPosMod1Gen2->GetXaxis()->SetTitle("cell X");
		fH2ClusterPosMod1Gen2->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 
		
		fH2ClusterPosMod3Gen2 = new TH2F("fH2ClusterPosMod3Gen2", "fH2ClusterPosMod3Gen2", 64, 0.5, 64.5, 56, 0.5, 56.5);
		fH2ClusterPosMod3Gen2->SetTitle("cluster positions module 1 Generator 1");
		fH2ClusterPosMod3Gen2->GetXaxis()->SetTitle("cell X");
		fH2ClusterPosMod3Gen2->GetYaxis()->SetTitle("cell Z");
		TotalNBins += 64*56; 
	}
   
	if(fFillNTupelClusterE)
	{
		fNTupelClusterEnergyMod1 = new TNtuple("fNTupelClusterEnergyMod1", "Cluster Energy in Module 1", "energy");
		fNTupelClusterEnergyMod2 = new TNtuple("fNTupelClusterEnergyMod2", "Cluster Energy in Module 2", "energy");
		fNTupelClusterEnergyMod3 = new TNtuple("fNTupelClusterEnergyMod3", "Cluster Energy in Module 3", "energy");
	}

	// ********** END Create Histograms ********** //



	// ********** START Add Histograms to the Output ********** //

	fOutput->Add(fH1NEvents);
	fOutput->Add(fH1NEventsNamed);
	fOutput->Add(fH1NClusters);
	fOutput->Add(fH1NClustersPHOS);	
	fOutput->Add(fH1NClustersPHOSafterCuts);
	fOutput->Add(fH2NCellsPerClusterVsClusterEnergy);
	fOutput->Add(fH1Zvtx);
	fOutput->Add(fH1Mass);
	fOutput->Add(fH1MassMixed);
	fOutput->Add(fH1ClusterE);
	fOutput->Add(fH1ClusterEAfterCuts);
	fOutput->Add(fH2EAfterCutsVsModNum);
	fOutput->Add(fH1Pi0DecayModes);
	fOutput->Add(fH1Pi0DecayModesAfterCuts);
	fOutput->Add(fH1Pi0DecayModesIsPrimary);
	fOutput->Add(fH1Pi0DecayModesAddedPi0);
	fOutput->Add(fH1Pi0DecayModesAddedPi0PHOS);
	fOutput->Add(fH2NDaughtersVsPtSimPi0s);
	fOutput->Add(fH2NDaughtersVsdRSimPi0s);
	fOutput->Add(fH1MCpionVertDistToEventVert);
	fOutput->Add(fH1MCpionVertDistToEventVertLowDR);
	fOutput->Add(fH1MCpionVertDistToEventVertIsPrimary);
	fOutput->Add(fH1Pi0TruthPt);
	fOutput->Add(fH1K0Pi0TruthPt);
	fOutput->Add(fH1PriPi0TruthPt);
	fOutput->Add(fH1SecPi0TruthPt); 
	fOutput->Add(fH1K0Pi0TruthPtPhi2PiY05);
	fOutput->Add(fH1PriPi0TruthPtPhi2PiY05);
	fOutput->Add(fH1SecPi0TruthPtPhi2PiY05); 

	fOutput->Add(fH1Pi0TruthPtPhos);
	fOutput->Add(fH1K0Pi0TruthPtPhos);
	fOutput->Add(fH1PriPi0TruthPtPhos);
	fOutput->Add(fH1SecPi0TruthPtPhos); 
	fOutput->Add(fH1Pi0TruthPtPhi2PiY05);
	fOutput->Add(fH1Pi0TruthPtPhi2PiY05dR1);
	fOutput->Add(fH1PriPi0TruthPtPhi2PiY05dR1);
	fOutput->Add(fH1Pi0TruthPtPhi2piY03);

	fOutput->Add(fH2Pi0TruthPhiEta);
	fOutput->Add(fH1ElectronConversionR);
	fOutput->Add(fH1Pi0TruthPtPhotonsPhos);
	fOutput->Add(fH1K0Pi0TruthPtPhotonsPhos);  
	fOutput->Add(fH1PriPi0TruthPtPhotonsPhos);
	fOutput->Add(fH1SecPi0TruthPtPhotonsPhos); 

	fOutput->Add(fH2YVsPhiGenPi0);
	fOutput->Add(fH2YVsPhiGenEta);
	fOutput->Add(fH2YVsPhiAddedPi0);
	fOutput->Add(fH2YVsPhiAddedEta);
	fOutput->Add(fH2YVsPhiAddedPi0PHOS);
	fOutput->Add(fH2YVsPhiAddedEtaPHOS);
	fOutput->Add(fH1YAddedPi0);
	fOutput->Add(fH1EtaAddedPi0);
	fOutput->Add(fH1YAddedPi0PHOS);
	fOutput->Add(fH1EtaAddedPi0PHOS);
	fOutput->Add(fH1ClusterEAddedPi0);
	fOutput->Add(fH1ClusterEAddedPi0PHOS);
	
	fOutput->Add(fH1DecGammAddPi0Eta);
	fOutput->Add(fH1DecGammAddPi0Y);
	fOutput->Add(fH1DecGammAddPi0Phi);
	fOutput->Add(fH1DecGammAddPi0E);
	fOutput->Add(fH1DecGammInPHOSAddPi0E);
	fOutput->Add(fH1DecGammAddPi0Asymm);
	fOutput->Add(fH1DecGammAddPi0OpAngle);
	fOutput->Add(fH1DecGammAddPi0ConvR);
	fOutput->Add(fH1DecGammAddPi0ConvRate);
	fOutput->Add(fH1TruthPtAddedPi0GammasInPHOS);
	fOutput->Add(fH1DecGammAddPi0PHOSEta);
	fOutput->Add(fH1DecGammAddPi0PHOSY);
	fOutput->Add(fH1DecGammAddPi0PHOSPhi);
	fOutput->Add(fH1DecGammAddPi0PHOSE);
	fOutput->Add(fH1DecGammInPHOSAddPi0PHOSE);
	fOutput->Add(fH1DecGammAddPi0PHOSAsymm);
	fOutput->Add(fH1DecGammAddPi0PHOSOpAngle);
	fOutput->Add(fH1DecGammAddPi0PHOSConvR);
	fOutput->Add(fH1DecGammAddPi0PHOSConvRate);
	fOutput->Add(fH1TruthPtAddedPi0PHOSGammasInPHOS);
	
	fOutput->Add(fH1TruthPtGenPi0);
	fOutput->Add(fH1TruthPtGenEta);
	fOutput->Add(fH1TruthPtAddedPi0);
	fOutput->Add(fH1TruthPtAddedPi02PiY05);
	fOutput->Add(fH1TruthPtAddedPi0MesonInPHOS);
	fOutput->Add(fH1TruthPtAddedPi0PHOSMesonInPHOS);
	fOutput->Add(fH1TruthPtAddedEta);
	fOutput->Add(fH1TruthPtAddedPi0PHOS);
	fOutput->Add(fH1TruthPtAddedEtaPHOS);
	fOutput->Add(fH2MPtAddedPi0);
	fOutput->Add(fH2MPtAddedPi0_unweighed);
	fOutput->Add(fH2MPtAddedPi0PHOS);
	fOutput->Add(fH2MPtAddedPi0PHOS_unweighed);
	fOutput->Add(fH2MPtAddedPi0PHOSSingleContr);
	fOutput->Add(fH2MPtAddedPi0PHOSMultContr);	
	fOutput->Add(fH2MEnergyDiffAddedPi0PHOS);
	fOutput->Add(fH2PtRecVsPtTruthAddedPi0PHOS);
	fOutput->Add(fH2MPtAddedEta);
	fOutput->Add(fH2MPtAddedEtaPHOS);
	fOutput->Add(fH2ClusterEVSPhotonE);
	
	fOutput->Add(fH2HitmapEtaVsPhi);
	fOutput->Add(fH3MPtAsymm);
	fOutput->Add(fH3MPtModules);
	fOutput->Add(fH3MPtTiming);
	fOutput->Add(fH3MPtAsymmMix);
	fOutput->Add(fH3MPtModulesMix);
	fOutput->Add(fH3MPtTimingMix);
	fOutput->Add(fH2DphiDeta);
	fOutput->Add(fH2DphiDetaMix);
	fOutput->Add(fH1ClusterM02);
	fOutput->Add(fH2Ellipse);
	fOutput->Add(fH1ClusterDisp);
	fOutput->Add(fH2CellsM02);
	fOutput->Add(fH1Chi2);
	fOutput->Add(fH1NTrkMatch);
	
	fOutput->Add(fH1NPrimVertContribut);
	fOutput->Add(fH1DistPileUpPrimVert);
	fOutput->Add(fH1nSPDPileUpVtxs);
	fOutput->Add(fH1NClustersVsCuts);
	fOutput->Add(fTProfMeanClusterEnergyVsCuts);
	fOutput->Add(fH2CellIdVsE);
	fOutput->Add(fH2LocalMaxCellsIdVsE);
	fOutput->Add(fH2ClusterPosCellsIdVsE);
	fOutput->Add(fH1ClusterTOFWeightedWithE);
	fOutput->Add(fH2ClusterTOFVsE);
	fOutput->Add(fH1CellMCLabel); 
	fOutput->Add(fH1AppliedClusterCuts); 
	fOutput->Add(fH1DistanceToBadChannel); 
	fOutput->Add(fH2ClusterPositionsMod1);
	fOutput->Add(fH2ClusterPositionsMod2);
	fOutput->Add(fH2ClusterPositionsMod3);
	fOutput->Add(fH2ClusterPosMod1Gen1);
	fOutput->Add(fH2ClusterPosMod3Gen1);
	fOutput->Add(fH2ClusterPosMod1Gen2);
	fOutput->Add(fH2ClusterPosMod3Gen2);
	fOutput->Add(fNTupelClusterEnergyMod1);
	fOutput->Add(fNTupelClusterEnergyMod2);
	fOutput->Add(fNTupelClusterEnergyMod3);

	// ********** END Add Histograms to the Output ********** //

	// Post data for ALL output slots >0 here, 
	// To get at least an empty histogram 
	// 1 is the outputnumber of a certain weg of task 1  
	PostData(1, fOutput); 
}  

//________________________________________________________________________
void AliAnalysisTaskPHOSNeutralMeson::UserExec(Option_t *) {


	// ************* START Main Loop Called for each Event ************* //

	AliMCEvent *mcEvent = MCEvent();
	Bool_t isMC = bool(mcEvent);
	
	TRandom3 randy; randy.SetSeed(0);
	TLorentzVector ParentMix;


	// ************* START Load Event ************* //

	AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();    

	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
	if (!aodH && !esdH)  AliError("Could not get ESD or AODInputHandler");

	if (esdH){ 
		fAnyEv = dynamic_cast<AliESDEvent*> (esdH->GetEvent());
	}
	else if (aodH){ 
		fAnyEv = dynamic_cast<AliAODEvent*> (aodH->GetEvent());  
	}
	else{
		AliFatal("Neither ESD nor AOD event found");
		return;
	}

	fH1NEventsNamed->Fill(2); //All

	// ************* END Load Event ************* //


	// ************** START Remove Fast Cluster Events *************** //
   
	TString trigClasses = fAnyEv->GetFiredTriggerClasses();
			
	// remove non-CINT1 events
	if(fUseOnlyCINT1Events){
		if(!trigClasses.Contains("CINT1")){
			return;
		}
	}
		
	// remove "fast cluster events"
	if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")){
		return;
	}

	fH1NEventsNamed->Fill(3); //no FAST && !ALL

	// ************** END Remove Fast Cluster Events ***************** //

	// ************* START Remove Events that have no Primary Vertex  ************* //

	if (esdH) {
		//AliVVertex has no function GetStatus
		if (!(dynamic_cast<AliESDEvent*>(fAnyEv)->GetPrimaryVertex()->GetStatus())){
			return; 
		}
	}
	else if (aodH) {
		//use NContributors instead of GetStatus to check if there is a primary vertex
		if (!((dynamic_cast<AliAODEvent*>(fAnyEv)->GetPrimaryVertex()->GetNContributors())>0)) {
			return; 
		}
	}

	fH1NEventsNamed->Fill(4); //prim Vertex

	// ************* END Remove Events that have no Primary Vertex ************* //

	
	// ************* START Remove Events that are PileUpFromSPD ***************** //

	Float_t primVertZ, pileUpVertZ, nPrimVertContributors, nPileUpVtxices;
	TClonesArray* SPDPileupVertices;
   
	nPrimVertContributors = fAnyEv->GetPrimaryVertex()->GetNContributors();
	fH1NPrimVertContribut->Fill(nPrimVertContributors);
				
	if(fAnyEv->IsPileupFromSPD(3,.8,3.,2.,5.)) {//default values for AliESDEvent and AliAODEvent
		if (esdH) SPDPileupVertices = dynamic_cast<AliESDEvent*>(fAnyEv)->GetPileupVerticesSPD();
		else if (aodH) return; //reject AOD Pileup Events
		nPileUpVtxices = SPDPileupVertices->GetEntries();
		fH1nSPDPileUpVtxs->Fill(nPileUpVtxices);
		primVertZ   = fAnyEv->GetPrimaryVertex()->GetZ();
	  	
		for (Int_t kk = 0; kk<nPileUpVtxices; kk++) {
			pileUpVertZ = ((AliESDVertex*) SPDPileupVertices->At(kk))->GetZ();
			fH1DistPileUpPrimVert->Fill(primVertZ-pileUpVertZ);
		}
		return; //reject ESD Pile Up Events!
	}
	
	fH1NEventsNamed->Fill(5); //no Pileup from SPD

	// ************* END Remove Events that are PileUpFromSPD ***************** //

	// ********************** START zVertex Cut ********************** //
	
	Double_t vertZ;
	vertZ = fAnyEv->GetPrimaryVertex()->GetZ();  
	
	if(fDoZvertexCut){
		if(fUseIsVertexSelected2013pA) { // vertex cut has to be done different for p-Pb data
			fUtils->SetMaxVtxZ(fZvertexCut);
			if(!fUtils->IsVertexSelected2013pA(fAnyEv)) {return;}
		} 
		else {	
			if(fabs(vertZ)>fZvertexCut) {return;}
		}
	}

	fH1Zvtx->Fill(vertZ);
	fH1NEventsNamed->Fill(6); //ZVertex
	
	// ********************** END zVertex Cut ********************** //

	// ********************* START Setting PHOS matrix   ****************************** //

	Int_t runNumber = 0;
	runNumber = fAnyEv->GetRunNumber();

	if (fPHOSGeo==0) {

		fPHOSGeo = AliPHOSGeometry::GetInstance() ;

		if(!fPHOSGeo){ //Geometry not yet constructed with Tender
			fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP","");

			AliOADBContainer geomContainer("phosGeo");
			geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
			TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
			for(Int_t mod=0; mod<5; mod++) {
				if(!matrixes->At(mod)) continue;
				fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
				printf("....TASK.....Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
				((TGeoHMatrix*)matrixes->At(mod))->Print() ;
			}
		}
	}
	
	if(fEventCounter == 0) { // Only done for the first Event

		if(fFillCellIdVsE) {
			Int_t recoPass = -1;
			TTree * t = am->GetTree();
			if(t){  
				TFile * f = t->GetCurrentFile() ;
				if(f){  
					TString fname(f->GetName());
					if(fname.Contains("pass1"))
						recoPass=1;
					else 
						if(fname.Contains("pass2"))
							recoPass=2;
						else 
							if(fname.Contains("pass3")) 
								recoPass=3;
							else 
								if(fname.Contains("pass4")) 
									recoPass=4;
				}
			}
			if(recoPass<0){
				AliError("Can not find pass number from file name, is set to -1");
			}

			//Load recalibration data
			AliOADBContainer calibContainer("phosRecalibration");
			calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
			TObjArray *recalib = (TObjArray*)calibContainer.GetObject(runNumber,"PHOSRecalibration");
			if(!recalib){
				AliFatal(Form("Can not read calibrations for run %d.\n",runNumber)) ;
			}
			else{
				fPHOSCalibData = (AliPHOSCalibData*)recalib->At(recoPass-1) ;
				if(!fPHOSCalibData) {
					AliFatal(Form("Can not find calibration for run %d, pass %d \n",runNumber, recoPass)) ;
				}
			}
		} //if(fFillCellIdVsE)
	} //if(fEventCounter == 0)

	// ********************* END Setting PHOS matrix for ESDs ******************** //

	// **************** START fill nEvent histogram ****************** //

	// NUMBER OF EVENTS FOR NORMALIZING RAW YIELD:
	fH1NEvents->Fill(3); //All Events with PrimVertex and not IsPileupFromSPD
	fH1NEventsNamed->Fill(1); //good Events
	
	// **************** END fill nEvent histogram ********************* //
 


	//******************** START Cell QA histograms ************************ //
	
	if(fFillCellIdVsE) {
		AliVCaloCells *cells = fAnyEv->GetPHOSCells();
		Int_t nPHOSCell=0;
		nPHOSCell = cells->GetNumberOfCells();
		Float_t cellEnergy;
		Short_t cellNumber;
		
		for (Int_t c1=0; c1<nPHOSCell; c1++) {
			cellNumber = cells->GetCellNumber(c1);
			cellEnergy = cells->GetCellAmplitude(cellNumber);	
			
			if(fRecalibrateModuleWise){
				Double_t recalib[3] = {fRecalFactorMod1, fRecalFactorMod2, fRecalFactorMod3 }; //recalFactors are set in AddTask!
				Int_t    modNrThisCell, relId_thiscell[4];
				fPHOSGeo->AbsToRelNumbering(cellNumber,relId_thiscell);
				modNrThisCell  = relId_thiscell[0];
				cellEnergy = cellEnergy*recalib[modNrThisCell-1];
			}
			fH2CellIdVsE->Fill(cellNumber, cellEnergy);
		}
	}
	
	//******************* END Cell QA histograms ************************ //
   
   //**************************************************************************************
	//***************************** START MC part ****************************************** //

	if(isMC){
	if(fFillMCHistos){				
	
		// Counters for mc particle stack numbering  // store start and end of pythia particles, added signals,
		fIHijingMin = 0;	// first hijing particle index      (will remain 0)
		fIHijingMax = 0;	// last hijing particle index
		//fIPi0Min =0;   	   // first added pi0 (particle index)
		//fIPi0Max = 0; 		// last  added pi0 (particle index)
		//fIEtaMin = 0;      // first added eta (particle index)
		//fIEtaMax = 0; 	   // last  added eta (particle index)
		//fIPi0EMCMin = 0;	// first particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
		//fIPi0EMCMax = 0;	// last  particle in pi0EMC_3 (particle index)  (usually 1 pi0 and its 2 decay photons)
		//fIPi0PHOSMin = 0;	// first particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
		//fIPi0PHOSMax = 0;	// last  particle in pi0PHS_4 (particle index)  (usually 1 pi0 and its 2 decay photons)
		//fIEtaEMCMin = 0;	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)
		//fIEtaEMCMax = 0;	// first particle in etaEMC_5 (particle index)  (usually 1 eta and its 2 decay photons)	
		//fIEtaPHOSMin = 0;	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
		//fIEtaPHOSMax = 0;	// first particle in etaPHS_6 (particle index)  (usually 1 eta and its 2 decay photons)
		fNAddedSignalsTotal = 0;	
		
		
		if (!mcEvent){
			AliError("Event is not an MC event"); 
			return;
		}
		
		const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
		if (!evtVtx) return;

		//**************************************************************************************
		// *************** START Added Signal Stuff ********************************************
		if(fAnalyseAddedSignals) {
			
			//cout<<" ANALYSIING ADDED SIGNALS!"<<endl;
			
			TList *genHeaders = mcEvent->GetCocktailList();
			 
			if(!genHeaders) 	{
				AliError("ERROR: Could not retrieve genHeaders"); 
				return; 
			}
				
			Int_t nGenerators = genHeaders->GetEntries();

			//printf("N generators %d \n", nGenerators);

			// igen==1 and igen==2 are pi0 and eta added signals, respectively, when looking at lhc13ef
			// check eventHeader2->GetName(); for your dataset
			for(Int_t igen = 0; igen < nGenerators; igen++) {
				AliGenEventHeader* eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
				TString name = eventHeader2->GetName();
				Int_t nProduced = eventHeader2->NProduced();
				
				//cout<<"nProduced "<<nProduced<<"  name: "<<name<<endl;
				
				if (name.Contains("hijing",TString::kIgnoreCase) || name.Contains("DPMJET",TString::kIgnoreCase)) {
					fIHijingMax=nProduced-1;
				} else {
					fNAddedSignalsTotal += nProduced;
					
					////Not really needed anymore, because one can just use mcP->GetGeneratorIndex(); 
					//if (igen==1) {                   //name: pi0_1
						//fIPi0Min=fIHijingMax+1;
						//fIPi0Max=fIPi0Min+nProduced-1;
					//}
					//if (igen==2) {							//name: eta_2
						//fIEtaMin=fIPi0Max+1;
						//fIEtaMax=fIEtaMin+nProduced-1;
					//}
					//if (igen==3) {							//name: pi0EMC_3
						//fIPi0EMCMin=fIEtaMax+1;
						//fIPi0EMCMax=fIPi0EMCMin+nProduced-1;
					//}
					//if (igen==4) {							//name: pi0PHS_4
						//fIPi0PHOSMin=fIPi0EMCMax+1;
						//fIPi0PHOSMax=fIPi0PHOSMin+nProduced-1;
					//}
					//if (igen==5) {							//name: etaEMC_5
						//fIEtaEMCMin=fIPi0PHOSMax+1;
						//fIEtaEMCMax=fIEtaEMCMin+nProduced-1;
					//}
					//if (igen==6) {							//name: etaPHS_6
						//fIEtaPHOSMin=fIEtaEMCMax+1;
						//fIEtaPHOSMax=fIEtaPHOSMin+nProduced-1;
					//}
				}
			} //loop over generators
			
			// // Comment in to test and check numbering scheme !
			//Printf("Hijing        IDs: min = %d,   max = %d",fIHijingMin, fIHijingMax);  
			//Printf("Pi0 flat      IDs: min = %d, max = %d",fIPi0Min,    fIPi0Max);
			//Printf("Eta flat      IDS: m<"##"<<endl;
										//cout<<"##"<<endl;
										//cout<<in = %d, max = %d",fIEtaMin,    fIEtaMax);			
			//Printf("Pi0 EMC  flat IDS: min = %d, max = %d",fIPi0EMCMin, fIPi0EMCMax);			
			//Printf("Pi0 PHOS flat IDS: min = %d, max = %d",fIPi0PHOSMin,fIPi0PHOSMax);			
			//Printf("Eta EMC  flat IDS: min = %d, max = %d",fIEtaEMCMin, fIEtaEMCMax);			
			//Printf("Eta PHOS flat IDS: min = %d, max = %d",fIEtaPHOSMin,fIEtaPHOSMax);	
			//Printf("   fNAddedSignalsTotal = %d",fNAddedSignalsTotal);
			//Printf("fIEtaPHOSMax-fIPi0Min+1 = %d",(fIEtaPHOSMax - fIPi0Min +1));
					
		} //if(fAddedSignal)
		
		//*****************END  Added Signal Stuff *************************************************
		//**************************************************************************************
		//cout<<"weights: "<<fAddedSignalWeight[0]<<" "<<fAddedSignalWeight[1]<<" "<<fAddedSignalWeight[2]<<" "<<fAddedSignalWeight[3]<<" "<<endl;


		//mcEvent->PreReadAll();    
		Int_t nTracksMC  = mcEvent->GetNumberOfTracks();
		Int_t nPTracksMC = mcEvent->GetNumberOfPrimaries();
		

		if(fAnalyseAddedSignals) nPTracksMC -= fNAddedSignalsTotal;  //subtract number of AddedSignals from primary tracks
		
		// All Hijing generated + all Added Particles  = NumberOfPrimaties!
		Bool_t bAddPi0 = kFALSE;
		Bool_t bAddEta = kFALSE;
	
		//cout<<"nTracksMC "<<nTracksMC<<endl;
		//cout<<"nPTracksMC "<<nPTracksMC<<endl;
		//cout<<"fIHijingMax "<<fIHijingMax<<endl;
		//cout<<nPTracksMC<<" "<<fIHijingMax<<endl;
			
		for (Int_t iTrack = 0; iTrack<nTracksMC; ++iTrack) {
			
			Bool_t isPrimary     = kFALSE;      // in general particles from 0 to mcEvent->GetNumberOfPrimaries();
															// except in datasets with added signals. Then its mcEvent->GetNumberOfPrimaries() - fNAddedSignalsTotal
			Bool_t isSecondary  = kFALSE;		// in general the oposite of isPrimary
															// in datasets with addedSignals its all particles after the last added particle
			
			Bool_t isAddedParticle  = kFALSE;		
			
			// optimized for lhc13e7.   (its impossible to write this in genral. One needs to check the AliGenEventHeader
			// 									of the specific dataset and see what generators where used and what their names are)
			Bool_t isHijing		 = kFALSE;	//Particle from Hijing (or the min bias genearator in general)
			Bool_t isAddedPi01	 = kFALSE; 	//Pi0 added from the second AliGenEventHeader (pi0_1)					
			Bool_t isAddedPi0PHOS = kFALSE; 	//Pi0 added from the fifth AliGenEventHeader shooting at PHOS (pi0PHOS_4)
			Bool_t isAddedEta1	 = kFALSE; 	//Eta added from the third addtional AliGenEventHeader (eta_2)
			Bool_t isAddedEtaPHOS = kFALSE; 	//Eta added from the seventh addtional AliGenEventHeader shooting at PHOS (etaPHOS_6)
			Bool_t isAddedAtEmcal = kFALSE; 	//Eta or Pi0 added from the seventh addtional AliGenEventHeader shooting at EMCAL (pi0EMC_3 or etaEMC_5)
			Bool_t isK0sDecay     = kFALSE;
			Bool_t isMaterialSec  = kFALSE;


			// get particle from stack
			AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
			if (!mcP) continue;
			
			
			//Short_t generatorIndex = mcP->GetGeneratorIndex();
			//if(generatorIndex != 0) {
				////cout<<generatorIndex<<endl;
				////mcP->Print();
			//}
			
			//cout<<iTrack<<" - "<<mcP->Label()<<" = "<<iTrack-mcP->Label()<<endl;

				// Older implementation to skip non pi0s and etas
			//// Skip all particles that arent pi0s or etas
			//if(fAnalyseMCPi0 && fAnalyseMCEta) {  //(analyse both) 
				//if        (mcP->PdgCode() == 111) {		//pi0
				//} 	else if(mcP->PdgCode() == 221) {		//eta
				//} 	else
					//continue;
			//// OR: Skip all particles that arent pi0s
			//} else if(fAnalyseMCPi0 && !fAnalyseMCEta) {
				//if        (mcP->PdgCode() == 111) {		//pi0
				//} 	else
					//continue;
			//// OR: Skip all particles that arent etas
			//} else if(!fAnalyseMCPi0 && fAnalyseMCEta) {
				//if        (mcP->PdgCode() == 221) {		//eta
				//} 	else
					//continue;
			//}
			
			
			//Bool used to skip particles that arent pi0 or eta
			Bool_t particleIsPi0OrEta = kFALSE;
			
			// If analysing MCPi0s and particle is pi0 -> keep particle
			if(fAnalyseMCPi0) {
				if(mcP->PdgCode() == 111) {
					particleIsPi0OrEta = kTRUE;	//pi0
				}
			}
			// If analysing MCEtas and particle is eta -> keep particle
			if(fAnalyseMCEta) {
				if(mcP->PdgCode() == 221) {
					particleIsPi0OrEta = kTRUE;	//eta
				}
			}
			//Skip all particles that didnt meet above requirement									
			if(!particleIsPi0OrEta) continue;
			
			
			Short_t generatorIndex = mcP->GetGeneratorIndex();
			Int_t mcLabel = mcP->Label();			
			Int_t pdgCode = mcP->PdgCode();
			
			//nProduced 151  name: hijing_0
			//nProduced 1  name: pi0_1
			//nProduced 1  name: eta_2
			//nProduced 3  name: pi0EMC_3
			//nProduced 3  name: pi0PHS_4
			//nProduced 3  name: etaEMC_5
			//nProduced 2  name: etaPHS_6


			if(fAnalyseAddedSignals) {
				// is it generator particle or added signal?
				isHijing = kTRUE;
				if(mcLabel > fIHijingMax){
					isHijing = kFALSE;
					
					if(generatorIndex == 1 && pdgCode == 111) isAddedPi01    = kTRUE;  //pi0_1
					if(generatorIndex == 2 && pdgCode == 221) isAddedEta1    = kTRUE;  //eta_2
					if(generatorIndex == 3 && pdgCode == 111) isAddedAtEmcal = kTRUE;  //pi0EMC_3
					if(generatorIndex == 4 && pdgCode == 111) isAddedPi0PHOS = kTRUE;  //pi0PHS_4
					if(generatorIndex == 5 && pdgCode == 221) isAddedAtEmcal = kTRUE;  //etaEMC_5
					if(generatorIndex == 6 && pdgCode == 221) isAddedEtaPHOS = kTRUE;  //etaPHS_6
					
					//if(mcLabel <= fIPi0Max){
						//isAddedPi01 = kTRUE;
						//////Printf("Found the added pi0.  ID:  %d   pdg: %d        generatorIndex: %d",mcLabel,mcP->PdgCode(),generatorIndex); 
					//} else if(mcLabel > fIPi0Max && mcLabel <= fIEtaMax){
						//isAddedEta1 = kTRUE;
						////Printf("Found the added eta.  ID: %d   pdg: %d     generatorIndex: %d",mcLabel,mcP->PdgCode(),generatorIndex);
					//} else if(mcLabel > fIEtaMax && mcLabel <= fIPi0EMCMax){
					   //isAddedAtEmcal = kTRUE;
					//} else if(mcLabel > fIPi0EMCMax && mcLabel <= fIPi0PHOSMax){
						//isAddedPi0PHOS = kTRUE;
						////Printf("Found the added PHOSpi0.  ID:  %d   pdg: %d    generatorIndex: %d",mcLabel,mcP->PdgCode(),generatorIndex); 
					//} else if(mcLabel > fIPi0PHOSMax && mcLabel <= fIEtaEMCMax){
						//isAddedAtEmcal = kTRUE;
					//} else if(mcLabel > fIEtaEMCMax && mcLabel <= fIEtaPHOSMax){
						//isAddedEtaPHOS = kTRUE;
						////Printf("Found the added PHOSeta.  ID: %d   pdg: %d    generatorIndex: %d",mcLabel,mcP->PdgCode(),generatorIndex); 
						////Printf("Found the added PHOSeta.  ID: %d   pdg: %d    generatorIndex: %d",mcLabel,mcP->PdgCode(),generatorIndex); 
					//}
				}
				
				//Usefull to study which particles are not included in the first AliGenEventHeader but are still Primary
				//for DPMJET they are pi0s and etas from strong deacys
				//for hijing the strong decays are in the first AliGenEventHeader
				//if(!isHijing && !bAddPi0 && !bAddEta && !isAddedEtaPHOS && !isAddedEtaPHOS && isPrimary) {
					//if(mcP->GetMother()>-1){
					//Int_t pdgcode = mcP->PdgCode();
					//Int_t pdgcodeMother = ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode();
					////cout<<"Particle that is !bGen && isPrimary pdgCode: "<<pdgcode<<"  mother "<<pdgcodeMother<<endl;
					////mainly pi0s from  221 213 223 -213 (strong decays)
					//}
				//}
								
				//Fill Histos for added signals	
				if(isHijing) {
					if(pdgCode== 111) {
						fH2YVsPhiGenPi0->Fill(mcP->Phi(),mcP->Y());
						fH1TruthPtGenPi0->Fill(mcP->Pt());
					} else if (pdgCode == 221) {
						fH2YVsPhiGenEta->Fill(mcP->Phi(),mcP->Y());
						fH1TruthPtGenEta->Fill(mcP->Pt());
					}	
				}
				Float_t wgtPi0 = CalculateAddedSignalWeightPi0Exp(mcP->Pt());
				
				// ****** ****** ****** Pi0 Generator wide range ****** ****** ***** ****** 
				if(isAddedPi01) {
					//cout<<"mcLabel label von pi0_1 pion "<<mcLabel<<endl;
					fH2YVsPhiAddedPi0->Fill(mcP->Phi(),mcP->Y());
					fH1TruthPtAddedPi0->Fill(mcP->Pt(),wgtPi0);
					if(mcP->Y()>-0.5 && mcP->Y()<0.5) {
						fH1TruthPtAddedPi02PiY05->Fill(mcP->Pt(),wgtPi0);
					}
					
					FillHistoWithDaughterInfoPi0(mcP, fH1Pi0DecayModesAddedPi0, mcEvent);
					 
					
					Double_t mcAddedPEta = mcP->Eta();
					Double_t mcAddedPhi  = mcP->Phi();
								
					if(mcAddedPEta >= fEtaAccMin && mcAddedPEta <= fEtaAccMax && mcAddedPhi >= fPhiAccMin && mcAddedPhi <= fPhiAccMax) {
						fH1TruthPtAddedPi0MesonInPHOS->Fill(mcP->Pt(),wgtPi0);
						
						if(fFillDecayPhotonInfoAddedSig) {
						FillDecayGammaHistos( mcP, mcEvent, wgtPi0,
													fH1DecGammAddPi0Eta,
													fH1DecGammAddPi0Y,
													fH1DecGammAddPi0Phi,
													fH1DecGammAddPi0E,
													fH1DecGammInPHOSAddPi0E,
													fH1DecGammAddPi0Asymm,
													fH1DecGammAddPi0OpAngle,
													fH1DecGammAddPi0ConvR,
													fH1DecGammAddPi0ConvRate,
													fH1TruthPtAddedPi0GammasInPHOS);
						}
					}
					
					fH1YAddedPi0->Fill(mcP->Y());
					fH1EtaAddedPi0->Fill(mcP->Eta());
					
					
				}
				// ****** ****** ****** ****** ****** ****** ****** ****** ****** ****** 
				
				if(isAddedEta1) {
					fH2YVsPhiAddedEta->Fill(mcP->Phi(),mcP->Y());
					fH1TruthPtAddedEta->Fill(mcP->Pt());
				}
				
				// ****** ****** ****** Pi0 Generator directly at PHOS ****** ****** ****** 
				if(isAddedPi0PHOS) {
					fH2YVsPhiAddedPi0PHOS->Fill(mcP->Phi(),mcP->Y());
					fH1TruthPtAddedPi0PHOS->Fill(mcP->Pt(),wgtPi0);
					
					FillHistoWithDaughterInfoPi0(mcP, fH1Pi0DecayModesAddedPi0PHOS, mcEvent);

					
					Double_t mcAddedPEta = mcP->Eta();
					Double_t mcAddedPhi  = mcP->Phi();
								
					if(mcAddedPEta >= fEtaAccMin && mcAddedPEta <= fEtaAccMax && mcAddedPhi >= fPhiAccMin && mcAddedPhi <= fPhiAccMax) {
						fH1TruthPtAddedPi0PHOSMesonInPHOS->Fill(mcP->Pt(),wgtPi0);
						
						if(fFillDecayPhotonInfoAddedSig) {
						FillDecayGammaHistos( mcP, mcEvent, wgtPi0,
													fH1DecGammAddPi0PHOSEta,
													fH1DecGammAddPi0PHOSY,
													fH1DecGammAddPi0PHOSPhi,
													fH1DecGammAddPi0PHOSE,
													fH1DecGammInPHOSAddPi0PHOSE,
													fH1DecGammAddPi0PHOSAsymm,
													fH1DecGammAddPi0PHOSOpAngle,
													fH1DecGammAddPi0PHOSConvR,
													fH1DecGammAddPi0PHOSConvRate,
													fH1TruthPtAddedPi0PHOSGammasInPHOS);
						}
					}
					
					fH1YAddedPi0PHOS->Fill(mcP->Y());
					fH1EtaAddedPi0PHOS->Fill(mcP->Eta());
					
					
				}
				// ****** ****** ****** ****** ****** ****** ****** ****** ****** ****** 
				
				if(isAddedEtaPHOS) {
					fH2YVsPhiAddedEtaPHOS->Fill(mcP->Phi(),mcP->Y());
					fH1TruthPtAddedEtaPHOS->Fill(mcP->Pt());
				}
				
				if(isAddedEta1 || isAddedEtaPHOS || isAddedAtEmcal || isAddedPi01 || isAddedPi0PHOS) {
					isAddedParticle = kTRUE;
				}
			} //if(fAnalyseAddedSignals)
			
			if(!isAddedParticle) {
						//check if its a primary particle
				if(iTrack<nPTracksMC)  isPrimary   = kTRUE;
				else                   isSecondary = kTRUE;
			}
			
			//if(!(isPrimary^isSecondary^isAddedParticle)) cout<<"isPrimary XOR isSecondary XOR isAddedParticle is false. something is wrong. (One and only one must be true.)"<<endl;
			
			//TODO: Either implement a way to set pdg code (eta or pi0) from AddTask
			//		  or implement all truthPt histos for eta as well an fill both etas and pi0 using if conditions
			if (mcP->PdgCode() != 111) continue; //also skip etas now. (keep only pi0)
			
			if (isAddedParticle) continue; // skip all added particles from here.
							
							
			//check if Particle is a K0 decays
			isK0sDecay = 0;
			if(mcP->GetMother()>-1) {
				if( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  310 ||
						((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -310  )
					isK0sDecay = 1;
			}
			
			// and it's close enough to the event vertex
			//Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) +
					//(mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
					
					
			Double_t dR = TMath::Sqrt(	(mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) +
												(mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()) +
												(mcP->Zv()-evtVtx->GetZ())*(mcP->Zv()-evtVtx->GetZ())
											);
					
			//TVector3 v;
			//v.SetXYZ(mcP->Xv(),mcP->Yv(),mcP->Zv());
	      //TVector3 v0;
			//v0.SetXYZ(evtVtx->GetX(),evtVtx->GetY(),evtVtx->GetZ());
			//TVector3 v0v = v-v0;
			//cout<<dR<<"  "<<v0v.Mag()<<endl;
			
			
			fH1MCpionVertDistToEventVert->Fill(dR);
			fH1MCpionVertDistToEventVertLowDR->Fill(dR);
			
			if(isPrimary) fH1MCpionVertDistToEventVertIsPrimary->Fill(dR);
			
			// fFillFeedDownHistos is not just for filling feed down, but neccessary for all mc anaylses at this point
			// this bool is not implemented well, needs fixing.  (2. Jan 2016)
			//if(fFillFeedDownHistos) {
		
			Int_t daughter[2] = {-1,-1};
			daughter[0] = mcP->GetDaughterFirst();
			daughter[1] = mcP->GetDaughterLast();
			
			
			// ************************** For not added pi0, check how they decay ***************************
		
			FillHistoWithDaughterInfoPi0(mcP, fH1Pi0DecayModes, mcEvent);
			if(isPrimary) FillHistoWithDaughterInfoPi0(mcP, fH1Pi0DecayModesIsPrimary, mcEvent);
			
			Int_t nDaughters = GetNDaughtersOfMCParticle(mcP);				
			
			if(fFillNDaughtersVsPtAndR) {
				fH2NDaughtersVsPtSimPi0s->Fill(nDaughters,mcP->Pt());
				fH2NDaughtersVsdRSimPi0s->Fill(nDaughters,dR);
			}
			
			Bool_t mcPIsPhysicalDecay = MCParticleIsPhysicalDecay(mcP, mcEvent);
			
			if(!mcPIsPhysicalDecay) continue; //skip the particles rejected by MCParticleIsPhysicalDecay
			
					
	      
			//if( v0v.Mag()>1 ) continue;
			
			//if(dR > 1) continue; //skip particles too far away from event vertex.
											//different primary definition

			// Make sure that the pion has 2 daughters
			//if (daughter[0]<0)  {
				//continue;			// if the pion does not decay...
			//}
			//if (daughter[1]<0) {
				//daughter[1]=daughter[0];     // if the pion only has one daughter ...?   
			//}
			//if (daughter[1]-daughter[0] != 1) {
			  //continue;	// if the Pion has more than 2 daughters (Dalitz-decay)
			//}
			
			
			// ******************* from HERE ******************************************
			//      this part of code only works if only two daughters are accepted
			//			now also 0, 1 or 3 are possible		  
			//			parts of this code where for investigating converted photons
			//			parts where for checking if all decay products where going towards detector acceptance
			//
			
			//Int_t eIndexofConvertedPhoton[2] = {-1,-1};

			//bool bacc = true; 	
			//bool binp = true;
			//Double_t eta_d[2] = {0.0,0.0};
			//Double_t phi_d[2] = {0.0,0.0};

			//for (Int_t daughter_index=0; daughter_index<nDaughters; daughter_index++){
			
				//const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(daughter[daughter_index]));
				//eta_d[daughter_index] = dmc->Eta();
				//phi_d[daughter_index] = dmc->Phi();
		  
				//if(!(dmc->PdgCode()==22))	  binp = false;
				//if(!(dmc->PdgCode()==22 && 
				//eta_d[daughter_index]>fEtaAccMin && eta_d[daughter_index]<fEtaAccMax && 
				//phi_d[daughter_index]>fPhiAccMin && phi_d[daughter_index]<fPhiAccMax))   bacc = false;	
		  
				//if(dmc->GetFirstDaughter()>0 && dmc->GetLastDaughter()>0) {
					
					//// get the photons's daughters... 
					//const AliMCParticle *dmcd1 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetFirstDaughter()));
					//const AliMCParticle *dmcd2 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetLastDaughter()));
					//Double_t productionR1 = TMath::Sqrt(dmcd1->Xv()*dmcd1->Xv() + dmcd1->Yv()*dmcd1->Yv());
					//if(bacc)  fH1ElectronConversionR->Fill(productionR1);
					
					//// check if this is a conversion... 
					//if( (dmcd1->PdgCode()== -1.0*dmcd2->PdgCode()) &&
					//(dmcd1->PdgCode()==11 || dmcd1->PdgCode()==-11) &&
					//productionR1<460.0){ //460 PHOS distance to beamline 
						
						////find the conv e with highest energy, assign it to be that photon decay product.
						//if( dmcd1->E() > dmcd2->E() )
							//eIndexofConvertedPhoton[daughter_index] = dmc->GetFirstDaughter();
						//else
							//eIndexofConvertedPhoton[daughter_index] = dmc->GetLastDaughter();
					//}
				//}
			//}
			// ***************** until HERE **************************************************
			
			
			//if(binp!=true) continue; //Skip particles that did not decay into 2 photons
			//if(binp==true) 
			
			isMaterialSec = 0;

			if(isSecondary){
				
				if(((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2212 || //proton
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2212 || //anti-proton
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2112 || //neutron
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2112 || //anti-neutron
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  321  || //K+
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -321  || //K-
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  211  || //pi+
					((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -211     //pi-
				  ) 
				  isMaterialSec = 1;
			}
			
			
			fH2Pi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta()); 
			fH1Pi0TruthPt->Fill(mcP->Pt());
			if(isK0sDecay)  		fH1K0Pi0TruthPt->Fill(mcP->Pt());
			if(isPrimary)	   	fH1PriPi0TruthPt->Fill(mcP->Pt());
			if(isSecondary && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPt->Fill(mcP->Pt());
			
			
			FillHistoWithDaughterInfoPi0(mcP, fH1Pi0DecayModesAfterCuts, mcEvent);	
			
			
			// If you get "floating point exeption" in the testtrain use a newer version of aliroot there might be a bug for pythia
			if(mcP->Y()>-0.5 && mcP->Y()<0.5){
														fH1Pi0TruthPtPhi2PiY05->Fill(mcP->Pt()); 
				if(isK0sDecay)  					fH1K0Pi0TruthPtPhi2PiY05->Fill(mcP->Pt());
				if(isPrimary)   					fH1PriPi0TruthPtPhi2PiY05->Fill(mcP->Pt());
				if(dR < 1)							fH1Pi0TruthPtPhi2PiY05dR1->Fill(mcP->Pt());
				if(dR < 1 && isPrimary)			fH1PriPi0TruthPtPhi2PiY05dR1->Fill(mcP->Pt());
				
				if(isSecondary && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPtPhi2PiY05->Fill(mcP->Pt());
			}
			
			
			if(mcP->Y()>-0.3 && mcP->Y()<0.3) {
				fH1Pi0TruthPtPhi2piY03->Fill(mcP->Pt()); 
			}
			
			Bool_t mesonInAcceptance = kFALSE;
			Double_t mcPEta = mcP->Eta();
			Double_t mcPPhi = mcP->Phi();
						
			if(mcPEta >= fEtaAccMin && mcPEta <= fEtaAccMax && mcPPhi >= fPhiAccMin && mcPPhi <= fPhiAccMax) {
				mesonInAcceptance = kTRUE;
			}
			
			if(mesonInAcceptance) {
											fH1Pi0TruthPtPhos->Fill(mcP->Pt());
				if(isK0sDecay)  		fH1K0Pi0TruthPtPhos->Fill(mcP->Pt());
				if(isPrimary)	   	fH1PriPi0TruthPtPhos->Fill(mcP->Pt());
				if(isSecondary && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPtPhos->Fill(mcP->Pt());
			}
				
			//}  //fFillFeedDownHistos)
		}//for(nTracksMC) (loop over all MC particles) 
	}//if(fFillMCHistos)
	}//if(isMC)


	//***************************** END MC part ****************************************** //
	//**************************************************************************************


	// ***************** START Declaration of Cluster Variables *************** //

	Int_t nclusters = 0;
	nclusters = fAnyEv->GetNumberOfCaloClusters(); //number of all (EMCAL and PHOS) clusters
	fH1NClusters->Fill(nclusters); //nclusters (EMCAL and PHOS) after event cuts

	Int_t izvtx = GetZvtxBin(vertZ);		 //Important for event Mixing
	Int_t imult = GetMultBin(nclusters); //Important for event Mixing
	TLorentzVector Photon, Parent;
	
	Double_t vertex[3];
	if (esdH) dynamic_cast<AliESDEvent*>(fAnyEv)->GetVertex()->GetXYZ(vertex);
	else if (aodH) dynamic_cast<AliAODEvent*>(fAnyEv)->GetVertex(0)->GetXYZ(vertex); // function GetVertex needs argument for AODs
	
	Double_t clusterE;

	// ***************** END Declaration of Cluster Variables *************** //



	// ***************** START Cluster loop for data (ESD/AOD) *************** //

	Int_t nclustersPHOS = 0;
	Int_t nclustersPHOSafterCuts = 0; 

	for(Int_t i=0; i<nclusters; i++) {

		AliVCluster*	virCluster;
		if (esdH)	  	virCluster = dynamic_cast<AliESDCaloCluster*>(fAnyEv->GetCaloCluster(i)); //pointer to Phos cluster
		else if (aodH)	virCluster = dynamic_cast<AliAODCaloCluster*>(fAnyEv->GetCaloCluster(i)); //pointer to Phos cluster
		if(!virCluster) { 
			AliError(Form("ERROR: Could not retrieve any (ESD or AOD) Cluster %d",i)); 
			continue; 
		}
		
		if(virCluster->IsPHOS()) {		
			
			fH1NClustersVsCuts->Fill(1); 
			fTProfMeanClusterEnergyVsCuts->Fill(1,virCluster->E());
			
			// Reasons why Clusters can be cut away:
			// Number of cells per cluster below defined value	(fClusterMinCells)
			// Cluster energy below defined value				(fClusterMinE)
			// Cluster M02 below defined value					(fClusterMinM02)
			// Position of cluster in badCell					(clusterPosBad = true)
			// Cluster-Tof is not within defined window			(badTiming = true)

			// All Cuts are false on default. (Set to true if cut is in use and applies to current cluster)
			Bool_t clusterPosBad = false;
			Bool_t badTiming = false;


			// ************* START Definition of cell variables **************** //

			Float_t cellEnergy = 0.0; 
			Float_t maxEnergy = 0.0;
			Int_t maxID = -1;

			Int_t cellMCLabel = -1;
			UShort_t* CellsID = virCluster->GetCellsAbsId();
			
			// ************* END Definition of cell variables **************** //
			
			
			
			// *************** START Getting Cluster Relative ID ****************** //
			
			Float_t pos[3] = {0,0,0};
			virCluster->GetPosition(pos);
			TVector3 vpos(pos);

			Int_t    modNrClusterPos, relId[4], cellX, cellZ; //cellAbsId
			fPHOSGeo->GlobalPos2RelId(vpos,relId);    
			modNrClusterPos  = relId[0];
			cellX = relId[2];
			cellZ = relId[3];
			
			// *************** END Getting Cluster Relative ID ****************** //
			
			
			
			// *************** START Cluster Recalibration ****************** //
			
			clusterE = virCluster->E();
			
			if(fRecalibrateModuleWise) {
				Double_t recalib[3] = {fRecalFactorMod1, fRecalFactorMod2, fRecalFactorMod3 };  //recalFactors are set in AddTask!
				virCluster->SetE(virCluster->E()*recalib[modNrClusterPos-1]);
			}
			
			if(fDoClusterEnergyRecalibration) {
				clusterE = RecalibratePHOSClusterEnergy(fRecalibrationOption, clusterE, runNumber); 
				virCluster->SetE(clusterE); 
			}
							
			if(isMC && fDoPeakSmearing) {
				TF1 *f1_peak_width = new TF1("f1_peak_width", "gaus", 0., 2.);
				f1_peak_width->SetParameters(1., 1., fSmearFactor);
				Float_t new_cluster_energy = f1_peak_width->GetRandom();
				clusterE = clusterE*new_cluster_energy; 
				virCluster->SetE(clusterE); 
			}				
			
			// *************** END Cluster Recalibration ****************** //

		
	
			// ************* START Fill Cell MC Label *********************** //
			
			if (isMC) {
				AliVCaloCells *vcells = fAnyEv->GetPHOSCells();
				
				for(Int_t kk=0; kk < virCluster->GetNCells(); kk++) {

					cellMCLabel = vcells->GetCellMCLabel(CellsID[kk]); 
					fH1CellMCLabel->Fill(cellMCLabel);
				}
			}
			
			// ************* END Fill Cell MC Label ************************* //

		
		
			// ************* START Timing Cut *********************************** //

			Float_t clusterTOF;
			clusterTOF   = virCluster->GetTOF();
			
			if(fDoTimingCut) {

				// Check Cluster Timing 
				if(clusterTOF < fTimingCutMin || clusterTOF > fTimingCutMax) {   
					badTiming = true;
				}    
			}

			// *************** END timing cut *********************************** //


					
			// *************** START Cluster Position/Bad Channel Cut (Bad Map) ****************** //

			if (modNrClusterPos < 1 || modNrClusterPos > 3) {
				AliError(Form("Wrong module number %d",modNrClusterPos));
				continue;
			}
			
			if(fApplyBadMapManually) {
				if(!fDoDistToBadCellCutOnCellLevel){	// default: no distance-to-BadCell-Cut
					if(!IsGoodChannel("PHOS",modNrClusterPos,cellX,cellZ)) {
						clusterPosBad = true;
					}
			  	}
				else {
			    	if(!IsGoodChannelDistanceOnCellLevel("PHOS",modNrClusterPos,cellX,cellZ,fDistToBadCellOnCellLevel)){
						clusterPosBad = true;
			    	}
				}
			}
			
			if (!fDoDistToBadCellCut) {
				fDistToBadCell = 0.;  // cut is disabled by setting value to zero
			}

			// *************** END Cluster Position/Bad Channel Cut (Bad Map) ****************** //
	   


			// *************** START Fill Cell Id vs Energy histogram (needed for Bad Maps) *************************** //
			
			if(fFillCellIdVsE) {
			
				// Find Max Contributing Cell
				AliVCaloCells *vcells = fAnyEv->GetPHOSCells();
				
				for(Int_t kk=0; kk < virCluster->GetNCells(); kk++) {
					cellEnergy = vcells->GetCellAmplitude(CellsID[kk]); 
					if(cellEnergy > maxEnergy) {
						maxEnergy = cellEnergy;
						maxID = CellsID[kk];
					}	     	      	   
				}
			
				// Cell ID vs. energy for cells that are at the position of a cluster
				Int_t cellClusterPosAbsID;
				fPHOSGeo->RelToAbsNumbering(relId, cellClusterPosAbsID); //Get AbsID of Cell in which the ClusterPosition is
				Double_t cellClusterPosE = vcells->GetCellAmplitude(cellClusterPosAbsID); //Get Energy of that cell
				
				// Cell ID vs. energy for cells that are the maximum of a cluster
				if(fRecalibrateModuleWise) {
					Double_t recalib[3] = {fRecalFactorMod1, fRecalFactorMod2, fRecalFactorMod3}; //recalFactors are set in AddTask!
					
					cellClusterPosE = cellClusterPosE*recalib[modNrClusterPos-1];
					
					Int_t    modNrMaxCell, relId_maxcell[4];
					fPHOSGeo->AbsToRelNumbering(maxID,relId_maxcell);
					modNrMaxCell  = relId_maxcell[0];
					maxEnergy = maxEnergy*recalib[modNrMaxCell-1];
				}

				fH2LocalMaxCellsIdVsE->Fill(maxID, maxEnergy);
				
				fH2ClusterPosCellsIdVsE->Fill(cellClusterPosAbsID,cellClusterPosE);
			}
			
			// *************** END Fill Cell Id vs Energy histogram (needed for Bad Maps) *************************** //



			// *************** START Fill Cluster shape, matched tracks, energy and number of cells per cluster histograms *************************** //

			if(virCluster->GetM02() != 0) fH1ClusterM02->Fill(virCluster->GetM02());	      
			fH2Ellipse->Fill(virCluster->GetM20(),virCluster->GetM02());
			fH1ClusterDisp->Fill(virCluster->GetDispersion());
			fH2CellsM02->Fill(virCluster->GetNCells(),virCluster->GetM02());
			fH1Chi2->Fill(virCluster->Chi2());
			fH1NTrkMatch->Fill(virCluster->GetNTracksMatched());

			fH1ClusterE->Fill(virCluster->E());

			fH2NCellsPerClusterVsClusterEnergy->Fill(virCluster->E(),virCluster->GetNCells());
			
			// *************** END Fill Cluster shape, matched tracks, energy and number of cells per cluster histograms *************************** //
		
		
			fH1NClustersVsCuts->Fill(2.);  
			fTProfMeanClusterEnergyVsCuts->Fill(2,virCluster->E());
			
			nclustersPHOS++;  
		

			// ********* START Applying cuts ********** //

			if (!(virCluster->GetNCells()<fClusterMinCells)) {
				fH1NClustersVsCuts->Fill(3.);
				fTProfMeanClusterEnergyVsCuts->Fill(3,virCluster->E());
				
				
				if (!(virCluster->E()<fClusterMinE)) {
					fH1NClustersVsCuts->Fill(4.);
					fTProfMeanClusterEnergyVsCuts->Fill(4,virCluster->E());
					
					
					if (!(virCluster->GetM02()<fClusterMinM02)) {
						fH1NClustersVsCuts->Fill(5.);	
						fTProfMeanClusterEnergyVsCuts->Fill(5,virCluster->E());
						
												
						if (!clusterPosBad) { 
							fH1NClustersVsCuts->Fill(6.);
							fTProfMeanClusterEnergyVsCuts->Fill(6,virCluster->E());

							fH1DistanceToBadChannel->Fill(virCluster->GetDistanceToBadChannel());
							
							if(virCluster->GetDistanceToBadChannel() > fDistToBadCell) {    
								fH1NClustersVsCuts->Fill(7.);
								fTProfMeanClusterEnergyVsCuts->Fill(7,virCluster->E());
								
								if(fFillTimingHistos) {
									fH1ClusterTOFWeightedWithE->Fill(clusterTOF, virCluster->E());
									fH2ClusterTOFVsE->Fill(clusterTOF, virCluster->E()); 
								}
												
								if (!badTiming) {
									fH1NClustersVsCuts->Fill(8.);
									fTProfMeanClusterEnergyVsCuts->Fill(8,virCluster->E());
						
									// ********* END Applying cuts ********** //

									// **************************** START Good PHOS Clusters ********************************************** //

									nclustersPHOSafterCuts++; 
					
									//Necessary for Added Signal Analysis
									//These can only be set to true if fAnalyseAddedSignals is True
									//Used later to fill cluster in correct fPhotonsAdded vector.
									//Only clusters from the same added particle can be mixed
									Bool_t motherIsHijing 		 = kFALSE;  // particle from hijing (or min bias generator in general)
									Bool_t motherIsAddedPi01	 = kFALSE; 	//Pi0 added from the second AliGenEventHeader (pi0_1)					
									Bool_t motherIsAddedPi0PHOS = kFALSE; 	//Pi0 added from the fifth AliGenEventHeader shooting at PHOS (pi0PHOS_4)
									Bool_t motherIsAddedEta1	 = kFALSE; 	//Eta added from the third addtional AliGenEventHeader (eta_2)
									Bool_t motherIsAddedEtaPHOS = kFALSE; 	//Eta added from the seventh addtional AliGenEventHeader shooting at PHOS (etaPHOS_6)
									Bool_t motherIsAddedAtEmcal = kFALSE; 	//Eta or Pi0 added from the seventh addtional AliGenEventHeader shooting at EMCAL (pi0EMC_3 or etaEMC_5)
									Bool_t motherIsFromAddedSignal  = kFALSE; 	//
									Bool_t clusHasSingleContrib = kFALSE;
									
									
									// ************* START Added Signal Analysis - Determine if cluster is from Added Signal *****************
									
									//
									Float_t clusterEminusParticleE = 0.0;
									
									if(fAnalyseAddedSignals) {
										
										// The label in the AOD is the index in the TClonesArray of AliAODMCParticles
										Int_t ilabel = -1;
										ilabel = virCluster->GetLabel();    //CusterLabel = MCParticle-Index
										
										//cout<<"##"<<endl;
										//cout<<"##"<<endl;
										//cout<<"#######"<<endl;
										
										//for(Int_t i = 0; i < virCluster->GetNLabels(); i++) {
											//AliAODMCParticle* particle = (AliAODMCParticle*)mcEvent->GetTrack((virCluster->GetLabels())[i]);
											//cout<<"         Other contributing primary: "<<particle->GetName()<< "; Energy "<<particle->E()<<endl;
										//}
										
										if(virCluster->GetNLabels() == 1) clusHasSingleContrib = kTRUE;
										
										//if(clusHasSingleContrib) cout<<"clusHasSingleContrib"<<endl;
										//if(!clusHasSingleContrib) cout<<"MULTIContrib"<<endl;
										
										
										//if(ilabel == -1) cout<<"as los?"<<en~dl;
										//cout<<"LABEL "<<ilabel<<eqndl;
										//if(0){
										if(ilabel != -1) {
											
											//get mcParticle that made the cluster
											AliMCParticle *mcPCluster = static_cast<AliMCParticle*>(mcEvent->GetTrack(ilabel)); 
											
											clusterEminusParticleE = virCluster->E() - mcPCluster->E();
											
											//cout<<mcPCluster->PdgCode()<<endl;
											
											if(mcPCluster->PdgCode() == 22 ) { // cluster is from photon
												
												fH2ClusterEVSPhotonE->Fill(mcPCluster->E(),virCluster->E());
											
												if(mcPCluster->GetMother() > 0) { //there is a mother
													
													//get mother
													AliMCParticle *mcPMother = ((AliMCParticle*)mcEvent->GetTrack(mcPCluster->GetMother()));
													
													Bool_t motherPi0OrEta = kFALSE;
													
													// Skip all particles that arent pi0s or etas
													if(fAnalyseMCPi0) {
														if(mcPMother->PdgCode() == 111) {
															motherPi0OrEta = kTRUE;	//pi0
														}
													}
													if(fAnalyseMCEta) {
														if(mcPMother->PdgCode() == 221) {
															motherPi0OrEta = kTRUE;	//eta
														}
													}
													 														
													if(motherPi0OrEta) {
														
														if(DecayedToGammaGamma(mcPMother, mcEvent)) {  //analyse only pi0->gammagamma
														
															Int_t mcLabelMother			  = mcPMother->Label();
															Short_t generatorIndexMother = mcPMother->GetGeneratorIndex();
															Int_t pdgCodeMother			  = mcPMother->PdgCode();
															
															//cout<<mcLabelMother<<endl;
															
															//check if mother is added signal
															motherIsHijing = kTRUE;
															if(mcLabelMother > fIHijingMax){
																motherIsHijing = kFALSE;
																
																//  ******** Added Signal generator pi0_1  ******** ********** ******** 
																if(generatorIndexMother == 1 && pdgCodeMother == 111) { //pi0_1
																	motherIsAddedPi01    = kTRUE;  
																	fAddedSignalWeight[0] = CalculateAddedSignalWeightPi0Exp(mcPMother->Pt());
																	fDifferenceClusterParticleE[0].push_back(clusterEminusParticleE);
																	fClusterHasSingleContr[0].push_back(clusHasSingleContrib);
																	
																	//Fill hitmaps
																	if(fFillClusterHitmapsAddedSignals) { // hitmap in local relative coordinates
																		switch(modNrClusterPos)
																		{
																			case 1: 
																				fH2ClusterPosMod1Gen1->Fill(cellX, cellZ);
																				break; 
																			case 3: 
																				fH2ClusterPosMod3Gen1->Fill(cellX, cellZ);
																				break; 
																		}
																	}	
																	
																fH1ClusterEAddedPi0->Fill(clusterE,fAddedSignalWeight[0]);
																
													
																}
																//******** ******** ******** ******** ******** ******** ******** ******** 
																
																if(generatorIndexMother == 2 && pdgCodeMother == 221) { //eta_2
																	motherIsAddedEta1    = kTRUE;  
																	
																	fDifferenceClusterParticleE[2].push_back(clusterEminusParticleE);
																	fClusterHasSingleContr[2].push_back(clusHasSingleContrib);
																}
																if(generatorIndexMother == 3 && pdgCodeMother == 111) { //pi0EMC_3
																	motherIsAddedAtEmcal = kTRUE;  
																}
																
																//  ******** Added Signal generator pi0PHOS_4  ******** ************ ********
																if(generatorIndexMother == 4 && pdgCodeMother == 111) { //pi0PHS_4
																	motherIsAddedPi0PHOS = kTRUE; 
																	fAddedSignalWeight[1] = CalculateAddedSignalWeightPi0Exp(mcPMother->Pt());
																	fTruthPtAddedPi0PHOS = mcPMother->Pt();
																	fDifferenceClusterParticleE[1].push_back(clusterEminusParticleE);
																	fClusterHasSingleContr[1].push_back(clusHasSingleContrib);
																	
																	//Fill hitmaps
																	if(fFillClusterHitmapsAddedSignals) { // hitmap in local relative coordinates
																		switch(modNrClusterPos)
																		{
																			case 1: 
																				fH2ClusterPosMod1Gen2->Fill(cellX, cellZ);
																				break; 
																			case 3: 
																				fH2ClusterPosMod3Gen2->Fill(cellX, cellZ);
																				break; 
																		}
																	}
																	
																	fH1ClusterEAddedPi0PHOS->Fill(clusterE,fAddedSignalWeight[1]);
																	
																}
																//******** ******** ******** ******** ******** ******** ******** ********
																
																if(generatorIndexMother == 5 && pdgCodeMother == 221) { //etaEMC_5
																	motherIsAddedAtEmcal = kTRUE;  																	
																}
																if(generatorIndexMother == 6 && pdgCodeMother == 221) { //etaPHS_6
																	motherIsAddedEtaPHOS = kTRUE;  
																	
																	fDifferenceClusterParticleE[3].push_back(clusterEminusParticleE);
																	fClusterHasSingleContr[3].push_back(clusHasSingleContrib);
																}
																
																//if(mcLabelMother <= fIPi0Max){    
																	//motherIsAddedPi01 = kTRUE;
																	//////fAddedSignalWeight[0] = CalculateAddedSignalWeightPi0(mcPMother->Pt());
																	////Printf("Found cluster from added pi0.  ID:  %d   pdg: %d",mcLabelMother,mcPMother->PdgCode()); 
																//} else if(mcLabelMother > fIPi0Max && mcLabelMother <= fIEtaMax){
																	//motherIsAddedEta1 = kTRUE;
																	////fAddedSignalWeight[2] = CalculateAddedSignalWeightEta(mcPMother->Pt());
																	////Printf("Found cluster from added eta.  ID: %d   pdg: %d",mcLabelMother,mcPMother->PdgCode()); 
																//} else if(mcLabelMother > fIEtaMax && mcLabelMother <= fIPi0EMCMax){
																	//motherIsAddedAtEmcal = kTRUE;
																//} else if(mcLabelMother > fIPi0EMCMax && mcLabelMother <= fIPi0PHOSMax){
																	////Printf("Found cluster from added PHOSpi0.  ID:  %d   pdg: %d",mcLabelMother,mcPMother->PdgCode());
																	//motherIsAddedPi0PHOS = kTRUE;
																	////if(fAddedSignalWeight[1] == CalculateAddedSignalWeightPi0(mcPMother->Pt())) cout<<"second from thad generator"<<endl;
																	////fAddedSignalWeight[1] = CalculateAddedSignalWeightPi0(mcPMother->Pt());
																//} else if(mcLabelMother > fIPi0PHOSMax && mcLabelMother <= fIEtaEMCMax){
																	//motherIsAddedAtEmcal = kTRUE;
																//} else if(mcLabelMother > fIEtaEMCMax && mcLabelMother <= fIEtaPHOSMax){
																	//motherIsAddedEtaPHOS = kTRUE;
																	////cout<<"####"<<endl;
																	////cout<<"####"<<endl;
																	////cout<<"####"<<endl;
																	////cout<<"generator of mcPCluster "<<mcPCluster->GetGeneratorIndex()<<endl;
																	////mcPCluster->Print();
																	////cout<<"generator of mcPMother "<<mcPMother->GetGeneratorIndex()<<endl;
																	////mcPMother->Print();
																	//fAddedSignalWeight[3] = CalculateAddedSignalWeightEta(mcPMother->Pt());
																	////Printf("Found cluster from added PHOSeta.  ID: %d   pdg: %d",mcLabelMother,mcPMother->PdgCode()); 
																//}
															}
															//cout<<"WEIGHTS: "<<fAddedSignalWeight[0]<<" "<<fAddedSignalWeight[1]<<" "<<fAddedSignalWeight[2]<<" "<<fAddedSignalWeight[3]<<" "<<endl;
															
															if(motherIsAddedPi01 || motherIsAddedPi0PHOS || motherIsAddedEta1 || motherIsAddedEtaPHOS || motherIsAddedAtEmcal) {
																motherIsFromAddedSignal = kTRUE;
																//if(mcPCluster->PdgCode() != 22 ) cout<<"mother is pi0 or eta, but cluster isnt from photon!"<<endl;
															}
															
															//for the pi0 shot at phos directly:
															//find the types of decays of the mothers of the particles that made a cluster in PHOS.
															//if(motherIsAddedPi01) {
																//FillHistoWithDaughterInfoPi0(mcPMother, fH1Pi0DecayModesAddedPi0PHOS, mcEvent);
															//} //if(motherIsAddedPi0PHOS)
															
														} //Decayed into 2 photons
													}  //if(motherPi0OrEta) 
												}  //else {if(0) cout<<"imother <= 0 ?!"<<endl;}//no mother 
											}//Cluster is from photon
										} //clusterlabel != -1
									} //if(faddedSignals)
									// ************* END Added Signal Analysis - Determine if cluster is from Added Signal *****************


									// ************* START Fill Cluster Hitmaps **************** //

									fH2HitmapEtaVsPhi->Fill(vpos.Eta(), vpos.Phi()); // eta-phi hitmap, global coordinates
									
									if(fFillClusterHitmaps) { // hitmap in local relative coordinates
										switch(modNrClusterPos)
										{
											case 1: 
												fH2ClusterPositionsMod1->Fill(cellX, cellZ);
												break; 
											case 2: 
												fH2ClusterPositionsMod2->Fill(cellX, cellZ);
												break; 
											case 3: 
												fH2ClusterPositionsMod3->Fill(cellX, cellZ);
												break; 
										}
									}	

									// ************* END Fill Cluster Hitmaps **************** //



									// ************ START Fill Energy Histograms/nTupel *********** //

									fH1ClusterEAfterCuts->Fill(clusterE);
									fH2EAfterCutsVsModNum->Fill(clusterE,modNrClusterPos);

									
									if(fFillNTupelClusterE){
										if     (modNrClusterPos == 1) {fNTupelClusterEnergyMod1->Fill(clusterE);}
										else if(modNrClusterPos == 2) {fNTupelClusterEnergyMod2->Fill(clusterE);}
										else if(modNrClusterPos == 3) {fNTupelClusterEnergyMod3->Fill(clusterE);}           			
									}
									
									// ************ END Fill Energy Histograms/nTupel *********** //

									
									
									// ************ START Fill Photons Vector for Reconstruction *********** //
									
									Photon.SetE(clusterE);
									
									// calculate momentum of the cluster from vertex and cluster position 
									// and cluster energy and save it in TLorentzVector of the photon
									// also saves cluster energy as energy of the TLorentzVector
									virCluster->GetMomentum(Photon,vertex); 
									
									//Skip Added Signals (motherIsFromAddedSignal can only be set to true if fAnalyseAddedSignals is true (meaning added signals are analysed))
									if(!motherIsFromAddedSignal) {
										//Add reconstructed cluster ('photon') to array of vectors. Used later to reconstruct Mass Vs Pt distribution
										fPhotons[0][izvtx][imult].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE) );
										fModNumber[0][izvtx][imult].push_back(modNrClusterPos);
									}
									
									if(fFillHMassPtTiming){
										// Check Cluster Timing 
										clusterTOF   = virCluster->GetTOF();
										if(clusterTOF < fTimingCutMin || clusterTOF > fTimingCutMax) {   
											//If the cluster has bad timing variable is set to TRUE
											badTiming = true;
										}   
										fClusterTiming[0][izvtx][imult].push_back(badTiming);
									}
																		
									
									if(fAnalyseAddedSignals) {
										if(motherIsAddedPi01)
											fPhotonsAdded[0].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE));
										if(motherIsAddedPi0PHOS)
											fPhotonsAdded[1].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE));
										if(motherIsAddedEta1)
											fPhotonsAdded[2].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE));
										if(motherIsAddedEtaPHOS)
											fPhotonsAdded[3].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE));
									}
									
									
									// ************ END Fill Photons Array for Reconstruction *********** //

									// ************ END Good PHOS Clusters ************ //

								} //(badTiming)
							} //fDistToBadCell
						} //(clusterPosBad)
					} //(M02)
				} //(E>MinE)
			} //(nCells)
		} //if(isPHOS)
	}//loop over nclusters

	// ***************** END Cluster loop for data (ESD/AOD) *************** //



	fH1NClustersPHOS->Fill(nclustersPHOS);
	fH1NClustersPHOSafterCuts->Fill(nclustersPHOSafterCuts);


	
	// ***************** START Same Event *************** //

	for(UInt_t i=0; i<fPhotons[0][izvtx][imult].size(); i++){
	    for(UInt_t j=i+1; j<fPhotons[0][izvtx][imult].size(); j++){
			Parent = fPhotons[0][izvtx][imult][i] + fPhotons[0][izvtx][imult][j];
			Double_t deltaphi = GetDeltaPhi(fPhotons[0][izvtx][imult][i],fPhotons[0][izvtx][imult][j]);
			Double_t deltaeta = GetDeltaEta(fPhotons[0][izvtx][imult][i],fPhotons[0][izvtx][imult][j]);
			Double_t pairasym = fabs(fPhotons[0][izvtx][imult][i].Pt()-fPhotons[0][izvtx][imult][j].Pt())/
				 (fPhotons[0][izvtx][imult][i].Pt()+fPhotons[0][izvtx][imult][j].Pt());
			
			fH1Mass        ->Fill(Parent.M());
			fH2DphiDeta->Fill(deltaphi,deltaeta);   
			
			if(fFillHMassPtModules)
			{
				if(fModNumber[0][izvtx][imult][i] == 1 && fModNumber[0][izvtx][imult][j] == 1) 	{fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),1);}
			  else	if(fModNumber[0][izvtx][imult][i] == 1 && fModNumber[0][izvtx][imult][j] == 2) {fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),2);}
			  else 	if(fModNumber[0][izvtx][imult][i] == 1 && fModNumber[0][izvtx][imult][j] == 3) {fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),3);}
			  else 	if(fModNumber[0][izvtx][imult][i] == 2 && fModNumber[0][izvtx][imult][j] == 2) {fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),4);}
			  else 	if(fModNumber[0][izvtx][imult][i] == 2 && fModNumber[0][izvtx][imult][j] == 3) {fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),5);}
			  else 	if(fModNumber[0][izvtx][imult][i] == 3 && fModNumber[0][izvtx][imult][j] == 3) {fH3MPtModules ->Fill(Parent.M(),Parent.Pt(),6);}
			}
			else if(fFillHMassPtTiming)
			{
				if     (fClusterTiming[0][izvtx][imult][i] == false && fClusterTiming[0][izvtx][imult][j] == false) {fH3MPtTiming->Fill(Parent.M(),Parent.Pt(),1);} // both Clusters are within the timing-cut
				else if(fClusterTiming[0][izvtx][imult][i] == false || fClusterTiming[0][izvtx][imult][j] == false) {fH3MPtTiming->Fill(Parent.M(),Parent.Pt(),2);} // one Cluster is within the timing-cut
				else if(fClusterTiming[0][izvtx][imult][i] == true  && fClusterTiming[0][izvtx][imult][j] == true)  {fH3MPtTiming->Fill(Parent.M(),Parent.Pt(),3);} // both Clusters are outside of the timing-cut 
			}
			else  if(!fFillNewAsymmClasses)		//default
			{
			  Int_t asymCut = 0;
			  if     (pairasym<0.1)  asymCut = 1;
			  else if(pairasym<0.7)  asymCut = 2;
			  else                   asymCut = 3;
			  fH3MPtAsymm ->Fill(Parent.M(),Parent.Pt(),asymCut);
			}
			else {fH3MPtAsymm ->Fill(Parent.M(),Parent.Pt(),pairasym);}
		}
	}

	// ***************** END Same Event *************** //



	// ***************** START Mixed Event *************** //

	for(UInt_t i=0; i<fPhotons[0][izvtx][imult].size(); i++){
		for(UInt_t ipool=1; ipool<fgkPoolDepth; ipool++){
			for(UInt_t j=0; j<fPhotons[ipool][izvtx][imult].size(); j++){
				Parent = fPhotons[0][izvtx][imult][i]+fPhotons[ipool][izvtx][imult][j];
				Double_t deltaphi = GetDeltaPhi(fPhotons[0][izvtx][imult][i],fPhotons[ipool][izvtx][imult][j]);
				Double_t deltaeta = GetDeltaEta(fPhotons[0][izvtx][imult][i],fPhotons[ipool][izvtx][imult][j]);
				Double_t pairasym = fabs(fPhotons[0][izvtx][imult][i].Pt()-fPhotons[ipool][izvtx][imult][j].Pt())/
					(fPhotons[0][izvtx][imult][i].Pt()+fPhotons[ipool][izvtx][imult][j].Pt());
				fH1MassMixed        ->Fill(Parent.M());	
				fH2DphiDetaMix->Fill(deltaphi,deltaeta);

				if(fFillHMassPtModules)
				{
					if     (fModNumber[0][izvtx][imult][i] == 1 && fModNumber[ipool][izvtx][imult][j] == 1) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),1);}
					else if(fModNumber[0][izvtx][imult][i] == 1 && fModNumber[ipool][izvtx][imult][j] == 2) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),2);}
					else if(fModNumber[0][izvtx][imult][i] == 1 && fModNumber[ipool][izvtx][imult][j] == 3) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),3);}
					else if(fModNumber[0][izvtx][imult][i] == 2 && fModNumber[ipool][izvtx][imult][j] == 2) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),4);}
					else if(fModNumber[0][izvtx][imult][i] == 2 && fModNumber[ipool][izvtx][imult][j] == 3) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),5);}
					else if(fModNumber[0][izvtx][imult][i] == 3 && fModNumber[ipool][izvtx][imult][j] == 3) {fH3MPtModulesMix ->Fill(Parent.M(),Parent.Pt(),6);}
				}
				else if(fFillHMassPtTiming)
				{ 
					if     (fClusterTiming[0][izvtx][imult][i] == false && fClusterTiming[ipool][izvtx][imult][j] == false) {fH3MPtTimingMix->Fill(Parent.M(),Parent.Pt(),1);}
					else if(fClusterTiming[0][izvtx][imult][i] == false || fClusterTiming[ipool][izvtx][imult][j] == false) {fH3MPtTimingMix->Fill(Parent.M(),Parent.Pt(),2);}
					else if(fClusterTiming[0][izvtx][imult][i] == true  && fClusterTiming[ipool][izvtx][imult][j] == true)  {fH3MPtTimingMix->Fill(Parent.M(),Parent.Pt(),3);}
				}
				else if(!fFillNewAsymmClasses) 
				{
					Int_t asymCut = 0;
					if     (pairasym<0.1)  asymCut = 1;
					else if(pairasym<0.7)  asymCut = 2;
					else                   asymCut = 3;
					fH3MPtAsymmMix ->Fill(Parent.M(),Parent.Pt(),asymCut);
				}
				else {fH3MPtAsymmMix->Fill(Parent.M(),Parent.Pt(),pairasym);}
			}
		} 
	}

	// ***************** END Mixed Event *************** //
	
	// ****************** START Combining Cluster from Added Sginals ************** //
	if(fAnalyseAddedSignals) {
		//cout<<" "<<endl;
		//cout<<" ##########"<<endl;
		//cout<<"SIZES ARE: "<<fPhotonsAdded[0].size()<<" "<<fPhotonsAdded[1].size()<<" "<<fPhotonsAdded[2].size()<<" "<<fPhotonsAdded[3].size()<<endl;
		//cout<<"SIZES ARE: "<<fDifferenceClusterParticleE[0].size()<<" "<<fDifferenceClusterParticleE[1].size()<<" "<<fDifferenceClusterParticleE[2].size()<<" "<<fDifferenceClusterParticleE[3].size()<<endl;
		//cout<<"SIZES ARE: "<<fClusterHasSingleContr[0].size()<<" "<<fClusterHasSingleContr[1].size()<<" "<<fClusterHasSingleContr[2].size()<<" "<<fClusterHasSingleContr[3].size()<<endl;
		//cout<<"VALUES ARE: fDifferenceClusterParticleE[i][0]  "<<fDifferenceClusterParticleE[0][0]<<" "<<fDifferenceClusterParticleE[1][0]<<" "<<fDifferenceClusterParticleE[2][0]<<" "<<fDifferenceClusterParticleE[3][0]<<endl;
		//cout<<"VALUES ARE: fDifferenceClusterParticleE[i][1]  "<<fDifferenceClusterParticleE[0][1]<<" "<<fDifferenceClusterParticleE[1][1]<<" "<<fDifferenceClusterParticleE[2][1]<<" "<<fDifferenceClusterParticleE[3][1]<<endl;
		//cout<<"VALUES ARE: fClusterHasSingleContr[i][0]  "<<fClusterHasSingleContr[0][0]<<" "<<fClusterHasSingleContr[1][0]<<" "<<fClusterHasSingleContr[2][0]<<" "<<fClusterHasSingleContr[3][0]<<endl;
		//cout<<"VALUES ARE: fClusterHasSingleContr[i][1]  "<<fClusterHasSingleContr[0][1]<<" "<<fClusterHasSingleContr[1][1]<<" "<<fClusterHasSingleContr[2][1]<<" "<<fClusterHasSingleContr[3][1]<<endl;
		
		//if(fPhotonsAdded[1].size() > 0) {
			//cout<<"VALUES ARE: fDifferenceClusterParticleE[1][0]  "<<fDifferenceClusterParticleE[1][0]<<endl;
			//cout<<"VALUES ARE: fClusterHasSingleContr[1][0]  "<<fClusterHasSingleContr[1][0]<<endl;
			//cout<<"VALUES ARE: fAddedSignalWeight[1]  "<<fAddedSignalWeight[1]<<endl;
			//if(fPhotonsAdded[1].size() > 1) {
				//cout<<"VALUES ARE: fDifferenceClusterParticleE[1][1]  "<<fDifferenceClusterParticleE[1][1]<<endl;
				//cout<<"VALUES ARE: fClusterHasSingleContr[1][1]  "<<fClusterHasSingleContr[1][1]<<endl;
			//}
		//}
			
		// If two photons from athe same generator hit PHOS, we can combine them and reconstruct the added signal
		// This is only that easy because there is only one added particle per generator. 
		// To implement the Hijing particles (or for a dataset with more than one particle per generator), one would also need to
		// store the mclabel of the mother for every cluster, to only pair particles from the same mother
		
		// ******************* RECONSTRUCT ADDED Pi0s *************************************
		//Fill MassPT histos for reconstructed signals from pi0_1
		if(fPhotonsAdded[0].size() > 1) {
			//if(fPhotonsAdded[0].size() < 3) {
				//Two clusters from Added Pi0	
			Parent = fPhotonsAdded[0][0]+fPhotonsAdded[0][1];
			//cout<<"Reconstrcted Pi0 with (M,P) = ("<<Parent.M()<<","<<Parent.P()<<")"<<endl;
			fH2MPtAddedPi0->Fill(Parent.M(),Parent.Pt(),fAddedSignalWeight[0]);
			fH2MPtAddedPi0_unweighed->Fill(Parent.M(),Parent.Pt());
			//} else AliError("fPhotonsAdded[0] vector to big, should only be size 2!");			
		}
		
		//Fill MassPT histos for reconstructed signals from pi0PHS_4
		if(fPhotonsAdded[1].size() > 1) {
			//Two clusters from Added Pi0-PHOS
			Parent = fPhotonsAdded[1][0]+fPhotonsAdded[1][1];
			//cout<<"Reconstrcted Pi0PHOS with (M,P) = ("<<Parent.M()<<","<<Parent.P()<<")"<<endl;
			fH2MPtAddedPi0PHOS->Fill(Parent.M(),Parent.Pt(),fAddedSignalWeight[1]);	
			fH2MPtAddedPi0PHOS_unweighed->Fill(Parent.M(),Parent.Pt());	
			fH2PtRecVsPtTruthAddedPi0PHOS->Fill(Parent.Pt(),fTruthPtAddedPi0PHOS);
			//cout<<Parent.Pt()<<" "<<fTruthPtAddedPi0PHOS<<endl;
			
			if(fFillMPtForSingleOrMultContrClus) {
				if(fClusterHasSingleContr[1][0] && fClusterHasSingleContr[1][1]) {    //both clusters have only 1 contributer
					//Fill MassPT for SingleContrib
					fH2MPtAddedPi0PHOSSingleContr->Fill(Parent.M(),Parent.Pt(),fAddedSignalWeight[1]);
				}
				if(!fClusterHasSingleContr[1][0] && !fClusterHasSingleContr[1][1]) {  //both clusters dont have more than 1 contributer
					//Fill MassPT for MultiContrib
					fH2MPtAddedPi0PHOSMultContr->Fill(Parent.M(),Parent.Pt(),fAddedSignalWeight[1]);
				}
			}
			
			Float_t maxDifference = 0.0;
			if(abs(fDifferenceClusterParticleE[1][0]) > abs(fDifferenceClusterParticleE[1][1])) {
				maxDifference = fDifferenceClusterParticleE[1][0];
			} else maxDifference = fDifferenceClusterParticleE[1][1];

			//Fill massVsMaxDifference
			fH2MEnergyDiffAddedPi0PHOS->Fill(Parent.M(),maxDifference);
		}
		//**********************************************************************************
		
		
		//******************* RECONSTRUCT ADDED Etas *************************************
		//Fill MassPT histos for reconstructed signals from eta_2
		if(fPhotonsAdded[2].size() > 1) {
			//Two clusters from Added Eta
			Parent = fPhotonsAdded[2][0]+fPhotonsAdded[2][1];
			//cout<<"Reconstrcted Eta with (M,P) = ("<<Parent.M()<<","<<Parent.P()<<")"<<endl;
			fH2MPtAddedEta->Fill(Parent.M(),Parent.Pt());	
		}
		
		//Fill MassPT histos for reconstructed signals from etaPHS_6
		if(fPhotonsAdded[3].size() > 1) {
			//Two clusters from Added Eta-PHOS
			Parent = fPhotonsAdded[3][0]+fPhotonsAdded[3][1];
			//cout<<"Reconstrcted EtaPHOS with (M,P) = ("<<Parent.M()<<","<<Parent.P()<<")"<<endl;
			fH2MPtAddedEtaPHOS->Fill(Parent.M(),Parent.Pt());	
		}
		//**********************************************************************************
		
		
		
		// CLEAR all added Signal Photons (no event mixing needed)
		for(Int_t i = 0; i<4; i++) {
			fPhotonsAdded[i].clear();
			fAddedSignalWeight[i] = 1.0;
			fDifferenceClusterParticleE[i].clear();  
			fClusterHasSingleContr[i].clear();
		}
		fTruthPtAddedPi0PHOS = 0.0;
	}
	// **************************************************************************** //


	// ***************** START Clearing Mixing Variables *************** //
	
	for(Int_t ipool=fgkPoolDepth-1; ipool>0; ipool--){
		fPhotons[ipool][izvtx][imult].clear();
		fModNumber[ipool][izvtx][imult].clear();
		fClusterTiming[ipool][izvtx][imult].clear();
		for(UInt_t i=0; i<fPhotons[ipool-1][izvtx][imult].size(); i++){
			fPhotons[ipool][izvtx][imult].push_back(fPhotons[ipool-1][izvtx][imult][i]);     
			fModNumber[ipool][izvtx][imult].push_back(fModNumber[ipool-1][izvtx][imult][i]);
			fClusterTiming[ipool][izvtx][imult].push_back(fModNumber[ipool-1][izvtx][imult][i]);
		}
	}
	
	fPhotons[0][izvtx][imult].clear();
	fModNumber[0][izvtx][imult].clear();
	fClusterTiming[0][izvtx][imult].clear();

	// ***************** END Clearing Mixing Variables *************** //

	// ***************** END Main Loop Called for each Event *************** //

	// NEW HISTO should be filled before this point, as PostData puts the
	// information for this iteration of the UserExec in the container
	PostData(1, fOutput);
	fEventCounter++;
} //END of UserExec() 


//________________________________________________________________________
void AliAnalysisTaskPHOSNeutralMeson::Terminate(Option_t *) {//specify what you want to have done
	// Called once at the end of the query.
}

//________________________________________________________________________
Int_t AliAnalysisTaskPHOSNeutralMeson::GetZvtxBin(Double_t vertZ)
{
	Int_t izvtx = -1;

	if     (vertZ<-35)
		izvtx=0;
	else if(vertZ<-30)
		izvtx=1;
	else if(vertZ<-25)
		izvtx=2;
	else if(vertZ<-20)
		izvtx=3;
	else if(vertZ<-15)
		izvtx=4;
	else if(vertZ<-10)
		izvtx=5;
	else if(vertZ< -5)
		izvtx=6;
	else if(vertZ<  0)
		izvtx=7;
	else if(vertZ<  5)
		izvtx=8;
	else if(vertZ< 10)
		izvtx=9;
	else if(vertZ< 15)
		izvtx=10;
	else if(vertZ< 20)
		izvtx=11;
	else if(vertZ< 25)
		izvtx=12;
	else if(vertZ< 30)
		izvtx=13;
	else if(vertZ< 35)
		izvtx=14;
	else
		izvtx=15;

	return izvtx;  
}

//________________________________________________________________________
Int_t AliAnalysisTaskPHOSNeutralMeson::GetMultBin(Int_t mult){

	// translates cluster-multiplicity into a multiplicity bin for event mixing
	// set in AddTask:	fMaxMultForMultBinning (default is 60)
	//					fNMultBins			   (default is 30)  !!! maximum value is 30, otherwise array is out of range !!!
	
	Float_t binningFactor = (Float_t) fMaxMultForMultBinning / (Float_t) fNMultBins; 
	Int_t imult = -1;	

	if     (mult<2)
		imult=0;
	else if(mult<fMaxMultForMultBinning)
		imult=mult-2;
	else
		imult=fMaxMultForMultBinning-1;
		
	imult = (int) (imult/binningFactor);  //cast to int to cut away decimal numbers to get discrete bins
	return imult;  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSNeutralMeson::IsGoodChannel(const Char_t * det, Int_t mod, Int_t ix, Int_t iz)
{
	//Check if the Cell at the Position of the Cluster belongs to the good ones
	if(strcmp(det,"PHOS")==0){
		if(mod>5 || mod<1){
			AliError(Form("No bad map for PHOS module %d ",mod)) ;
			return kTRUE ;
		}
		if(!fPHOSBadMap[mod]){
			AliError(Form("No Bad map for PHOS module %d",mod)) ;
			return kTRUE ;
		}
		if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
			return kFALSE ;
		else
			return kTRUE ;
	}
	else{
		AliError(Form("Can not find bad channels for detector %s ",det)) ;
	}
	return kTRUE ;
}

//_________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPHOSNeutralMeson::IsGoodChannelDistanceOnCellLevel(const Char_t * det, Int_t mod, Int_t ix, Int_t iz, Int_t dist)
{
	//Check if the Cells at the Position of the Cluster and in an dist*dist area belong to the good ones
	Bool_t goodCellSurroundings = kTRUE; 

	if(strcmp(det,"PHOS")==0){
		if(mod>5 || mod<1){
			AliError(Form("No bad map for PHOS module %d ",mod)) ;
			return kTRUE ;
		}
		if(!fPHOSBadMap[mod]){
			AliError(Form("No Bad map for PHOS module %d",mod)) ;
			return kTRUE ;
		}
		
		for(Int_t ix0=ix-dist; ix0<=ix+dist; ix0++)
		{
			if(ix0<1 || ix0>=65){continue;}
			else {
		    	for(Int_t iz0=iz-dist; iz0<=iz+dist; iz0++) {
					if(iz0<1 || iz0>=57){continue;}
					else {
						if(fPHOSBadMap[mod]->GetBinContent(ix0,iz0)>0) {goodCellSurroundings = kFALSE;}
					}
				}
			}
		}
	}	
	else{
		AliError(Form("Can not find bad channels for detector %s ",det)) ;
	}
	
	return goodCellSurroundings;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::GetDeltaPhi(TLorentzVector p1, TLorentzVector p2){
	Double_t dphi = p1.Phi() - p2.Phi();

	if(dphi<0.5*TMath::Pi())  
		dphi = dphi + 2.0*TMath::Pi();

	if(dphi>1.5*TMath::Pi())  
		dphi = dphi - 2.0*TMath::Pi();

	return dphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::GetDeltaEta(TLorentzVector p1, TLorentzVector p2){

	Double_t deta = p1.PseudoRapidity() - p2.PseudoRapidity();
	return deta;
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::RecalibratePHOSClusterEnergy(TString calibOption, Double_t E, Int_t run)
{
	//naming convention: first: dataset to be calibrated, second: what it should be closer to (doesn't have to be another dataset) 

	Double_t E_recalib; 

	if(calibOption == "lhc12d_12h") { // move the peaks of 12d closer to 12h
		if(run >= 184964 && run <= 185031) {
			E_recalib = E*1.009155;
		}
		if(run >= 184135 && run <= 184208) {
			E_recalib = E*1.00869; 
		}
		if(run == 185284) {
			E_recalib = E*1.00869; 
		}
	}
	
	if( calibOption == "lhc12d_14e2c") { // moves the peaks of the MC-data closer to lhc12d
		E_recalib = ((1+0.02/(1+TMath::Power(2*E/1.5,2)))*0.985)*E; 
	}

	return E_recalib;
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::CalculateAddedSignalWeightPi0Tsalis(Double_t pt)
{
	//Calculates the value of a tsalis fit for a given pT.
	//Can be used to calculate filling weight for reconstructed (m_inv, pT). 
	//Use truth pT of simulated particle as argument!
	Double_t weight = 1.0;
	//TF1 *fTsalis = new TF1("fTsalis","2.0*3.14159*x*[0]/pow(exp(-[1]*x-[2]*x*x) +x/[3],[4])",0.0,30);
	//fTsalis->SetParameters(fTsalisPi0Param1,fTsalisPi0Param2,fTsalisPi0Param3,fTsalisPi0Param4,fTsalisPi0Param5);
	//weight = fTsalis->Eval(pT);
	//delete fTsalis;
	return weight;
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::CalculateAddedSignalWeightPi0Exp(Double_t pt)
{
	//Calculates the value of an exponential fit for a given pT.
	//Can be used to calculate filling weight for reconstructed (m_inv, pT). 
	//Use truth pT of simulated particle as argument!
	Double_t weight = 1.0;
	Double_t par0 = fExpParam1;
	Double_t par1 = fExpParam2;
	Double_t par2 = fExpParam3;
	Double_t par3 = fExpParam4;
	
	weight = (par0*pow(par1/((par1*exp(-par3*pt)+pt)),par2));
	
	return weight;
}


//_________________________________________________________________________________
Double_t AliAnalysisTaskPHOSNeutralMeson::CalculateAddedSignalWeightEta(Double_t pt)
{
	
	// Needs to be implemented. So far only a copy of pi0 function, but partly commented out.
	// Simply returns 1.0 so far.
	// Will probably work the same way as for pi0 but needs its own set of parameters as member variables
	
	//Calculates the value of a tsalis fit for a given pT.
	//Can be used to calculate filling weight for reconstructed (m_inv, pT). 
	//Use truth pT of simulated particle as argument!
	Double_t weight = 1.0;
	//TF1 *fTsalis = new TF1("fTsalis","2.0*3.14159*x*[0]/pow(exp(-[1]*x-[2]*x*x) +x/[3],[4])",0.0,30);
	//fTsalis->SetParameters(fTsalisPi0Param1,fTsalisPi0Param2,fTsalisPi0Param3,fTsalisPi0Param4,fTsalisPi0Param5);
	//weight = fTsalis->Eval(pT);
	//delete fTsalis;
	return weight;
}

Int_t AliAnalysisTaskPHOSNeutralMeson::GetNDaughtersOfMCParticle(AliMCParticle* mcP) {
	
	//Returns the number of daughters for an AliMCParticle
	//(This function only exists for AliAODMCParticles - implementation was taken from that class)
	// virtual Int_t GetNDaughters  () const { return fDaughter[1]>0 ? fDaughter[1]-fDaughter[0]+1 : (fDaughter[0]>0 ? 1:0 ) ;}
	
	//MC Particles save the mc-label of their their first and their last daughter
	//all daughters are on the stack between those labels (including of course)
	//The default value for first and last daughter are different for AliMCParticle (-1) and AliAODMCParticle (0)
	
	Int_t daughter[2] = {-1,-1};
	daughter[0] = mcP->GetDaughterFirst();
	daughter[1] = mcP->GetDaughterLast();
	return daughter[1]>0 ? daughter[1]-daughter[0]+1 : (daughter[0]>0 ? 1:0 );
	
}

Bool_t AliAnalysisTaskPHOSNeutralMeson::DecayedToGammaGamma(AliMCParticle * mcP, AliMCEvent* mcEvent) {
	//
	// Returns bool whether particle has exaclty 2 daughters and these daughters are both photons

	Int_t daughter[2] = {-1,-1};
	daughter[0] = mcP->GetDaughterFirst();
	daughter[1] = mcP->GetDaughterLast();
	Int_t nDaughters = daughter[1]>0 ? daughter[1]-daughter[0]+1 : (daughter[0]>0 ? 1:0 );
		  
	if(nDaughters == 2) {
		AliMCParticle *daughter1 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]));
		AliMCParticle *daughter2 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[1]));
		if(daughter1->PdgCode() == 22 && daughter2->PdgCode() == 22) { //gamma gamma decay
			return kTRUE;
		} else return kFALSE;
	} else return kFALSE;

}

Bool_t AliAnalysisTaskPHOSNeutralMeson::MCParticleIsPhysicalDecay(AliMCParticle * mcP, AliMCEvent* mcEvent) {

	//It is ok if the pi0 has 0 or 1 daughters. 
	//Because not all daughters are stored in ALICE stack, only those which interacted with any ALICE active volume in given simulation 
	//=> pi0s without daughters or only 1 daughter are technically possible.
	
	//Not ok are pi0s that have  (has so far only been observed in lhc13b2_efix (DPMJET)):
	//  2 daughters but not gamma gamma
	//  3 daughters but not dalitz
	//  more than 3 daughters
	
	
	Int_t daughter[2] = {-1,-1};
	daughter[0] = mcP->GetDaughterFirst();
	daughter[1] = mcP->GetDaughterLast();
	Int_t nDaughters = daughter[1]>0 ? daughter[1]-daughter[0]+1 : (daughter[0]>0 ? 1:0 );
	
	// default return value is true
	// only set to false if any of the conditions explained above are met
	Bool_t isOk = kTRUE;
	
	
	if(nDaughters == 2) { // gamma gamma
		AliMCParticle *daughter1 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]));
		AliMCParticle *daughter2 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[1]));
		if(daughter1->PdgCode() == 22 && daughter2->PdgCode() == 22) { //gamma gamma decay
		} else {
			isOk = kFALSE;	//2 daughters but not gamma gamma
		}
	} else if(nDaughters == 3) { //dalitz ??
		AliMCParticle *daughter1 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]));
		AliMCParticle *daughter2 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]+1));
		AliMCParticle *daughter3 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[1]));
		Int_t pdgCode1 = daughter1->PdgCode();
		Int_t pdgCode2 = daughter2->PdgCode();
		Int_t pdgCode3 = daughter3->PdgCode();
		Bool_t oneElectron = kFALSE;
		Bool_t onePositron = kFALSE;
		Bool_t oneGamma = kFALSE;
		if(pdgCode1 == 11  ^ pdgCode2 == 11  ^ pdgCode3 == 11)  oneElectron = kTRUE;  //exactly 1 electron
		if(pdgCode1 == -11 ^ pdgCode2 == -11 ^ pdgCode3 == -11) onePositron = kTRUE;
		if(pdgCode1 == 22  ^ pdgCode2 == 22  ^ pdgCode3 == 22)  oneGamma = kTRUE;
		if(oneElectron && onePositron && oneGamma) { //dalitz decay
		} else {
			isOk = kFALSE;	//3 daughters but not dalitz
		}
	} else if(nDaughters >= 4 ) {
		isOk = kFALSE;  //4 or more daughters
	}
	
	return isOk;

}


void AliAnalysisTaskPHOSNeutralMeson::FillHistoWithDaughterInfoPi0(AliMCParticle * mcP, TH1F* histo, AliMCEvent* mcEvent) {
	//Can be used to fill histos such as fH1Pi0DecayModes or fH1Pi0DecayModesAddedPi0PHOS
	//Currently only works for AODs
	//Histos need to implemented like this:
	//fH1Pi0DecayModes = new TH1F("fH1Pi0DecayModes","fH1Pi0DecayModes",10,0.5,10.5);
   //fH1Pi0DecayModes->GetXaxis()->SetTitle("type of decay / reaction");
   //fH1Pi0DecayModes->GetYaxis()->SetTitle("##pi^{0}");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(1,"0");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(2,"1");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(3,"#gamma #gamma");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(4,"2 (not #gamma #gamma)");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(6,"3 (not e^{+}e^{-}#gamma)");
   //fH1Pi0DecayModes->GetXaxis()->SetBinLabel(7,">=4");
	
	Int_t daughter[2] = {-1,-1};
	daughter[0] = mcP->GetDaughterFirst();
	daughter[1] = mcP->GetDaughterLast();
	Int_t nDaughters = daughter[1]>0 ? daughter[1]-daughter[0]+1 : (daughter[0]>0 ? 1:0 );
	
	if(nDaughters == 0) histo->Fill(1);
	if(nDaughters == 1) histo->Fill(2);   //these are gamma, or e+ or e-  
	if(nDaughters == 2) {
		AliMCParticle *daughter1 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]));
		AliMCParticle *daughter2 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[1]));
		if(daughter1->PdgCode() == 22 && daughter2->PdgCode() == 22) { //gamma gamma decay
			histo->Fill(3);
		} else {
			histo->Fill(4); 
			//cout<<"################################################"<<endl;
			//cout<<"Found a pi0 with 2 daughters, but not gamma gamma: "<<endl;
			//mcP->Print();
			//cout<<"~~~~~ These are the daughters: ~~~~~ "<<endl;
			//daughter1->Print();
			//daughter2->Print();			
		}
	}
	if(nDaughters == 3) { //dalitz ??
		AliMCParticle *daughter1 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]));
		AliMCParticle *daughter2 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[0]+1));
		AliMCParticle *daughter3 = static_cast< AliMCParticle *>(mcEvent->GetTrack(daughter[1]));
		Int_t pdgCode1 = daughter1->PdgCode();
		Int_t pdgCode2 = daughter2->PdgCode();
		Int_t pdgCode3 = daughter3->PdgCode();
		Bool_t oneElectron = kFALSE;
		Bool_t onePositron = kFALSE;
		Bool_t oneGamma = kFALSE;
		if(pdgCode1 == 11  ^ pdgCode2 == 11  ^ pdgCode3 == 11)  oneElectron = kTRUE;  //exactly 1 electron
		if(pdgCode1 == -11 ^ pdgCode2 == -11 ^ pdgCode3 == -11) onePositron = kTRUE;
		if(pdgCode1 == 22  ^ pdgCode2 == 22  ^ pdgCode3 == 22)  oneGamma = kTRUE;
		if(oneElectron && onePositron && oneGamma) { //dalitz decay
			histo->Fill(5);
			//cout<<"################################################"<<endl;
			//cout<<"Found a pi0 with 3 daughters, and it is a dalitz: "<<endl;
			//mcP->Print();
			//cout<<"~~~~~ These are the daughters: ~~~~~ "<<endl;
			//daughter1->Print();
			//daughter2->Print();
			//daughter3->Print();
		} else {
			histo->Fill(6);
			//cout<<"##########################################################"<<endl;
			//cout<<"Found a pi0 with 3 daughters, but its NOT a dalitz decay: "<<endl;
			//mcP->Print();
			//cout<<"~~~~~ These are the daughters: ~~~~~ "<<endl;
			//daughter1->Print();
			//daughter2->Print();
			//daughter3->Print();
		}
	}
	if(nDaughters >= 4) histo->Fill(7);		
}

void AliAnalysisTaskPHOSNeutralMeson::FillDecayGammaHistos(AliMCParticle * mcP, AliMCEvent* mcEvent,  Float_t wgtPi0,
																					TH1F* fH1DecGammAddPi0Eta,
																					TH1F* fH1DecGammAddPi0Y,
																					TH1F* fH1DecGammAddPi0Phi,
																					TH1F* fH1DecGammAddPi0E,
																					TH1F* fH1DecGammInPHOSAddPi0E,
																					TH1F* fH1DecGammAddPi0Asymm,
																					TH1F* fH1DecGammAddPi0OpAngle,
																					TH1F* fH1DecGammAddPi0ConvR,
																					TH1F* fH1DecGammAddPi0ConvRate,
																					TH1F* fH1TruthPtAddedPi0GammasInPHOS) {

	//Call this for the pi0s from added signal generators that met the MesonInPHOS condition
	//This is meant to help compare addedSignalGenerators on the decay-photon level
	
	//fH1DecGammAddPi0Eta
	//fH1DecGammAddPi0Y
	//fH1DecGammAddPi0Phi
	//fH1DecGammAddPi0E
	//fH1DecGammInPHOSAddPi0E
	//fH1DecGammAddPi0Asymm
	//fH1DecGammAddPi0OpAngle
	//fH1DecGammAddPi0ConvR
	//fH1DecGammAddPi0ConvRate
	//fH1TruthPtAddedPi0GammasInPHOS
	
	Int_t daughter[2] = {-1,-1};
	daughter[0] = mcP->GetDaughterFirst();
	daughter[1] = mcP->GetDaughterLast();
	Int_t nDaughters = daughter[1]>0 ? daughter[1]-daughter[0]+1 : (daughter[0]>0 ? 1:0 );
		
	if(nDaughters == 2) {
		
		Bool_t bothGammasInPHOS = kTRUE;
		Double_t eta_d[2] = {0.0,0.0};
		Double_t y_d[2] = {0.0,0.0};
		Double_t phi_d[2] = {0.0,0.0};
		Double_t E_d[2] = {0.0,0.0};
		TLorentzVector gammaVector[2]; 		
		
		for (Int_t dIndex=0; dIndex<nDaughters; dIndex++){	
			
			const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(daughter[dIndex]));
			eta_d[dIndex] = dmc->Eta();
			y_d[dIndex] = dmc->Y();
			phi_d[dIndex] = dmc->Phi();
			E_d[dIndex] = dmc->E();
			gammaVector[dIndex] = TLorentzVector(dmc->Px(),dmc->Py(),dmc->Pz(),dmc->E());
			
			fH1DecGammAddPi0Eta->Fill(eta_d[dIndex]);
			fH1DecGammAddPi0Y	 ->Fill(y_d[dIndex]);
			fH1DecGammAddPi0Phi->Fill(phi_d[dIndex]*(180/TMath::Pi()));
			fH1DecGammAddPi0E	 ->Fill(E_d[dIndex],wgtPi0);
			
			if(eta_d[dIndex] >= fEtaAccMin && eta_d[dIndex] <= fEtaAccMax && phi_d[dIndex] >= fPhiAccMin && phi_d[dIndex] <= fPhiAccMax) {
				fH1DecGammInPHOSAddPi0E->Fill(E_d[dIndex],wgtPi0);
			} else {
				bothGammasInPHOS = kFALSE;
			}
			
			
			Int_t grandDaughter[2] = {-1,-1};
			grandDaughter[0] = dmc->GetDaughterFirst();
			grandDaughter[1] = dmc->GetDaughterLast();
			Int_t nGrandDaughters = grandDaughter[1]>0 ? grandDaughter[1]-grandDaughter[0]+1 : (grandDaughter[0]>0 ? 1:0 );
			
			Bool_t didNotConvertBeforePHOS = kFALSE;
			
			if(nGrandDaughters > 0) {
				//Get the granddaughters
				const AliMCParticle *gdmc1 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetDaughterFirst()));
				const AliMCParticle *gdmc2 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetDaughterLast()));
				Double_t productionR1 = TMath::Sqrt(gdmc1->Xv()*gdmc1->Xv() + gdmc1->Yv()*gdmc1->Yv());
				fH1DecGammAddPi0ConvR->Fill(productionR1);
				
				if((gdmc1->PdgCode()== -1.0*gdmc2->PdgCode()) &&
					(gdmc1->PdgCode()==11 || gdmc1->PdgCode()==-11)) {
					
				}  //else {
				//	cout<<"photon has daughters but not e+e-"<<endl;
				//}
			
				//check if this is gamma->e+e- conversion
				if((gdmc1->PdgCode()== -1.0*gdmc2->PdgCode()) &&
					(gdmc1->PdgCode()==11 || gdmc1->PdgCode()==-11) &&
					productionR1<460.0) { //460 PHOS distance to beamline
				
				fH1DecGammAddPi0ConvRate->Fill(0); //yes
				} else {
					didNotConvertBeforePHOS = kTRUE;
				}
			} //if(nGrandDaughters > 0)
			
			if(didNotConvertBeforePHOS) {
				fH1DecGammAddPi0ConvRate->Fill(1); //no
			}
			
		} //end of loop over gammas
		
		if(bothGammasInPHOS) {
			fH1TruthPtAddedPi0GammasInPHOS->Fill(mcP->Pt(),wgtPi0);
		}
		
		Double_t pairasymm = fabs(E_d[0]-E_d[1])/(E_d[0]+E_d[1]);
		fH1DecGammAddPi0Asymm->Fill(pairasymm);


		Float_t openingAngle = gammaVector[0].Angle(gammaVector[1].Vect());
		
		fH1DecGammAddPi0OpAngle->Fill(openingAngle);

	} //if(nDaughters == 2)

}	//::FillDecayGammaHistos
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
