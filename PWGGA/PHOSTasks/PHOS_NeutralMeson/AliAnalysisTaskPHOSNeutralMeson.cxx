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

// EMCAL includes
#include "AliEMCALRecoUtils.h"

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
	fApplyBadMapManually(false),
	fMaxMultForMultBinning(60),   //default makes sense for pPb. Set differently in your AddTask if needed.
	fNMultBins(30),
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
	fH1MCpionVertDistToEventVert(0),
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
	fH1Pi0TruthPtPhi2piY03(0),
	fH2Pi0TruthPhiEta(0), 
	fH1ElectronConversionR(0),
	fH1Pi0TruthPtPhotonsPhos(0), 
	fH1K0Pi0TruthPtPhotonsPhos(0), 
	fH1PriPi0TruthPtPhotonsPhos(0),
	fH1SecPi0TruthPtPhotonsPhos(0),  
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
	fNTupelClusterEnergyMod1(0),
	fNTupelClusterEnergyMod2(0),
	fNTupelClusterEnergyMod3(0),
	fPHOSGeo(0),
	fPHOSCalibData(0), //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
	fUtils(0)
{
    // Dummy constructor ALWAYS needed for I/O.

    // Initialize the PHOS geometry
    // Set bad channel map
    Char_t key[55] ;
    for(Int_t i=0; i<6; i++){
		snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
		fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
	}

	fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");     
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
	fH1MCpionVertDistToEventVert(0),
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
	fH1Pi0TruthPtPhi2piY03(0),
	fH2Pi0TruthPhiEta(0), 
	fH1ElectronConversionR(0),
	fH1Pi0TruthPtPhotonsPhos(0),
	fH1K0Pi0TruthPtPhotonsPhos(0),
	fH1PriPi0TruthPtPhotonsPhos(0),
	fH1SecPi0TruthPtPhotonsPhos(0), 
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
	fNTupelClusterEnergyMod1(0),
	fNTupelClusterEnergyMod2(0),
	fNTupelClusterEnergyMod3(0),
	fPHOSGeo(0),
	fPHOSCalibData(0), //neccesary for cell by cell calibration, before filling CellID_vs_E histos.
	fUtils(0)
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
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");    
   
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

	cout << "__________AliAnalysisTaskFraNeutralMeson: Input settings__________" << endl;
	cout << " fMcMode:             " << fMcMode   << endl;
	cout << " number of zvtx bins: " << fgkZvtxBins << endl;
	cout << " number of mult bins: " << fNMultBins << endl;
	cout << " poolDepth:           " << fgkPoolDepth << endl;
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

    Int_t Mbins = 3000;
    Float_t Mlow = 0.0, Mup = 3.0;
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

    Int_t ptbins = 200;
    Float_t ptlow = 0.0, ptup = 40.0;
    Int_t Ebins = 1000;
    Float_t Elow = 0.0, Eup = 20.0;
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

    fH1MCpionVertDistToEventVert = new TH1F("fH1MCpionVertDistToEventVert", "fH1MCpionVertDistToEventVert", 50000, 0, 5000);
    fH1MCpionVertDistToEventVert->GetXaxis()->SetTitle("dR");
    fH1MCpionVertDistToEventVert->GetYaxis()->SetTitle("counts");
    fH1MCpionVertDistToEventVert->SetMarkerStyle(kFullCircle);
    TotalNBins+=10000;

    fH1Pi0TruthPt = new TH1F("fH1Pi0TruthPt", "P_{T} distribution for Truth Pi0's", ptbins, ptlow, ptup);
    fH1Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fH1Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fH1Pi0TruthPt->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;

    if(fFillFeedDownHistos)
	{
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
    }

    fH1Pi0TruthPtPhos = new TH1F("fH1Pi0TruthPtPhos", "P_{T} distribution for Truth Pi0's (hit Phos)", ptbins, ptlow, ptup);
    fH1Pi0TruthPtPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fH1Pi0TruthPtPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fH1Pi0TruthPtPhos->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;

    if(fFillFeedDownHistos)
    {
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
    }

    fH1Pi0TruthPtPhi2PiY05 = new TH1F("fH1Pi0TruthPtPhi2PiY05", "P_{T} for Truth Pi0's [|y_{#pi^{0}}|<0.5 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
    fH1Pi0TruthPtPhi2PiY05->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fH1Pi0TruthPtPhi2PiY05->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fH1Pi0TruthPtPhi2PiY05->SetMarkerStyle(kFullCircle);
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

    if(fFillFeedDownHistos)
    {
		fH1ElectronConversionR = new TH1F("fH1ElectronConversionR", "conversion point (radius)", 600,0,600);
		fH1ElectronConversionR->GetXaxis()->SetTitle("production radius [cm]");
		fH1ElectronConversionR->SetMarkerStyle(kFullCircle);
		TotalNBins+=600; 

		fH1Pi0TruthPtPhotonsPhos = new TH1F("fH1Pi0TruthPtPhotonsPhos", "P_{T} distribution for Pions with both photons (in Phos)", ptbins, ptlow, ptup);
    	fH1Pi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    	fH1Pi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    	fH1Pi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
    	TotalNBins+=ptbins;

    	fH1K0Pi0TruthPtPhotonsPhos = new TH1F("fH1K0Pi0TruthPtPhotonsPhos", "P_{T} distribution for Pions from K0 with both photons (in Phos)", ptbins, ptlow, ptup);
    	fH1K0Pi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    	fH1K0Pi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    	fH1K0Pi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
    	TotalNBins+=ptbins;

    	fH1PriPi0TruthPtPhotonsPhos = new TH1F("fH1PriPi0TruthPtPhotonsPhos", "P_{T} distribution for primary Pions with both photons (in Phos)", ptbins, ptlow, ptup);
    	fH1PriPi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    	fH1PriPi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    	fH1PriPi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
    	TotalNBins+=ptbins;

    	fH1SecPi0TruthPtPhotonsPhos = new TH1F("fH1SecPi0TruthPtPhotonsPhos", "P_{T} distribution for secondary Pions with both photons (in Phos)", ptbins, ptlow, ptup);
		fH1SecPi0TruthPtPhotonsPhos->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    	fH1SecPi0TruthPtPhotonsPhos->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    	fH1SecPi0TruthPtPhotonsPhos->SetMarkerStyle(kFullCircle);
		TotalNBins+=ptbins;
   }

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
	fOutput->Add(fH1MCpionVertDistToEventVert);
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
	fOutput->Add(fH1Pi0TruthPtPhi2piY03);

	fOutput->Add(fH2Pi0TruthPhiEta);
	fOutput->Add(fH1ElectronConversionR);
	fOutput->Add(fH1Pi0TruthPtPhotonsPhos);
	fOutput->Add(fH1K0Pi0TruthPtPhotonsPhos);  
	fOutput->Add(fH1PriPi0TruthPtPhotonsPhos);
	fOutput->Add(fH1SecPi0TruthPtPhotonsPhos); 

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



	// ***************** START zVertex Cut ***************** //

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

	// ***************** END zVertex Cut ***************** //



	// ********************* START Setting PHOS matrix   ****************************** //

	Int_t runNumber = 0;
	runNumber = fAnyEv->GetRunNumber();

	if (fPHOSGeo==0) {


	
	}

	if(fEventCounter == 0) { // Only done for the first Event
		AliOADBContainer geomContainer("phosGeo");
		geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
		
		TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
		for(Int_t mod=0; mod<5; mod++) {
			if(!matrixes->At(mod)) continue;
			fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
			printf("....TASK.....Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
			((TGeoHMatrix*)matrixes->At(mod))->Print() ;
		}

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
   


	//******************* START MC part ************************ //

	if(isMC){
		 		 
		Int_t isPrimary     = 0;
		Int_t isK0sDecay    = 0;
		Int_t isMaterialSec    = 0;

		if (!mcEvent){
			AliError("Event is not an MC event"); 
			return;
		}

		const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
		if (!evtVtx) return;


		mcEvent->PreReadAll();    
		Int_t nTracksMC  = mcEvent->GetNumberOfTracks();
		Int_t nPTracksMC = mcEvent->GetNumberOfPrimaries();

		// Loop over all MC-particles 
		for (Int_t iTrack = 0; iTrack<nTracksMC; ++iTrack) {

			AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
			if (!mcP) continue;

			//check if its a primary particle
			if(iTrack<nPTracksMC)  isPrimary = 1;
			else                   isPrimary = 0;

			//check if Particle is a K0 decays
			isK0sDecay = 0;
			if(mcP->GetMother()>-1){
				if( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  310 ||
						((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -310  )
					isK0sDecay = 1;
			}
			
			// it's a pion !! 
			if(mcP->PdgCode() != 111) continue; // for pi0s
			//if(mcP->PdgCode() != 221) continue; // for etas 
			
			// and it's close enough to the event vertex
			Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) +
					(mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));

			fH1MCpionVertDistToEventVert->Fill(dR);

			if(dR > 0.01) continue; //skip particles too far away from event vertex.

			// Feed-Down
			if(fFillFeedDownHistos) {
			
				Int_t daughter[2] = {-1,-1};
				daughter[0] = mcP->GetFirstDaughter();
				daughter[1] = mcP->GetLastDaughter();

				// Make sure that the pion has 2 daughters
				if (daughter[0]<0)  continue;			// if the pion does not decay...
				if (daughter[1]<0)  daughter[1]=daughter[0];   // if the pion only has one daughter ...?   
				if (daughter[1]-daughter[0] != 1)  continue;	// if the Pion has more than 2 daughters (Dalitz-decay)

				Int_t eIndexofConvertedPhoton[2] = {-1,-1};

				bool bacc = true; 	
				bool binp = true;
				Double_t eta_d[2] = {0.0,0.0};
				Double_t phi_d[2] = {0.0,0.0};

				for (Int_t daughter_index=0; daughter_index<2; daughter_index++){
				
					const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(daughter[daughter_index]));
					eta_d[daughter_index] = dmc->Eta();
					phi_d[daughter_index] = dmc->Phi();
			  
					if(!(dmc->PdgCode()==22))	  binp = false;
					if(!(dmc->PdgCode()==22 && 
					eta_d[daughter_index]>fEtaAccMin && eta_d[daughter_index]<fEtaAccMax && 
					phi_d[daughter_index]>fPhiAccMin && phi_d[daughter_index]<fPhiAccMax))   bacc = false;	
			  
					if(dmc->GetFirstDaughter()>0 && dmc->GetLastDaughter()>0) {
						
						// get the photons's daughters... 
						const AliMCParticle *dmcd1 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetFirstDaughter()));
						const AliMCParticle *dmcd2 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetLastDaughter()));
						Double_t productionR1 = TMath::Sqrt(dmcd1->Xv()*dmcd1->Xv() + dmcd1->Yv()*dmcd1->Yv());
						if(bacc)  fH1ElectronConversionR->Fill(productionR1);
						
						// check if this is a conversion... 
						if( (dmcd1->PdgCode()== -1.0*dmcd2->PdgCode()) &&
						(dmcd1->PdgCode()==11 || dmcd1->PdgCode()==-11) &&
						productionR1<460.0){ //460 PHOS distance to beamline 
							
							//find the conv e with highest energy, assign it to be that photon decay product.
							if( dmcd1->E() > dmcd2->E() )
								eIndexofConvertedPhoton[daughter_index] = dmc->GetFirstDaughter();
							else
								eIndexofConvertedPhoton[daughter_index] = dmc->GetLastDaughter();
						}
					}
				}

				if(binp!=true) continue; 

				isMaterialSec = 0;

				if(isPrimary!=1){
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
				if(isK0sDecay)  					fH1K0Pi0TruthPt->Fill(mcP->Pt());
				if(isPrimary)	   					fH1PriPi0TruthPt->Fill(mcP->Pt());
				if(isPrimary!=1 && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPt->Fill(mcP->Pt());
				
				// If you get "floating point exeption" in the testtrain use a newer version of aliroot there might be a bug for pythia
				if(mcP->Y()>-0.5 && mcP->Y()<0.5){
												fH1Pi0TruthPtPhi2PiY05->Fill(mcP->Pt()); 
					if(isK0sDecay)  					fH1K0Pi0TruthPtPhi2PiY05->Fill(mcP->Pt());
					if(isPrimary)   					fH1PriPi0TruthPtPhi2PiY05->Fill(mcP->Pt());
					if(isPrimary!=1 && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPtPhi2PiY05->Fill(mcP->Pt());
				}
				
				if(mcP->Y()>-0.3 && mcP->Y()<0.3) {
					fH1Pi0TruthPtPhi2piY03->Fill(mcP->Pt()); 
				}
			
				if(bacc == 1) {
												fH1Pi0TruthPtPhotonsPhos->Fill(mcP->Pt());
					if(isK0sDecay)  					fH1K0Pi0TruthPtPhotonsPhos->Fill(mcP->Pt());
					if(isPrimary)	   					fH1PriPi0TruthPtPhotonsPhos->Fill(mcP->Pt());
					if(isPrimary!=1 && isMaterialSec!=1 && isK0sDecay !=1)	fH1SecPi0TruthPtPhotonsPhos->Fill(mcP->Pt());
				}
			}  
		}//for(nTracksMC) (end of the go-through-all-tracks loop) 
	}//if(isMC)

	//******************* END MC part ************************ //



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
									


									// ************ START Good PHOS Clusters ************ //

									nclustersPHOSafterCuts++; 



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
									
									fPhotons[0][izvtx][imult].push_back( TLorentzVector(Photon.Px(),Photon.Py(),Photon.Pz(),clusterE) );
									fModNumber[0][izvtx][imult].push_back(modNrClusterPos);
									
									if(fFillHMassPtTiming){
									
										// Check Cluster Timing 
										clusterTOF   = virCluster->GetTOF();
										if(clusterTOF < fTimingCutMin || clusterTOF > fTimingCutMax) {   
											//If the cluster has bad timing variable is set to TRUE
											badTiming = true;
										}   
										fClusterTiming[0][izvtx][imult].push_back(badTiming);
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
