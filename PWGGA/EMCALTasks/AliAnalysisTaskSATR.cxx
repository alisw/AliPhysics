/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// --- Root ---
#include <Riostream.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TList.h>

// --- AliRoot ---
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVCaloTrigger.h"
#include "AliVCluster.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"
#include "AliCaloCalibPedestal.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskSATR.h"

ClassImp(AliAnalysisTaskSATR)

//________________________________________________________________________
AliAnalysisTaskSATR::AliAnalysisTaskSATR() : 
  AliAnalysisTaskSE(),
  fOutput(0),
  fHistEclus(0),
  fHistEmaxClus(0),
  fHistEtavsPhiMaxClus(0),
  fHistEtavsEmaxClus(0),
  fHistPhivsEmaxClus(0),
  fHistTOFvsEclus(0),
  fHistTOFvsEclusC(0),
  fHistNcellsvsEclus(0),
  fHistAmpTClus(0),
  fHistAmpMaxTClus(0),
  fHistEtavsPhiMaxTClus(0),
  fHistEmaxClusvsAmpMaxTClus(0),
  fHistEmaxClusvsAmpMatchedTClus(0),
  fHistEmaxClusNotMatchingTClus(0),
  fHistEtavsPhiMaxClusNotMatchingTClus(0),
  fHistEmatchedClusvsAmpMaxTClus(0),
  fHistAmpMaxTClusNotMatchingClus(0),
  fHistEtavsPhiMaxTClusNotMatchingClus(0),
  fHistIdxMaxClusvsIdxMaxTClus(0),
  fHistPhiMaxClusvsPhiMaxTClus(0),
  fHistEtaMaxClusvsEtaMaxTClus(0),
  fHistTOFmaxClusvsTimeMaxTClus(0),
  fHistEmatchedClusvsAmpMatchedTClus(0),
  fHistEmatchedClus(0),
  fHistEmaxMatchedClus(0),
  fHistAmpL1TimeSum(0),
  fHistAmpMaxL1TimeSum(0),
  fHistAmpMaxL1TimeSumVScent(0),
  fHistAmpFastORvsAmpL1TimeSum(0),
  fHistAmpFastOR(0),
  fHistAmpMaxFastOR(0),
  fHistTimeFastOR(0),
  fHistEtavsPhiFastOR(0),
  fHistEtavsPhiMaxFastOR(0),
  fHistTimeDispFastOR(0),
  fHistTimevsL0TimeFastOR(0),
  fHistNtimesFastOR(0),
  fHistEcells(0),
  fHistEmaxCell(0),
  fHistTOFvsEcells(0),
  fHistTOFvsEcellsC(0),
  fHistEmaxCellvsAmpFastOR(0),
  fCaloClustersName("CaloClusters"),
  fTriggerClustersName("triggerClusters"),
  fMinCutL0Amp(-1),
  fMaxCutL0Amp(999),
  fMinCutClusEnergy(-1),
  fMaxCutClusEnergy(999),
  fTimeCutOn(0),
  fMinL0Time(-20),
  fMaxL0Time(20),
  fMinClusTime(-1),
  fMaxClusTime(1),
  fCheckDeadClusters(0),
  fPedestal(0),
  fLoadPed(0),
  fOCDBpath(""),
  fMinDistanceFromBadTower(0.5),
  fClusterizer(0),
  fTriggerClusterizer(0),
  fRun(-1)
{
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskSATR::AliAnalysisTaskSATR(const char *name) : 
  AliAnalysisTaskSE(name),
  fOutput(0),
  fHistEclus(0),
  fHistEmaxClus(0),
  fHistEtavsPhiMaxClus(0),
  fHistEtavsEmaxClus(0),
  fHistPhivsEmaxClus(0),
  fHistTOFvsEclus(0),
  fHistTOFvsEclusC(0),
  fHistNcellsvsEclus(0),
  fHistAmpTClus(0),
  fHistAmpMaxTClus(0),
  fHistEtavsPhiMaxTClus(0),
  fHistEmaxClusvsAmpMaxTClus(0),
  fHistEmaxClusvsAmpMatchedTClus(0),
  fHistEmaxClusNotMatchingTClus(0),
  fHistEtavsPhiMaxClusNotMatchingTClus(0),
  fHistEmatchedClusvsAmpMaxTClus(0),
  fHistAmpMaxTClusNotMatchingClus(0),
  fHistEtavsPhiMaxTClusNotMatchingClus(0),
  fHistIdxMaxClusvsIdxMaxTClus(0),
  fHistPhiMaxClusvsPhiMaxTClus(0),
  fHistEtaMaxClusvsEtaMaxTClus(0),
  fHistTOFmaxClusvsTimeMaxTClus(0),
  fHistEmatchedClusvsAmpMatchedTClus(0),
  fHistEmatchedClus(0),
  fHistEmaxMatchedClus(0),
  fHistAmpL1TimeSum(0),
  fHistAmpMaxL1TimeSum(0),
  fHistAmpMaxL1TimeSumVScent(0),
  fHistAmpFastORvsAmpL1TimeSum(0),
  fHistAmpFastOR(0),
  fHistAmpMaxFastOR(0),
  fHistTimeFastOR(0),
  fHistEtavsPhiFastOR(0),
  fHistEtavsPhiMaxFastOR(0),
  fHistTimeDispFastOR(0),
  fHistTimevsL0TimeFastOR(0),
  fHistNtimesFastOR(0),
  fHistEcells(0),
  fHistEmaxCell(0),
  fHistTOFvsEcells(0),
  fHistTOFvsEcellsC(0),
  fHistEmaxCellvsAmpFastOR(0),
  fCaloClustersName("CaloClusters"),
  fTriggerClustersName("triggerClusters"),
  fMinCutL0Amp(-1),
  fMaxCutL0Amp(999),
  fMinCutClusEnergy(-1),
  fMaxCutClusEnergy(999),
  fTimeCutOn(0),
  fMinL0Time(-20),
  fMaxL0Time(20),
  fMinClusTime(-1),
  fMaxClusTime(1),
  fCheckDeadClusters(0),
  fPedestal(0),
  fLoadPed(0),
  fOCDBpath(""),
  fMinDistanceFromBadTower(0.5),
  fClusterizer(0),
  fTriggerClusterizer(0),
  fRun(-1)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSATR::~AliAnalysisTaskSATR()
{
  delete fOutput;
}

//________________________________________________________________________
void AliAnalysisTaskSATR::UserCreateOutputObjects()
{
  // Create histograms
  
  if (fClusterizer)
    fCaloClustersName = fClusterizer->GetNewClusterArrayName();
  
  if (fTriggerClusterizer)
    fTriggerClustersName = fTriggerClusterizer->GetNewClusterArrayName();
  
  fOutput = new TList();
  fOutput->SetOwner();

  fHistEclus = new TH1F("fHistEclus","Energy Cluster spectrum", Ebins, Elow, Eup);
  fHistEclus->GetXaxis()->SetTitle("Cluster Energy [GeV]");
  fHistEclus->GetYaxis()->SetTitle("counts");
  
  fHistEmaxClus = new TH1F("fHistEmaxClus","Maximum Cluster Energy per event spectrum",Ebins, Elow, Eup);
  fHistEmaxClus->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
  fHistEmaxClus->GetYaxis()->SetTitle("counts");
  
  fHistEtavsPhiMaxClus = new TH2I("fHistEtavsPhiMaxClus","nEta vs nPhi of maximum cluster per event", nEtabins, nEtalow, nEtaup, nPhibins, nPhilow, nPhiup);
  fHistEtavsPhiMaxClus->GetXaxis()->SetTitle("nEta");
  fHistEtavsPhiMaxClus->GetYaxis()->SetTitle("nPhi");
	
  fHistEtavsEmaxClus = new TH2F("fHistEtavsEmaxClus","nEta vs Energy of maximum cluster per event", nEtabins, nEtalow, nEtaup, Ebins, Elow, Eup);
  fHistEtavsEmaxClus->GetYaxis()->SetTitle("Energy [GeV]");
  fHistEtavsEmaxClus->GetXaxis()->SetTitle("nEta");
  
  fHistPhivsEmaxClus = new TH2F("fHistPhivsEmaxClus","nPhi vs Energy of maximum cluster per event", nPhibins, nPhilow, nPhiup, Ebins, Elow, Eup);
  fHistPhivsEmaxClus->GetYaxis()->SetTitle("Energy [GeV]");
  fHistPhivsEmaxClus->GetXaxis()->SetTitle("nPhi");
  
  fHistTOFvsEclus = new TH2F("fHistTOFvsEclus","TOF vs Energy of clusters", TOFbins, TOFlow, TOFup, Ebins, Elow, Eup);
  fHistTOFvsEclus->GetXaxis()->SetTitle("TOF [s]");
  fHistTOFvsEclus->GetYaxis()->SetTitle("Energy [GeV]");
	
  fHistTOFvsEclusC = new TH2F("fHistTOFvsEclusC","TOF vs Energy of clusters (corrected)", TOFbins, TOFlow, TOFup, Ebins, Elow, Eup);
  fHistTOFvsEclusC->GetXaxis()->SetTitle("TOF");
  fHistTOFvsEclusC->GetYaxis()->SetTitle("E [GeV]");
  
  fHistNcellsvsEclus = new TH2F("fHistNcellsvsEclus","N_{cells} vs Energy of clusters", 20, 0, 20, Ebins, Elow, Eup);
  fHistNcellsvsEclus->GetXaxis()->SetTitle("N_{cells}");
  fHistNcellsvsEclus->GetYaxis()->SetTitle("E [GeV]");
  
  fHistAmpTClus = new TH1F("fHistAmpTClus","L0 Trigger Candidate Amplitude spectrum", L0Ampbins, L0Amplow, L0Ampup * 4);
  fHistAmpTClus->GetXaxis()->SetTitle("L0 Trigger Candidate Amplitude");
  fHistAmpTClus->GetYaxis()->SetTitle("counts");
	
  fHistAmpMaxTClus = new TH1F("fHistAmpMaxTClus","Maximum L0 Trigger Candidate Amplitude per event spectrum", L0Ampbins, L0Amplow, L0Ampup * 4);
  fHistAmpMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Amplitude");
  fHistAmpMaxTClus->GetYaxis()->SetTitle("counts");
  
  fHistEtavsPhiMaxTClus = new TH2I("fHistEtavsPhiMaxTClus","nEta vs nPhi of maximum L0 Trigger Candidate per event", nEtabins, nEtalow, nEtaup, nPhibins, nPhilow, nPhiup);
  fHistEtavsPhiMaxTClus->GetXaxis()->SetTitle("nEta");
  fHistEtavsPhiMaxTClus->GetYaxis()->SetTitle("nPhi");
  
  fHistEmaxClusvsAmpMaxTClus = new TH2F("fHistEmaxClusvsAmpMaxTClus","Maximum Cluster Energy vs Maximum L0 Trigger Candidate Amplitude", L0Ampbins, L0Amplow, L0Ampup * 4, Ebins, Elow, Eup);
  fHistEmaxClusvsAmpMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Amplitude");
  fHistEmaxClusvsAmpMaxTClus->GetYaxis()->SetTitle("Maximum Cluster Energy [GeV]");
  
  fHistEmaxClusvsAmpMatchedTClus = new TH2F("fHistEmaxClusvsAmpMatchedTClus","Maximum Cluster Energy vs Matched L0 Trigger Candidate Amplitude", L0Ampbins, L0Amplow, L0Ampup * 4, Ebins, Elow, Eup);
  fHistEmaxClusvsAmpMatchedTClus->GetXaxis()->SetTitle("Matched L0 Trigger Candidate Amplitude");
  fHistEmaxClusvsAmpMatchedTClus->GetYaxis()->SetTitle("Maximum Cluster Energy [GeV]");

  fHistEmaxClusNotMatchingTClus = new TH1F("fHistEmaxClusNotMatchingTClus","Maximum Cluster Energy per event not matching any fake ALTRO data", Ebins, Elow, Eup);
  fHistEmaxClusNotMatchingTClus->GetXaxis()->SetTitle("Maximum Cluster Energy [GeV]");
  fHistEmaxClusNotMatchingTClus->GetYaxis()->SetTitle("counts");

  fHistEtavsPhiMaxClusNotMatchingTClus = new TH2I("fHistEtavsPhiMaxClusNotMatchingTClus","nEta vs nPhi of maximum cluster per event not matching any fake ALTRO data", nEtabins, nEtalow, nEtaup, nPhibins, nPhilow, nPhiup);
  fHistEtavsPhiMaxClusNotMatchingTClus->GetXaxis()->SetTitle("nEta");
  fHistEtavsPhiMaxClusNotMatchingTClus->GetYaxis()->SetTitle("nPhi");
  
  fHistEmatchedClusvsAmpMaxTClus = new TH2F("fHistEmatchedClusvsAmpMaxTClus","Matched Cluster Energy vs Maximum L0 Trigger Candidate Amplitude", L0Ampbins, L0Amplow, L0Ampup * 4, Ebins, Elow, Eup);
  fHistEmatchedClusvsAmpMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Amplitude");
  fHistEmatchedClusvsAmpMaxTClus->GetYaxis()->SetTitle("Matched Cluster Energy [GeV]");
  
  fHistAmpMaxTClusNotMatchingClus = new TH1F("fHistAmpMaxTClusNotMatchingClus","Maximum L0 Trigger Candidate Amplitude not matching any FEE data", L0Ampbins, L0Amplow, L0Ampup * 4);
  fHistAmpMaxTClusNotMatchingClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Amplitude");
  fHistAmpMaxTClusNotMatchingClus->GetYaxis()->SetTitle("counts");
  
  fHistEtavsPhiMaxTClusNotMatchingClus = new TH2I("fHistEtavsPhiMaxTClusNotMatchingClus","nEta vs nPhi of maximum L0 Trigger Candidate Amplitude per event not matching any FEE data", nEtabins, nEtalow, nEtaup, nPhibins, nPhilow, nPhiup);
  fHistEtavsPhiMaxTClusNotMatchingClus->GetXaxis()->SetTitle("nEta");
  fHistEtavsPhiMaxTClusNotMatchingClus->GetYaxis()->SetTitle("nPhi");
  
  fHistIdxMaxClusvsIdxMaxTClus = new TH2I("fHistIdxMaxClusvsIdxMaxTClus","Maximum Cluster Index vs Maximum L0 Trigger Candidate Index", Indexesbins, Indexeslow, Indexesup, Indexesbins, Indexeslow, Indexesup);
  fHistIdxMaxClusvsIdxMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Index");
  fHistIdxMaxClusvsIdxMaxTClus->GetYaxis()->SetTitle("Maximum Cluster Index");
  
  fHistPhiMaxClusvsPhiMaxTClus = new TH2I("fHistPhiMaxClusvsPhiMaxTClus","Maximum Cluster nPhi vs Maximum L0 Trigger Candidate nPhi", nPhibins, nPhilow, nPhiup, nPhibins, nPhilow, nPhiup);
  fHistPhiMaxClusvsPhiMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate nPhi");
  fHistPhiMaxClusvsPhiMaxTClus->GetYaxis()->SetTitle("Maximum Cluster nPhi");

  fHistEtaMaxClusvsEtaMaxTClus = new TH2I("fHistEtaMaxClusvsEtaMaxTClus","Maximum Cluster nEta vs Maximum L0 Trigger Candidate nEta", nEtabins, nEtalow, nEtaup, nEtabins, nEtalow, nEtaup);
  fHistEtaMaxClusvsEtaMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate nEta");
  fHistEtaMaxClusvsEtaMaxTClus->GetYaxis()->SetTitle("Maximum Cluster nEta");
  
  fHistTOFmaxClusvsTimeMaxTClus = new TH2F("fHistTOFmaxClusvsTimeMaxTClus","Maximum Cluster TOF vs Maximum L0 Trigger Candidate Time", L0Timebins, L0Timelow, L0Timeup, TOFbins, TOFlow, TOFup);
  fHistTOFmaxClusvsTimeMaxTClus->GetXaxis()->SetTitle("Maximum L0 Trigger Candidate Time");
  fHistTOFmaxClusvsTimeMaxTClus->GetYaxis()->SetTitle("Maximum Cluster TOF");

  fHistEmatchedClusvsAmpMatchedTClus = new TH2F("fHistEmatchedClusvsAmpMatchedTClus","Matched Cluster Energy vs Matched L0 Trigger Candidate Amplitude", L0Ampbins, L0Amplow, L0Ampup * 4, Ebins, Elow, Eup);
  fHistEmatchedClusvsAmpMatchedTClus->GetXaxis()->SetTitle("Matched Trigger Candidate Amplitude");
  fHistEmatchedClusvsAmpMatchedTClus->GetYaxis()->SetTitle("Matched Cluster Energy [GeV]");
  
  fHistEmatchedClus = new TH1F("fHistEmatchedClus","Cluster Energy spectrum that matched L0 Trigger Candidate",Ebins, Elow, Eup);
  fHistEmatchedClus->GetXaxis()->SetTitle("Cluster Energy [GeV]");
  fHistEmatchedClus->GetYaxis()->SetTitle("counts");
  
  fHistEmaxMatchedClus = new TH1F("fHistEmaxMatchedClus","Maximum Cluster Energy per event that matched L0 Trigger Candidate",Ebins, Elow, Eup);
  fHistEmaxMatchedClus->GetXaxis()->SetTitle("Cluster Energy [GeV]");
  fHistEmaxMatchedClus->GetYaxis()->SetTitle("counts");
  
  fHistAmpL1TimeSum = new TH1F("fHistAmpL1TimeSum","L1 Time Sum Amplitude spectrum",L1Ampbins, L1Amplow, L1Ampup);
  fHistAmpL1TimeSum->GetXaxis()->SetTitle("L1 Time Sum Amplitude");
  fHistAmpL1TimeSum->GetYaxis()->SetTitle("counts");
  
  fHistAmpMaxL1TimeSum = new TH1F("fHistAmpMaxL1TimeSum","L1 Time Sum Amplitude per event spectrum",L1Ampbins, L1Amplow, L1Ampup);
  fHistAmpMaxL1TimeSum->GetXaxis()->SetTitle("Maximum L1 Time Sum Amplitude");
  fHistAmpMaxL1TimeSum->GetYaxis()->SetTitle("counts");

  fHistAmpMaxL1TimeSumVScent = new TH2F("fHistAmpMaxL1TimeSumVScent","L1 Time Sum Amplitude vs. centrality", 100, 0, 100, L1Ampbins, L1Amplow, L1Ampup);
  fHistAmpMaxL1TimeSumVScent->GetXaxis()->SetTitle("Centrality [%]");
  fHistAmpMaxL1TimeSumVScent->GetYaxis()->SetTitle("Maximum L1 Time Sum Amplitude");

  fHistAmpFastORvsAmpL1TimeSum = new TH2F("fHistAmpFastORvsAmpL1TimeSum","FastOR Amplitude vs L1 Time Sum Amplitude",L0Ampbins, L0Amplow, L0Ampup,L1Ampbins, L1Amplow, L1Ampup);
  fHistAmpFastORvsAmpL1TimeSum->GetXaxis()->SetTitle("FastOR Amplitude");
  fHistAmpFastORvsAmpL1TimeSum->GetYaxis()->SetTitle("L1 Time Sum Amplitude");
  
  fHistAmpFastOR = new TH1F("fHistAmpFastOR","FastOR Amplitude spectrum",L0Ampbins, L0Amplow, L0Ampup);
  fHistAmpFastOR->GetXaxis()->SetTitle("FastOR Amplitude");
  fHistAmpFastOR->GetYaxis()->SetTitle("counts");
	
  fHistAmpMaxFastOR = new TH1F("fHistAmpMaxFastOR","Maximum FastOR Amplitude per event spectrum",L0Ampbins, L0Amplow, L0Ampup);
  fHistAmpMaxFastOR->GetXaxis()->SetTitle("max FastOR amplitude");
  fHistAmpMaxFastOR->GetYaxis()->SetTitle("counts");
  
  fHistTimeFastOR = new TH1F("fHistTimeFastOR","FastOR L0 Time distribution", L0Timebins, L0Timelow, L0Timeup);
  fHistTimeFastOR->GetXaxis()->SetTitle("FastOR L0 Time");
  fHistTimeFastOR->GetYaxis()->SetTitle("counts");

  fHistEtavsPhiFastOR = new TH2I("fHistEtavsPhiFastOR","FastOR gCol vs gRow", ColTrgbins, ColTrglow, ColTrgup,  RowTrgbins, RowTrglow, RowTrgup);
  fHistEtavsPhiFastOR->GetXaxis()->SetTitle("FastOR gCol");
  fHistEtavsPhiFastOR->GetYaxis()->SetTitle("FastOR gRow");
	
  fHistEtavsPhiMaxFastOR = new TH2I("fHistEtavsPhiMaxFastOR","Maximum FastOR gCol vs gRow", ColTrgbins, ColTrglow, ColTrgup,  RowTrgbins, RowTrglow, RowTrgup);
  fHistEtavsPhiMaxFastOR->GetXaxis()->SetTitle("Maximum FastOR gCol");
  fHistEtavsPhiMaxFastOR->GetYaxis()->SetTitle("Maximum FastOR gRow");
  
  fHistTimeDispFastOR = new TH1F("fHistTimeDispFastOR", "Dispersion from the maximum FastOR of L0 Time per event", L0Timebins * 2, -L0Timeup + .5, L0Timeup - .5);
  fHistTimeDispFastOR->GetXaxis()->SetTitle("Disperion");
  fHistTimeDispFastOR->GetYaxis()->SetTitle("counts");
  
  fHistTimevsL0TimeFastOR = new TH2F("fHistTimevsL0TimeFastOR", "FastOR Time vs FastOR L0 Time", L0Timebins, L0Timelow, L0Timeup, L0Timebins, L0Timelow, L0Timeup);
  fHistTimevsL0TimeFastOR->GetXaxis()->SetTitle("Time");
  fHistTimevsL0TimeFastOR->GetYaxis()->SetTitle("L0 Time");
  
  fHistNtimesFastOR = new TH1I("fHistNtimesFastOR", "FastOR NL0Times distribution", 5, 0, 5);
  fHistNtimesFastOR->GetXaxis()->SetTitle("NL0Times");
  fHistNtimesFastOR->GetYaxis()->SetTitle("counts");
  
  fHistEcells = new TH1F("fHistEcells","Cell Energy spectrum", Ebins, Elow, Eup);
  fHistEcells->GetXaxis()->SetTitle("Energy [GeV]");
  fHistEcells->GetYaxis()->SetTitle("counts");
  
  fHistEmaxCell = new TH1F("fHistEmaxCell","Maximum Cell Energy per event spectrum", Ebins, Elow, Eup);
  fHistEmaxCell->GetXaxis()->SetTitle("Maximum Cell Energy [GeV]");
  fHistEmaxCell->GetYaxis()->SetTitle("counts");
  
  fHistTOFvsEcells = new TH2F("fHistTOFvsEcells","TOF vs Energy of cells", TOFbins, TOFlow, TOFup, Ebins, Elow, Eup);
  fHistTOFvsEcells->GetXaxis()->SetTitle("TOF [s]");
  fHistTOFvsEcells->GetYaxis()->SetTitle("Energy [GeV]");
	
  fHistTOFvsEcellsC = new TH2F("fHistTOFvsEcellsC","TOF vs Energy of cells (corrected)", TOFbins, TOFlow, TOFup, Ebins, Elow, Eup);
  fHistTOFvsEcellsC->GetXaxis()->SetTitle("TOF");
  fHistTOFvsEcellsC->GetYaxis()->SetTitle("E [GeV]");
  
  fHistEmaxCellvsAmpFastOR = new TH2F("fHistEmaxCellvsAmpFastOR","Maximum Cell Energy vs Maximum FastOR Amplitude", L0Ampbins, L0Amplow, L0Ampup, Ebins, Elow, Eup);
  fHistEmaxCellvsAmpFastOR->GetXaxis()->SetTitle("Maximum FastOR Amplitude");
  fHistEmaxCellvsAmpFastOR->GetYaxis()->SetTitle("Maximum Cell Energy [GeV]");
  
  fOutput->Add(fHistEclus);
  fOutput->Add(fHistEmaxClus);
  fOutput->Add(fHistEtavsPhiMaxClus);
  fOutput->Add(fHistEtavsEmaxClus);
  fOutput->Add(fHistPhivsEmaxClus);
  fOutput->Add(fHistTOFvsEclus);
  fOutput->Add(fHistTOFvsEclusC);
  fOutput->Add(fHistNcellsvsEclus);
  fOutput->Add(fHistAmpTClus);
  fOutput->Add(fHistAmpMaxTClus);
  fOutput->Add(fHistEtavsPhiMaxTClus);
  fOutput->Add(fHistEmaxClusvsAmpMaxTClus);
  fOutput->Add(fHistEmaxClusvsAmpMatchedTClus);
  fOutput->Add(fHistEmaxClusNotMatchingTClus);
  fOutput->Add(fHistEtavsPhiMaxClusNotMatchingTClus);
  fOutput->Add(fHistEmatchedClusvsAmpMaxTClus);
  fOutput->Add(fHistAmpMaxTClusNotMatchingClus);
  fOutput->Add(fHistEtavsPhiMaxTClusNotMatchingClus);
  fOutput->Add(fHistIdxMaxClusvsIdxMaxTClus);
  fOutput->Add(fHistPhiMaxClusvsPhiMaxTClus);
  fOutput->Add(fHistEtaMaxClusvsEtaMaxTClus);
  fOutput->Add(fHistTOFmaxClusvsTimeMaxTClus);
  fOutput->Add(fHistEmatchedClusvsAmpMatchedTClus);
  fOutput->Add(fHistEmatchedClus);
  fOutput->Add(fHistEmaxMatchedClus);
  fOutput->Add(fHistAmpL1TimeSum);
  fOutput->Add(fHistAmpMaxL1TimeSum);
  fOutput->Add(fHistAmpMaxL1TimeSumVScent);
  fOutput->Add(fHistAmpFastORvsAmpL1TimeSum);
  fOutput->Add(fHistAmpFastOR);
  fOutput->Add(fHistAmpMaxFastOR);
  fOutput->Add(fHistTimeFastOR);
  fOutput->Add(fHistEtavsPhiFastOR);
  fOutput->Add(fHistEtavsPhiMaxFastOR);
  fOutput->Add(fHistTimeDispFastOR);
  fOutput->Add(fHistTimevsL0TimeFastOR);
  fOutput->Add(fHistNtimesFastOR);
  fOutput->Add(fHistEcells);
  fOutput->Add(fHistEmaxCell);
  fOutput->Add(fHistTOFvsEcells);
  fOutput->Add(fHistTOFvsEcellsC);
  fOutput->Add(fHistEmaxCellvsAmpFastOR);
	
  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskSATR::Init() 
{
  if (fRun <= 0)
    return;
  
  AliCDBManager *cdb = AliCDBManager::Instance();
  
  if (!fPedestal && fLoadPed) {
    if (!cdb->IsDefaultStorageSet() && !fOCDBpath.IsNull())
      cdb->SetDefaultStorage(fOCDBpath);
    cdb->SetRun(fRun);
    AliCDBEntry *entry = static_cast<AliCDBEntry*>(AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals"));
    if (entry) 
      fPedestal =  static_cast<AliCaloCalibPedestal*>(entry->GetObject());
  }
}

//________________________________________________________________________
void AliAnalysisTaskSATR::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  // Create pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event){ 
    AliError("ERROR: Could not retrieve event"); 
    return; 
  }
  
  fRun = event->GetRunNumber();

  Float_t cent = 99;

  AliCentrality *centobj = InputEvent()->GetCentrality();
  if (centobj) {
    cent = centobj->GetCentralityPercentile("V0M");
  }
  
  Init();
  
  AliEMCALGeometry *fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  if (!fGeom) {
    AliError("Couldn't get geometry. Returning...");
    return;
  }
  
  TClonesArray *caloClusters = dynamic_cast<TClonesArray*>(event->FindListObject(fCaloClustersName));  

  if (!caloClusters){
    AliError("Troubles trying to load clusters!");
    caloClusters = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  }
 
  if(!caloClusters){
    AliError("Cannot get the CaloClusters!");
    return;
  }
  
  Int_t NumberOfCaloClusters = caloClusters->GetEntriesFast();
  
  TClonesArray *triggerClusters = dynamic_cast<TClonesArray*>(event->FindListObject(fTriggerClustersName));  
 
  Int_t NumberOfTriggerClusters = 0;
  
  if(triggerClusters){
    NumberOfTriggerClusters = triggerClusters->GetEntriesFast();
  }
  
  Float_t clusterTime_corr;
  if(event->GetBunchCrossNumber() % 4 < 2 )
    clusterTime_corr = -0.0000001;
  else 
    clusterTime_corr = 0;
  
  Float_t maxL0amp = -1;
  
  for (Int_t i = 0; i < NumberOfTriggerClusters; i++) {
    AliVCluster* trgCluster = dynamic_cast<AliVCluster*>(triggerClusters->At(i));
    
    Float_t L0amp = trgCluster->E();

    if (fCheckDeadClusters && trgCluster->GetDistanceToBadChannel() < 3)
      continue;
    
    if (fTimeCutOn && (trgCluster->GetTOF() < fMinL0Time || trgCluster->GetTOF() > fMaxL0Time))
      continue;
    
    if (maxL0amp < L0amp)
      maxL0amp = L0amp;
    } // trigger cluster loop
  
  if (maxL0amp > fMaxCutL0Amp || maxL0amp < fMinCutL0Amp) {
    AliWarning(Form("Event skipped because maximum trigger cluster amplitude = %f out of limits (%f, %f)", maxL0amp, fMinCutL0Amp, fMaxCutL0Amp));
    return; 
  }
  
  Float_t maxe = -1;
  
  for (Int_t iClusters = 0; iClusters < NumberOfCaloClusters; iClusters++) {
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(caloClusters->At(iClusters));
    if (!cluster){
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) 
      continue;
    
    if (fCheckDeadClusters && cluster->GetDistanceToBadChannel() < 3) {
      continue;
    }
    
    if (fTimeCutOn && (cluster->GetTOF() + clusterTime_corr < fMinClusTime || cluster->GetTOF() + clusterTime_corr > fMaxClusTime))
      continue;
    
    if (maxe < cluster->E()){
      maxe = cluster->E();
    }
  } //cluster loop 
    
  if (maxe > fMaxCutClusEnergy || maxL0amp < fMinCutClusEnergy) {
    AliWarning(Form("Event skipped because maximum cluster energy = %f out of limits (%f, %f)", maxe, fMaxCutClusEnergy, fMinCutClusEnergy));
    return;
  }
  
  AliVCaloTrigger *triggers = event->GetCaloTrigger("EMCAL");
  
  Float_t maxL0FastORamp = -1;
  Float_t maxL0FastORtime = -1;
  
  if (triggers && triggers->GetEntries() > 0) {
    
    triggers->Reset();
    Float_t L0FastORamp = 0;
    
    while (triggers->Next()) {
      
      triggers->GetAmplitude(L0FastORamp);
      
      if (L0FastORamp < 0)
        continue;
      
      if (maxL0FastORamp < L0FastORamp) {
	
        Int_t ntimes = 0;
        
        triggers->GetNL0Times(ntimes);
        
        if (fTimeCutOn && ntimes < 1) {
	  AliWarning(Form("FastOR removed from analysis because did not contribute to L0 trigger (no L0 time information)"));
          continue;
	}
        
        maxL0FastORamp = L0FastORamp;
        
        Int_t trgtimes[25];
        triggers->GetL0Times(trgtimes);
        
        maxL0FastORtime = trgtimes[0];
        
      }
    }
  }
  
  maxL0FastORamp = -1;
  Int_t maxL0FastORrow = -1;
  Int_t maxL0FastORcol = -1;
  Float_t maxL1amp = -1;
  Float_t TimeBinAmp = 0;
  
  if (triggers && triggers->GetEntries() > 0) {
    
    triggers->Reset();
    Float_t L0FastORamp = 0;
    Int_t L1amp = 0;
    
    while (triggers->Next()) {
      
      triggers->GetAmplitude(L0FastORamp);
      
      if (L0FastORamp < 0)
        continue;
      
      Int_t ntimes = 0;
      
      triggers->GetNL0Times(ntimes);
      
      triggers->GetTime(TimeBinAmp);
      
      fHistNtimesFastOR->Fill(ntimes);
      
      if (ntimes > 0) {
        Int_t trgtimes[25];
        triggers->GetL0Times(trgtimes);
        
        Int_t mintime = trgtimes[0];
        Int_t maxtime = trgtimes[0];
        
        for (Int_t i = 0; i < ntimes; ++i) {
          if (trgtimes[i] < mintime)
            mintime = trgtimes[i];
          if (maxtime < trgtimes[i])
            maxtime = trgtimes[i];
          
          fHistTimeFastOR->Fill(trgtimes[i]);
          fHistTimevsL0TimeFastOR->Fill(TimeBinAmp, trgtimes[i]);
          fHistTimeDispFastOR->Fill(trgtimes[i] - maxL0FastORtime);
        }
        
        if (fTimeCutOn && ((fMinL0Time > mintime) || (fMaxL0Time < maxtime)))
          continue;
      }
      
      if (fTimeCutOn && ntimes < 1)
        continue;
      
      Int_t gCol=0, gRow=0;
      triggers->GetPosition(gCol, gRow);
      
      Int_t find = -1;
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(gCol,gRow,find);
      
      if (find < 0)
	continue;
      
      Int_t cidx[4] = {-1};
      Bool_t ret = fGeom->GetCellIndexFromFastORIndex(find, cidx);
      
      if (!ret)
	continue;

      Int_t nSupMod = 0, nModule = 0, nIphi = 0, nIeta = 0, iphi = 0, ieta = 0;      

      Bool_t deadCluster = kFALSE;
      
      if (fCheckDeadClusters && fPedestal) {
        for (Int_t i = 0; i < 4; i++){
	  fGeom->GetCellIndex(cidx[i], nSupMod, nModule, nIphi, nIeta);
	  fGeom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta);
	  
	  Double_t d = fPedestal->GetDeadMap(nSupMod)->GetBinContent(ieta,iphi);
	  if (d == AliCaloCalibPedestal::kDead || d == AliCaloCalibPedestal::kHot){
	    AliWarning(Form("Dead/hot FastOR removed from analysis"));
	    deadCluster = kTRUE;
          }
	}
      }
      
      if (deadCluster)
	continue;

      fHistAmpFastOR->Fill(L0FastORamp);
      fHistEtavsPhiFastOR->Fill(gCol, gRow);
      
      if (maxL0FastORamp < L0FastORamp) {
	maxL0FastORamp = L0FastORamp;
	maxL0FastORcol = gCol;
	maxL0FastORrow = gRow;
      }
      
      triggers->GetL1TimeSum(L1amp);
      fHistAmpL1TimeSum->Fill(L1amp);
      fHistAmpFastORvsAmpL1TimeSum->Fill(L0FastORamp,L1amp);
      
      if (maxL1amp < L1amp) 
	maxL1amp = L1amp;
  
    }
  }
  
  if (maxL0FastORamp > -1) {
    fHistAmpMaxFastOR->Fill(maxL0FastORamp); 
    fHistEtavsPhiMaxFastOR->Fill(maxL0FastORcol, maxL0FastORrow);
  }
  
  if (maxL1amp > -1) {
    fHistAmpMaxL1TimeSum->Fill(maxL1amp);
    fHistAmpMaxL1TimeSumVScent->Fill(cent, maxL1amp);
  }
  
  maxL0amp = 0;
  Float_t EmatchMaxL0 = -1;
  Int_t maxL0index = -1;
  Float_t EmatchL0 = -1;
  AliVCluster* maxTcluster = 0;
 
  for (Int_t i = 0; i < NumberOfTriggerClusters; i++){
    AliVCluster* trgCluster = dynamic_cast<AliVCluster*>(triggerClusters->At(i));
    
    if (fCheckDeadClusters && trgCluster->GetDistanceToBadChannel() < 3) {
      AliWarning(Form("Trigger cluster %d removed from analysis because distance to bad channel = %f < 3", i, trgCluster->GetDistanceToBadChannel()));
      continue;
    }
    
    if (fTimeCutOn && (trgCluster->GetTOF() < fMinL0Time || trgCluster->GetTOF() > fMaxL0Time)) {
      AliWarning(Form("Trigger cluster %d removed from analysis because out of time %f, limits (%d, %d)", i, trgCluster->GetTOF(), fMinL0Time, fMaxL0Time));
      continue;
    }
    Float_t L0amp = trgCluster->E();
    
    fHistAmpTClus->Fill(L0amp);
      
    AliVCluster* matchedcluster = GetClusterFromId(caloClusters, trgCluster->GetID());
    
    if (matchedcluster) {
      EmatchL0 = matchedcluster->E();

      fHistEmatchedClus->Fill(EmatchL0);
      fHistEmatchedClusvsAmpMatchedTClus->Fill(L0amp, EmatchL0);
    }
    
    if (maxL0amp < L0amp) {
      maxL0amp = L0amp;
      maxL0index = trgCluster->GetID();
      EmatchMaxL0 = EmatchL0;
      maxTcluster = trgCluster;
    }
  }
  
  Int_t etaMaxTClus = 0;
  Int_t phiMaxTClus = 0;
  
  if (maxTcluster != 0) {
    fHistAmpMaxTClus->Fill(maxL0amp);
    fHistEmaxMatchedClus->Fill(EmatchMaxL0);
    
    if (fTriggerClusterizer) {
      Int_t nEtaDigitsSupMod = fGeom->GetNEta() * fGeom->GetNETAdiv(); // always 48?;
      Int_t nPhiDigitsSupMod = fGeom->GetNPhi() * fGeom->GetNPHIdiv(); // always 24?;
      
      Int_t nTRUPhi = 1;
      Int_t nTRUEta = 1;
      
      Int_t nEtaDigits = nEtaDigitsSupMod * fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
      Int_t nPhiDigits = nPhiDigitsSupMod * fGeom->GetNPhiSuperModule();    
      
      if (fTriggerClusterizer->GetTRUShift()) {
        nTRUPhi = fGeom->GetNPhiSuperModule() * 3;
        nTRUEta = fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
        nEtaDigits /= nTRUEta;
        nPhiDigits /= nTRUPhi;
      }
      
      Int_t nClusEtaNoShift = nEtaDigits / fTriggerClusterizer->GetnEta();
      Int_t nClusPhiNoShift = nPhiDigits / fTriggerClusterizer->GetnPhi();
      
      Int_t nClusters =  nClusEtaNoShift * nClusPhiNoShift * nTRUEta * nTRUPhi;
       
      etaMaxTClus = (maxTcluster->GetID() % nClusters) / (nClusPhiNoShift * nTRUPhi);
      phiMaxTClus = (maxTcluster->GetID() % nClusters) % (nClusPhiNoShift * nTRUPhi);
      
      fHistEtavsPhiMaxTClus->Fill(etaMaxTClus, phiMaxTClus);
    }
    
    AliVCluster *matchedClus = GetClusterFromId(caloClusters, maxTcluster->GetID());
    
    if (!matchedClus) {
      fHistAmpMaxTClusNotMatchingClus->Fill(maxL0amp);
      if (fClusterizer)
        fHistEtavsPhiMaxTClusNotMatchingClus->Fill(etaMaxTClus, phiMaxTClus);
    }
    else {
      fHistEmatchedClusvsAmpMaxTClus->Fill(maxL0amp, matchedClus->E());
    }
  }
    
  maxe = 0;
  Float_t maxClusterTime = 0;
	AliVCluster* maxcluster = 0;
  Int_t maxClusId = -1;
  Int_t maxiCluster = -1;
  
  for (Int_t iClusters = 0; iClusters < NumberOfCaloClusters; iClusters++) {
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(caloClusters->At(iClusters));
    if (!cluster){
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) 
      continue;
    
    if (fCheckDeadClusters && cluster->GetDistanceToBadChannel() < 3) {
      AliWarning(Form("Cluster %d removed from analysis because distance to bad channel = %f < 3", iClusters, cluster->GetDistanceToBadChannel()));
      continue;
    }
    
    if (fTimeCutOn && (cluster->GetTOF() < fMinL0Time || cluster->GetTOF() > fMaxL0Time)) {
      AliWarning(Form("Cluster %d removed from analysis because out of time %f, limits (%d, %d)", iClusters, cluster->GetTOF(), fMinL0Time, fMaxL0Time));
      continue;
    }
    
    if (maxe < cluster->E()){
      maxe = cluster->E();
      maxcluster = cluster;
      maxClusterTime = cluster->GetTOF() + clusterTime_corr;
      maxClusId = cluster->GetID();
      maxiCluster = iClusters;
    }
    
    fHistEclus->Fill(cluster->E());
    
    fHistTOFvsEclus->Fill(cluster->GetTOF(),cluster->E());
    
    fHistTOFvsEclusC->Fill(cluster->GetTOF() + clusterTime_corr,cluster->E());
    
    fHistNcellsvsEclus->Fill(cluster->GetNCells(),cluster->E());
    
  } //cluster loop 
  
  Int_t etaMaxClus = 0;
  Int_t phiMaxClus = 0;
  
  if (maxcluster != 0) {
    
    fHistEmaxClus->Fill(maxe);
    
    if (fClusterizer) {
      Int_t nEtaDigitsSupMod = fGeom->GetNEta() * fGeom->GetNETAdiv(); // always 48?;
      Int_t nPhiDigitsSupMod = fGeom->GetNPhi() * fGeom->GetNPHIdiv(); // always 24?;
      
      Int_t nTRUPhi = 1;
      Int_t nTRUEta = 1;
      
      Int_t nEtaDigits = nEtaDigitsSupMod * fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
      Int_t nPhiDigits = nPhiDigitsSupMod * fGeom->GetNPhiSuperModule();    
      
      if (fClusterizer->GetTRUShift()) {
        nTRUPhi = fGeom->GetNPhiSuperModule() * 3;
        nTRUEta = fGeom->GetNumberOfSuperModules() / fGeom->GetNPhiSuperModule();
        nEtaDigits /= nTRUEta;
        nPhiDigits /= nTRUPhi;
      }
      
      Int_t nClusEtaNoShift = nEtaDigits / fClusterizer->GetnEta();
      Int_t nClusPhiNoShift = nPhiDigits / fClusterizer->GetnPhi();
      
      Int_t nClusters =  nClusEtaNoShift * nClusPhiNoShift * nTRUEta * nTRUPhi;
      
      etaMaxClus = (maxcluster->GetID() % nClusters) / (nClusPhiNoShift * nTRUPhi);
      phiMaxClus = (maxcluster->GetID() % nClusters) % (nClusPhiNoShift * nTRUPhi);
      
      fHistEtavsPhiMaxClus->Fill(etaMaxClus, phiMaxClus);
      fHistEtavsEmaxClus->Fill(etaMaxClus, maxcluster->E());
      fHistPhivsEmaxClus->Fill(phiMaxClus, maxcluster->E());
    }
    
    AliVCluster *matchedTClus = GetClusterFromId(triggerClusters, maxcluster->GetID());
    
    if (!matchedTClus || (fTimeCutOn && (matchedTClus->GetTOF() < fMinL0Time || matchedTClus->GetTOF() > fMaxL0Time))) {
      fHistEmaxClusNotMatchingTClus->Fill(maxe);
      if (fClusterizer)
        fHistEtavsPhiMaxClusNotMatchingTClus->Fill(etaMaxClus, phiMaxClus);
    }
    else {
      fHistEmaxClusvsAmpMatchedTClus->Fill(matchedTClus->E(), maxe);
    }
  }
	
  if (maxTcluster != 0 && maxcluster != 0) {	
    fHistEmaxClusvsAmpMaxTClus->Fill(maxL0amp, maxe);    
    fHistTOFmaxClusvsTimeMaxTClus->Fill(maxTcluster->GetTOF(), maxcluster->GetTOF() + clusterTime_corr);
    fHistIdxMaxClusvsIdxMaxTClus->Fill(maxL0index, maxClusId);
    fHistPhiMaxClusvsPhiMaxTClus->Fill(phiMaxTClus, phiMaxClus);
    fHistEtaMaxClusvsEtaMaxTClus->Fill(etaMaxTClus, etaMaxClus);
  }
  
  AliVCaloCells *cells = event->GetEMCALCells();
  Int_t icells;
  Float_t maxCellEnergy = -1;
  for (icells = 0; icells < cells->GetNumberOfCells(); icells++) {
    if (fTimeCutOn && (cells->GetTime(icells) + clusterTime_corr < fMinClusTime || cells->GetTime(icells) + clusterTime_corr > fMaxClusTime))
      continue;
    
    fHistEcells->Fill(cells->GetAmplitude(icells));
    
    fHistTOFvsEcells->Fill(cells->GetTime(icells), cells->GetAmplitude(icells));
    fHistTOFvsEcellsC->Fill(cells->GetTime(icells) + clusterTime_corr, cells->GetAmplitude(icells));
    //cout << cells->GetTime(icells) << endl;
    
    if (maxCellEnergy < cells->GetAmplitude(icells)) 
      maxCellEnergy = cells->GetAmplitude(icells);
    
  }
  
  if (maxCellEnergy > -1) {
    fHistEmaxCell->Fill(maxCellEnergy);
  }
  
  if (maxL0FastORamp > -1 && maxCellEnergy > -1) {
    fHistEmaxCellvsAmpFastOR->Fill(maxL0FastORamp, maxCellEnergy);
  }
  
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
  
}


//________________________________________________________________________
void AliAnalysisTaskSATR::Terminate(Option_t *) 
{

}

AliVCluster* AliAnalysisTaskSATR::GetClusterFromId(TClonesArray *caloClusters, Int_t id)
{
  for (Int_t iClusters = 0; iClusters < caloClusters->GetEntriesFast(); iClusters++){
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(caloClusters->At(iClusters));
    if (!cluster){
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) 
      continue;
    
    if (id == cluster->GetID())
      return cluster;
  } //cluster loop 
  
  return 0;
}
