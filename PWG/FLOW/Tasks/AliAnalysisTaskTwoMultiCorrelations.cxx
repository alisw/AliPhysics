
/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//--------------------------------------------------------------------------------------//
// Analysis task for the computation of the multiparticle correlations for the flow     //
// harmonics v_1 to v_6. This version of the script compute the 2-, 4- and 6- particle  //
// correlations for all the useful combinations of these six harmonics. It can take     //
// Monte Carlo simulations data (e.g. HIJING), as well as the experimental Pb-Pb data   //
// taken by the ALICE experiment.                                                       //
// The method used to compute the multiparticle correlations is the Generic Framework   //
// based on Q-vectors. A setter lets open the possibility to cross-check the results    //
// with nested loops.                                                                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
// Version: 27.02.2019                                                                  //
//--------------------------------------------------------------------------------------//

#include "AliAnalysisTaskTwoMultiCorrelations.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "Riostream.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultSelection.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliMCVertex.h"
#include "AliAODMCParticle.h"
#include "TMath.h"
#include "TComplex.h"
#include "TFile.h"
#include "TList.h"
#include <vector>
#include "TH2I.h"
#include "TH1D.h"
#include "TH1I.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskTwoMultiCorrelations)

//######################################################################################//
// Mandatory methods for AliAnalysisTaskSE.
//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations() :
  AliAnalysisTaskSE(),
// General parameters of the analysis.
  fMainList(NULL),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fHighestFlowHarmonic(6),
  fMaxNumberOfParticlesInCorrelations(8),
  fProcessOnlyAOD(kFALSE),
  fProcessOnlyMC(kFALSE),
  fProcessBothMCandAOD(kFALSE),
  fComputeEtaGaps(kFALSE),
  fCrosscheckWithNestedLoops(kFALSE),
  fDoTDCorrelationHisto(kFALSE),
  fUseParticleWeights(kFALSE),
// Parameters related to the number of tracks.
  fMultiplicityList(NULL),
  fHistoInitialCentrality(NULL),
  fHistoFinalCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoInitialNumberOfTracksMain(NULL),
  fHistoInitialNumberOfTracksGlobal(NULL),
  fHistoFinalNumberOfTracksMain(NULL),
  fHistoFinalNumberOfTracksGlobal(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHistoInitialFilterCorrelations(NULL),
  fHistoFinalFilterCorrelations(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
  fNumberOfBinsHINOT(30000),
  fMaxBinHINOT(30000.),
  fNumberOfBinsHFNOT(30000),
  fMaxBinHFNOT(30000.),
  fNumberOfBinsHFC(5000),
  fMaxBinHFC(5000.),
// Parameters related to the event selection criteria.
  fEventSelectionList(NULL),
  fHistoInitialPVX(NULL),
  fHistoFinalPVX(NULL),
  fHistoInitialPVY(NULL),
  fHistoFinalPVY(NULL),
  fHistoInitialPVZ(NULL),
  fHistoFinalPVZ(NULL),
  fCutOnPVX(kFALSE),
  fCutOnPVY(kFALSE),
  fCutOnPVZ(kFALSE),
  fPVXMin(-44.),
  fPVXMax(44.),
  fPVYMin(-44),
  fPVYMax(44.),
  fPVZMin(-10.),
  fPVZMax(10),
// Parameters related to the track selection criteria.
  fTrackSelectionList(NULL),
  fHistoInitialPt(NULL),
  fHistoFinalPt(NULL),
  fHistoEfficiencyCorrection(NULL),
  fHistoInitialEta(NULL),
  fHistoFinalEta(NULL),
  fHistoInitialPhi(NULL),
  fHistoFinalPhi(NULL),
  fHistoInitialNumberOfTPC(NULL),
  fHistoFinalNumberOfTPC(NULL),
  fHistoInitialChiSquare(NULL),
  fHistoFinalChiSquare(NULL),
  fHistoInitialDCAxy(NULL),
  fHistoFinalDCAxy(NULL),
  fHistoInitialDCAz(NULL),
  fHistoFinalDCAz(NULL),
  fHistoInitialCharge(NULL),
  fHistoFinalCharge(NULL),
  fHistoInitialNumberOfITS(NULL),
  fHistoFinalNumberOfITS(NULL),
  fCutOnTDCorrelations(kFALSE),
  fCutOnPt(kFALSE),
  fCutOnEta(kFALSE),
  fCutOnNumberOfTPC(kFALSE),
  fCutOnChiSquarePInTPC(kFALSE),
  fCutOnDCAxy(kFALSE),
  fCutOnDCAz(kFALSE),
  fCutOnCharge(kFALSE),
  fMainFilter(768),
  fGlobalFilter(256),
  fMultiplicityMin(6),
  fMultiplicityMinA(0.),
  fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.),
  fMultiplicityMaxB(0.),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fNumberOfTPCMin(70),
  fChiSquarePInTPCMin(0.1),
  fChiSquarePInTPCMax(4.),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
  fCharge(0),
  fCutOnTracksMax(kFALSE),
  fNumberOfTracksMaxZero(0),
  fNumberOfTracksMaxFive(0),
  fNumberOfTracksMaxTen(0),
  fNumberOfTracksMaxTwenty(0),
  fNumberOfTracksMaxThirty(0),
  fNumberOfTracksMaxForty(0),
  fNumberOfTracksMaxFifty(0),
  fNumberOfTracksMaxSixty(0),
  fNumberOfTracksMaxSeventy(0),
// Parameters related to the multi-particle correlations.
  fMultiParticleCorrelationsList(NULL),
  fProfileTwoParticleCorrelations(NULL),
  fProfileFourParticleCorrelations(NULL),
  fProfileFourParticleCorrelationsCrossCheck(NULL),
  fProfileSixParticleCorrelations(NULL),
  fProfileTwoParticleCorrelationsNestedLoops(NULL),
  fProfileFourParticleCorrelationsNestedLoops(NULL),
// Parameters related to the 2-particle correlations with eta gaps.
  fTwoParticleCorrelationsWithEtaGapsList(NULL)
{
/* Dummy constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Initialise 'fQvectors' to zero.
  InitialiseArraysOfQvectors();

// Initialise the pointers of the TProfiles to NULL.
  InitialiseArraysOfTProfiles();

}

//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
  AliAnalysisTaskSE(name),
// General parameters of the analysis.
  fMainList(NULL),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fHighestFlowHarmonic(6),
  fMaxNumberOfParticlesInCorrelations(8),
  fProcessOnlyAOD(kFALSE),
  fProcessOnlyMC(kFALSE),
  fProcessBothMCandAOD(kFALSE),
  fComputeEtaGaps(kFALSE),
  fCrosscheckWithNestedLoops(kFALSE),
  fDoTDCorrelationHisto(kFALSE),
  fUseParticleWeights(kFALSE),
// Parameters related to the number of tracks.
  fMultiplicityList(NULL),
  fHistoInitialCentrality(NULL),
  fHistoFinalCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoInitialNumberOfTracksMain(NULL),
  fHistoInitialNumberOfTracksGlobal(NULL),
  fHistoFinalNumberOfTracksMain(NULL),
  fHistoFinalNumberOfTracksGlobal(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHistoInitialFilterCorrelations(NULL),
  fHistoFinalFilterCorrelations(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
  fNumberOfBinsHINOT(30000),
  fMaxBinHINOT(30000.),
  fNumberOfBinsHFNOT(30000),
  fMaxBinHFNOT(30000.),
  fNumberOfBinsHFC(5000),
  fMaxBinHFC(5000.),
// Parameters related to the event selection criteria.
  fEventSelectionList(NULL),
  fHistoInitialPVX(NULL),
  fHistoFinalPVX(NULL),
  fHistoInitialPVY(NULL),
  fHistoFinalPVY(NULL),
  fHistoInitialPVZ(NULL),
  fHistoFinalPVZ(NULL),
  fCutOnPVX(kFALSE),
  fCutOnPVY(kFALSE),
  fCutOnPVZ(kFALSE),
  fPVXMin(-44.),
  fPVXMax(44.),
  fPVYMin(-44),
  fPVYMax(44.),
  fPVZMin(-10.),
  fPVZMax(10.),
// Parameters related to the track selection criteria.
  fTrackSelectionList(NULL),
  fHistoInitialPt(NULL),
  fHistoFinalPt(NULL),
  fHistoEfficiencyCorrection(NULL),
  fHistoInitialEta(NULL),
  fHistoFinalEta(NULL),
  fHistoInitialPhi(NULL),
  fHistoFinalPhi(NULL),
  fHistoInitialNumberOfTPC(NULL),
  fHistoFinalNumberOfTPC(NULL),
  fHistoInitialChiSquare(NULL),
  fHistoFinalChiSquare(NULL),
  fHistoInitialDCAxy(NULL),
  fHistoFinalDCAxy(NULL),
  fHistoInitialDCAz(NULL),
  fHistoFinalDCAz(NULL),
  fHistoInitialCharge(NULL),
  fHistoFinalCharge(NULL),
  fHistoInitialNumberOfITS(NULL),
  fHistoFinalNumberOfITS(NULL),
  fCutOnTDCorrelations(kFALSE),
  fCutOnPt(kFALSE),
  fCutOnEta(kFALSE),
  fCutOnNumberOfTPC(kFALSE),
  fCutOnChiSquarePInTPC(kFALSE),
  fCutOnDCAxy(kFALSE),
  fCutOnDCAz(kFALSE),
  fCutOnCharge(kFALSE),
  fMainFilter(768),
  fGlobalFilter(256),
  fMultiplicityMin(6),
  fMultiplicityMinA(0.),
  fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.),
  fMultiplicityMaxB(0.),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fNumberOfTPCMin(70),
  fChiSquarePInTPCMin(0.1),
  fChiSquarePInTPCMax(4.),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
  fCharge(0),
  fCutOnTracksMax(kFALSE),
  fNumberOfTracksMaxZero(0),
  fNumberOfTracksMaxFive(0),
  fNumberOfTracksMaxTen(0),
  fNumberOfTracksMaxTwenty(0),
  fNumberOfTracksMaxThirty(0),
  fNumberOfTracksMaxForty(0),
  fNumberOfTracksMaxFifty(0),
  fNumberOfTracksMaxSixty(0),
  fNumberOfTracksMaxSeventy(0),
// Parameters related to the multi-particle correlations.
  fMultiParticleCorrelationsList(NULL),
  fProfileTwoParticleCorrelations(NULL),
  fProfileFourParticleCorrelations(NULL),
  fProfileFourParticleCorrelationsCrossCheck(NULL),
  fProfileSixParticleCorrelations(NULL),
  fProfileTwoParticleCorrelationsNestedLoops(NULL),
  fProfileFourParticleCorrelationsNestedLoops(NULL),
// Parameters related to the 2-particle correlations with eta gaps.
  fTwoParticleCorrelationsWithEtaGapsList(NULL)
{
/* Constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Create the mother list with the ownership of everything inside it.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE);

// Define the input and output slots.
  DefineOutput(1, TList::Class());

// Initialise 'fQvectors' to zero.
  InitialiseArraysOfQvectors();

// Initialise the pointers of the TProfiles to NULL.
  InitialiseArraysOfTProfiles();

/// TBA: how to deal with non-unit particle weights.
}

//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
/* Destructor of the class. */
  if (fMainList) {delete fMainList;}
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Define the outputs of the task at the beginning of the analysis. */
// Avoid name clashes.
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  // JEfficiency for NUA correction : DongJo
  fEfficiency = new AliJEfficiency();
  fEfficiency->SetMode(1) ; // 1: priod should work for you
  fEfficiency->SetDataPath( "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data" );
// Book all the lists.
  this->BookAllLists();

// Book the histograms in each daughter list.
  this->BookMultiplicityList();
  this->BookEventSelectionList();
  this->BookTrackSelectionList();
  this->BookMultiParticleCorrelationsList();
  this->BookTwoParticleCorrelationsWithEtaGapsList();

// Continue to avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1, fMainList);
  fFirstEvent = kTRUE;
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the analysis for each event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";

  AliAODEvent *currentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());  // Pointer to an AOD event.
  AliMCEvent  *currentMCEvent = MCEvent();                                  // Pointer to a Monte Carlo event.

// Select the analysis method according to the type of files.
  if ((Int_t)fProcessOnlyAOD + (Int_t)fProcessOnlyMC + (Int_t)fProcessBothMCandAOD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one 'fProcess' must be set to kTRUE.");
  }
  else if (fProcessOnlyAOD) {AnalyseAODevent(currentAODEvent);}
  else if (fProcessOnlyMC) {AnalyseMCevent(currentMCEvent);}
  else if (fProcessBothMCandAOD) {Fatal(sMethodName.Data(),"ERROR: TBA analysis with both AOD and MC files.");}

// PostData.
  PostData(1, fMainList);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
/* Save the outputs after the running over the events. */
// Access the mother list.
  fMainList = (TList*)GetOutputData(1);
  if (!fMainList) {exit(1);}

// Create the output file and save the mother list inside.
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fMainList->Write(fMainList->GetName(),TObject::kSingleKey);
  delete outputFile;
}

//######################################################################################//
// Methods called in the constructors.
//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfQvectors()
{
/* Initialise to zero all the elements in 'fQvectors'. */
  for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    }
  }
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfTProfiles()
{
/* Initialise to NULL the pointers to the TProfiles inside an array. */
  for (Int_t i = 0; i < 6; i++) {fProfileTwoParticleCorrelationsWithEtaGaps[i] = NULL;}
}

//######################################################################################//
// Methods called in 'UserCreateOutputObjects'.
//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists in fMainList. */
// Check if the mother list exists.
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()";
  if (!fMainList) {Fatal(sMethodName.Data(), "ERROR: 'fMainList' does not exist.");}

// Daughter list for the histograms containing the number of tracks.
  fMultiplicityList = new TList();
  fMultiplicityList->SetName("fMultiplicityList");
  fMultiplicityList->SetOwner(kTRUE);
  fMainList->Add(fMultiplicityList);

// Daughter list for the histograms containing the event selection criteria.
  fEventSelectionList = new TList();
  fEventSelectionList->SetName("fEventSelectionList");
  fEventSelectionList->SetOwner(kTRUE);
  fMainList->Add(fEventSelectionList);

// Daughter list for the histograms containing the track selection criteria.
  fTrackSelectionList = new TList();
  fTrackSelectionList->SetName("fTrackSelectionList");
  fTrackSelectionList->SetOwner(kTRUE);
  fMainList->Add(fTrackSelectionList);

// Daughter list with the multi-particle correlations.
  fMultiParticleCorrelationsList = new TList();
  fMultiParticleCorrelationsList->SetName("fMultiParticleCorrelationsList");
  fMultiParticleCorrelationsList->SetOwner(kTRUE);
  fMainList->Add(fMultiParticleCorrelationsList);

// Daughter list with the 2-p correlations computed with eta gaps.
  fTwoParticleCorrelationsWithEtaGapsList = new TList();
  fTwoParticleCorrelationsWithEtaGapsList->SetName("fTwoParticleCorrelationsWithEtaGapsList");
  fTwoParticleCorrelationsWithEtaGapsList->SetOwner(kTRUE);
  fMainList->Add(fTwoParticleCorrelationsWithEtaGapsList);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookMultiplicityList()
{
/* Book the histograms containing the number of tracks. */
// Distributions of the centrality.
  fHistoInitialCentrality = new TH1D("fHistoInitialCentrality",
    "Centrality before the event selection", 100, 0., 100.);
  fHistoInitialCentrality->SetStats(kTRUE);
  fHistoInitialCentrality->GetXaxis()->SetTitle("Centrality percentile");
  fHistoInitialCentrality->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoInitialCentrality);

  fHistoFinalCentrality = new TH1D("fHistoFinalCentrality",
    "Centrality after the event selection", 100, 0., 100.);
  fHistoFinalCentrality->SetStats(kTRUE);
  fHistoFinalCentrality->GetXaxis()->SetTitle("Centrality percentile");
  fHistoFinalCentrality->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoFinalCentrality);

// Distribution of the initial, unfiltered number of tracks.
  fHistoInitialNumberOfTracks = new TH1I("fHistoInitialNumberOfTracks",
    "Initial number of tracks without filter", fNumberOfBinsHINOT, 0., fMaxBinHINOT);
  fHistoInitialNumberOfTracks->SetStats(kTRUE);
  fHistoInitialNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fHistoInitialNumberOfTracks->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoInitialNumberOfTracks);

// Distributions of the number of tracks for each filter bit before the removal of the high multiplicity outliers.
  fHistoInitialNumberOfTracksMain = new TH1I("fHistoInitialNumberOfTracksMain",
    "Number of tracks before the outliers cuts for the main filter bit",
    fNumberOfBinsHINOT, 0., fMaxBinHINOT);
  fHistoInitialNumberOfTracksMain->SetStats(kTRUE);
  fHistoInitialNumberOfTracksMain->GetXaxis()->SetTitle("Number of tracks");
  fHistoInitialNumberOfTracksMain->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoInitialNumberOfTracksMain);

  fHistoInitialNumberOfTracksGlobal = new TH1I("fHistoInitialNumberOfTracksGlobal",
    "Number of tracks before the outlier cuts for the global filter bit",
    fNumberOfBinsHINOT, 0., fMaxBinHINOT);
  fHistoInitialNumberOfTracksGlobal->SetStats(kTRUE);
  fHistoInitialNumberOfTracksGlobal->GetXaxis()->SetTitle("Number of tracks");
  fHistoInitialNumberOfTracksGlobal->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoInitialNumberOfTracksGlobal);

// Distributions of the number of tracks for each filter bit after the removal of the high multiplicity outliers.
  fHistoFinalNumberOfTracksMain = new TH1I("fHistoFinalNumberOfTracksMain",
    "Number of tracks after the outliers cuts for the main filter bit",
    fNumberOfBinsHFNOT, 0., fMaxBinHFNOT);
  fHistoFinalNumberOfTracksMain->SetStats(kTRUE);
  fHistoFinalNumberOfTracksMain->GetXaxis()->SetTitle("Number of tracks");
  fHistoFinalNumberOfTracksMain->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoFinalNumberOfTracksMain);

  fHistoFinalNumberOfTracksGlobal = new TH1I("fHistoFinalNumberOfTracksGlobal",
    "Number of tracks after the outliers cuts for the global filter bit",
    fNumberOfBinsHFNOT, 0., fMaxBinHFNOT);
  fHistoFinalNumberOfTracksGlobal->SetStats(kTRUE);
  fHistoFinalNumberOfTracksGlobal->GetXaxis()->SetTitle("Number of tracks");
  fHistoFinalNumberOfTracksGlobal->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoFinalNumberOfTracksGlobal);

// Distribution of the number of tracks used in the analysis itself.
  fHistoFinalNumberOfTracks = new TH1I("fHistoFinalNumberOfTracks",
    "Number of tracks used in the analysis for the main filter bit",
    fNumberOfBinsHFNOT, 0., fMaxBinHFNOT);
  fHistoFinalNumberOfTracks->SetStats(kTRUE);
  fHistoFinalNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fHistoFinalNumberOfTracks->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoFinalNumberOfTracks);

// 2D correlation histograms between two filters.
  if (fDoTDCorrelationHisto) // Manually enabled only to decide the cuts, as it is memory-greedy.
  {
    fHistoInitialFilterCorrelations = new TH2I("fHistoInitialFilterCorrelations",
      "2D multiplicity correlations between filters before the outliers cuts",
      fNumberOfBinsHFC, 0., fMaxBinHFC, fNumberOfBinsHFC, 0., fMaxBinHFC);
    fHistoInitialFilterCorrelations->SetStats(kTRUE);
    fHistoInitialFilterCorrelations->GetXaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fGlobalFilter));
    fHistoInitialFilterCorrelations->GetYaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fMainFilter));
    fHistoInitialFilterCorrelations->SetMarkerSize(0.5);
    fHistoInitialFilterCorrelations->SetMarkerColor(kBlue);
    fHistoInitialFilterCorrelations->SetMarkerStyle(kFullCircle);
    fMultiplicityList->Add(fHistoInitialFilterCorrelations);

    fHistoFinalFilterCorrelations = new TH2I("fHistoFinalFilterCorrelations",
      "2D multiplicity correlations between filters after the outliers removal",
      fNumberOfBinsHFC, 0., fMaxBinHFC, fNumberOfBinsHFC, 0., fMaxBinHFC);
    fHistoFinalFilterCorrelations->SetStats(kTRUE);
    fHistoFinalFilterCorrelations->GetXaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fGlobalFilter));
    fHistoFinalFilterCorrelations->GetYaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fMainFilter));
    fHistoFinalFilterCorrelations->SetMarkerSize(0.5);
    fHistoFinalFilterCorrelations->SetMarkerColor(kBlue);
    fHistoFinalFilterCorrelations->SetMarkerStyle(kFullCircle);
    fMultiplicityList->Add(fHistoFinalFilterCorrelations);
  }
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookEventSelectionList()
{
/* Book the control histograms for the event selection criteria. */
// Distributions of the x-position of the PV.
  fHistoInitialPVX = new TH1D("fHistoInitialPVX",
    "PV_{x} before the event selection", 1000, -20., 20.);
  fHistoInitialPVX->SetStats(kTRUE);
  fHistoInitialPVX->GetXaxis()->SetTitle("PV_{x} [cm]");
  fHistoInitialPVX->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVX);

  fHistoFinalPVX = new TH1D("fHistoFinalPVX",
    "PV_{x} after the event selection", 1000, -20., 20.);
  fHistoFinalPVX->SetStats(kTRUE);
  fHistoFinalPVX->GetXaxis()->SetTitle("PV_{x} [cm]");
  fHistoFinalPVX->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVX);

// Distributions of the y-position of the PV.
  fHistoInitialPVY = new TH1D("fHistoInitialPVY",
    "PV_{y} before the event selection", 1000, -20., 20.);
  fHistoInitialPVY->SetStats(kTRUE);
  fHistoInitialPVY->GetXaxis()->SetTitle("PV_{y} [cm]");
  fHistoInitialPVY->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVY);

  fHistoFinalPVY = new TH1D("fHistoFinalPVY",
    "PV_{y} after the event selection", 1000, -20., 20.);
  fHistoFinalPVY->SetStats(kTRUE);
  fHistoFinalPVY->GetXaxis()->SetTitle("PV_{y} [cm]");
  fHistoFinalPVY->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVY);

// Distributions of the z-position of the PV.
  fHistoInitialPVZ = new TH1D("fHistoInitialPVZ",
    "PV_{z} before the event selection", 1000, -20., 20.);
  fHistoInitialPVZ->SetStats(kTRUE);
  fHistoInitialPVZ->GetXaxis()->SetTitle("PV_{z} [cm]");
  fHistoInitialPVZ->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVZ);

  fHistoFinalPVZ = new TH1D("fHistoFinalPVZ",
    "PV_{z} after the event selection", 1000, -20., 20.);
  fHistoFinalPVZ->SetStats(kTRUE);
  fHistoFinalPVZ->GetXaxis()->SetTitle("PV_{z} [cm]");
  fHistoFinalPVZ->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVZ);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookTrackSelectionList()
{
/* Book the control histograms for the track selection criteria. */
// Distributions of the transverse momentum and the efficiency correction.
  fHistoInitialPt = new TH1D("fHistoInitialPt",
    "p_{T} before the track selection", 1000, 0., 20.);
  fHistoInitialPt->SetStats(kTRUE);
  fHistoInitialPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoInitialPt->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialPt);

  fHistoFinalPt = new TH1D("fHistoFinalPt",
    "p_{T} after the track selection", 1000, 0., 20.);
  fHistoFinalPt->SetStats(kTRUE);
  fHistoFinalPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoFinalPt->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalPt);

  fHistoEfficiencyCorrection = new TH1D("fHistoEfficiencyCorrection",
    "Efficiency correction before inversion", 1000, 0., 5.);
  fHistoEfficiencyCorrection->SetStats(kTRUE);
  fHistoEfficiencyCorrection->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoEfficiencyCorrection->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoEfficiencyCorrection);

// Distributions of the pseudorapidity.
  fHistoInitialEta = new TH1D("fHistoInitialEta",
    "#eta before the track selection", 1000, -5.5, 5.5);
  fHistoInitialEta->SetStats(kTRUE);
  fHistoInitialEta->GetXaxis()->SetTitle("#eta");
  fHistoInitialEta->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialEta);

  fHistoFinalEta = new TH1D("fHistoFinalEta",
    "#eta after the track selection", 1000, -5.5, 5.5);
  fHistoFinalEta->SetStats(kTRUE);
  fHistoFinalEta->GetXaxis()->SetTitle("#eta");
  fHistoFinalEta->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalEta);

// Distributions of the azimuthal angles.
  fHistoInitialPhi = new TH1D("fHistoInitialPhi",
    "#phi before the track selection", 1000, 0., 6.3);
  fHistoInitialPhi->SetStats(kTRUE);
  fHistoInitialPhi->GetXaxis()->SetTitle("#phi");
  fHistoInitialPhi->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialPhi);

  fHistoFinalPhi = new TH1D("fHistoFinalPhi",
    "#phi after the track selection", 1000, 0., 6.3);
  fHistoFinalPhi->SetStats(kTRUE);
  fHistoFinalPhi->GetXaxis()->SetTitle("#phi");
  fHistoFinalPhi->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalPhi);

// Distributions of the number of TPC clusters before the track selection.
  fHistoInitialNumberOfTPC = new TH1I("fHistoInitialNumberOfTPC",
    "Number of TPC clusters before the track selection", 1000, 0., 170.);
  fHistoInitialNumberOfTPC->SetStats(kTRUE);
  fHistoInitialNumberOfTPC->GetXaxis()->SetTitle("Number of TPC clusters");
  fHistoInitialNumberOfTPC->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialNumberOfTPC);

  fHistoFinalNumberOfTPC = new TH1I("fHistoFinalNumberOfTPC",
    "Number of TPC clusters after the track selection", 1000, 0., 170.);
  fHistoFinalNumberOfTPC->SetStats(kTRUE);
  fHistoFinalNumberOfTPC->GetXaxis()->SetTitle("Number of TPC clusters");
  fHistoFinalNumberOfTPC->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalNumberOfTPC);

// Distributions of the chi^2 of the track momentum in TPC.
  fHistoInitialChiSquare = new TH1D("fHistoInitialChiSquare",
    "#chi^{2}/NDF in the TPC before the track selection", 1000, 0., 20.);
  fHistoInitialChiSquare->SetStats(kTRUE);
  fHistoInitialChiSquare->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoInitialChiSquare->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialChiSquare);

  fHistoFinalChiSquare = new TH1D("fHistoFinalChiSquare",
    "#chi^{2}/NDF in the TPC after the track selection", 1000, 0., 20.);
  fHistoFinalChiSquare->SetStats(kTRUE);
  fHistoFinalChiSquare->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoFinalChiSquare->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalChiSquare);

// Distributions of the xy-coordinate of the DCA.
  fHistoInitialDCAxy = new TH1D("fHistoInitialDCAxy",
    "DCA_{xy} before the track selection", 1000, 0., 10.);
  fHistoInitialDCAxy->SetStats(kTRUE);
  fHistoInitialDCAxy->GetXaxis()->SetTitle("DCA_{xy} [cm]");
  fHistoInitialDCAxy->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialDCAxy);

  fHistoFinalDCAxy = new TH1D("fHistoFinalDCAxy",
    "DCA_{xy} after the track selection", 1000, 0., 10.);
  fHistoFinalDCAxy->SetStats(kTRUE);
  fHistoFinalDCAxy->GetXaxis()->SetTitle("DCA_{xy} [cm]");
  fHistoFinalDCAxy->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalDCAxy);

// Distributions of the z-coordinate of the DCA.
  fHistoInitialDCAz = new TH1D("fHistoInitialDCAz",
    "DCA_{z} before the track selection", 1000, 0., 10.);
  fHistoInitialDCAz->SetStats(kTRUE);
  fHistoInitialDCAz->GetXaxis()->SetTitle("DCA_{z} [cm]");
  fHistoInitialDCAz->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialDCAz);

  fHistoFinalDCAz = new TH1D("fHistoFinalDCAz",
    "DCA_{z} after the track selection", 1000, 0., 10.);
  fHistoFinalDCAz->SetStats(kTRUE);
  fHistoFinalDCAz->GetXaxis()->SetTitle("DCA_{z} [cm]");
  fHistoFinalDCAz->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalDCAz);

// Distributions of the electric charge of the tracks.
  fHistoInitialCharge = new TH1I("fHistoInitialCharge",
    "Electric charge before the track selection", 2, -2, 2);
  fHistoInitialCharge->SetStats(kTRUE);
  fHistoInitialCharge->GetXaxis()->SetTitle("Charge");
  fHistoInitialCharge->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialCharge);

  fHistoFinalCharge = new TH1I("fHistoFinalCharge",
    "Electric charge after the track selection", 2, -2, 2);
  fHistoFinalCharge->SetStats(kTRUE);
  fHistoFinalCharge->GetXaxis()->SetTitle("Charge");
  fHistoFinalCharge->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalCharge);

// Distributions of the number of clusters in the ITS.
  fHistoInitialNumberOfITS = new TH1I("fHistoInitialNumberOfITS",
    "Number of ITS clusters before the track selection", 1000, 0., 20.);
  fHistoInitialNumberOfITS->SetStats(kTRUE);
  fHistoInitialNumberOfITS->GetXaxis()->SetTitle("Number of ITS clusters");
  fHistoInitialNumberOfITS->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoInitialNumberOfITS);

  fHistoFinalNumberOfITS = new TH1I("fHistoFinalNumberOfITS",
    "Number of ITS clusters after the track selection", 1000, 0., 20.);
  fHistoFinalNumberOfITS->SetStats(kTRUE);
  fHistoFinalNumberOfITS->GetXaxis()->SetTitle("Number of ITS clusters");
  fHistoFinalNumberOfITS->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalNumberOfITS);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookMultiParticleCorrelationsList()
{
/* Book the TProfiles with the multiparticle correlations. */
// 2-particle correlations.
  fProfileTwoParticleCorrelations = new TProfile("fProfileTwoParticleCorrelations",
    "#LT#LT2#GT#GT_{n,-n}", 6, 0., 6.);
  fProfileTwoParticleCorrelations->SetStats(kTRUE);
  fProfileTwoParticleCorrelations->Sumw2();
  fProfileTwoParticleCorrelations->GetXaxis()->SetTitle("n");
  fMultiParticleCorrelationsList->Add(fProfileTwoParticleCorrelations);

  if (fCrosscheckWithNestedLoops)
  {
    fProfileTwoParticleCorrelationsNestedLoops = new TProfile("fProfileTwoParticleCorrelationsNestedLoops",
      "#LT#LT2#GT#GT_{n,-n} with nested loops", 6, 0., 6.);
    fProfileTwoParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileTwoParticleCorrelationsNestedLoops->Sumw2();
    fProfileTwoParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("n");
    fMultiParticleCorrelationsList->Add(fProfileTwoParticleCorrelationsNestedLoops);
  }

// 4-particle correlations.
  fProfileFourParticleCorrelations = new TProfile("fProfileFourParticleCorrelations",
    "#LT#LT4#GT#GT_{m,n,-m,-n}", 21, 0., 21.);
  fProfileFourParticleCorrelations->SetStats(kTRUE);
  fProfileFourParticleCorrelations->Sumw2();
  fProfileFourParticleCorrelations->GetXaxis()->SetTitle("(m,n)");
  fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelations);

  fProfileFourParticleCorrelationsCrossCheck = new TProfile("fProfileFourParticleCorrelationsCrossCheck",
    "Cross-check of #LT#LT4#GT#GT_{m,n,-m,-n}", 15, 0., 15.);
  fProfileFourParticleCorrelationsCrossCheck->SetStats(kTRUE);
  fProfileFourParticleCorrelationsCrossCheck->Sumw2();
  fProfileFourParticleCorrelationsCrossCheck->GetXaxis()->SetTitle("(m,n)");
  fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelationsCrossCheck);

  if (fCrosscheckWithNestedLoops)
  {
    fProfileFourParticleCorrelationsNestedLoops = new TProfile("fProfileFourParticleCorrelationsNestedLoops",
      "#LT#LT4#GT#GT_{m,n,-m,-n} with nested loops", 21, 0., 21.);
    fProfileFourParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileFourParticleCorrelationsNestedLoops->Sumw2();
    fProfileFourParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("(m,n)");
    fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelationsNestedLoops);
  }

// 6-particle correlations.
  fProfileSixParticleCorrelations = new TProfile("fProfileSixParticleCorrelations",
    "#LT#LT6#GT#GT_{l,m,n,-l,-m,-n}", 20, 0., 20.);
  fProfileSixParticleCorrelations->SetStats(kTRUE);
  fProfileSixParticleCorrelations->Sumw2();
  fProfileSixParticleCorrelations->GetXaxis()->SetTitle("(l,m,n)");
  fMultiParticleCorrelationsList->Add(fProfileSixParticleCorrelations);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookTwoParticleCorrelationsWithEtaGapsList()
{
/* Book the TProfiles for the 2-p correlations with eta gaps if selected. */
  if (fComputeEtaGaps)
  {
    Double_t EtaGaps[11] = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.}; // Choice of values for the eta gap.
    for (Int_t i = 0; i < 6; i++)
    {
      fProfileTwoParticleCorrelationsWithEtaGaps[i] = new TProfile("", "", 11, 0., 11.);
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->SetName(Form("fProfileTwoParticleCorrelationsWithEtaGaps_v%d", i+1));
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->SetTitle(Form("#LT#LT2#GT#GT_{%d,-%d}", i+1, i+1));
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->SetStats(kTRUE);
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->Sumw2();
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->GetXaxis()->SetTitle("#eta gap");
      fTwoParticleCorrelationsWithEtaGapsList->Add(fProfileTwoParticleCorrelationsWithEtaGaps[i]);

    // Set bin labels.
      for (Int_t gap = 1; gap < 12; gap++)
      {
       fProfileTwoParticleCorrelationsWithEtaGaps[i]->GetXaxis()->SetBinLabel(gap, Form("%.1f", EtaGaps[gap-1]));
      }
    }
  }
}

//######################################################################################//
// Methods called in 'UserExec'.
//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *inputAODevent)
{
/* Execute the analysis for the given AOD event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *inputAODevent)";
  if (!inputAODevent) {Fatal(sMethodName.Data(), "ERROR: no AOD event found.");}

// Check if the event passes the event selection criteria.
  if(!ApplyEventSelectionAOD(inputAODevent)) {return;}

  // Efficiency : DongJo
  if( fFirstEvent ) {
    fEfficiency->SetRunNumber( inputAODevent->GetRunNumber() );
    fEfficiency->Load();
    fFirstEvent = kFALSE;
  }

// Obtain locally the centrality of the event for the efficiency correction.
  TString centralityEstimatorAOD = "centralityEstimatorAOD";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimatorAOD = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimatorAOD = "CL1";}
  AliMultSelection *multiplicitySelectionAOD = (AliMultSelection*)inputAODevent->FindListObject("MultSelection");
  if (!multiplicitySelectionAOD) {return;}

  Double_t eventCentralityAOD = multiplicitySelectionAOD->GetMultiplicityPercentile(Form("%s", centralityEstimatorAOD.Data()));

// Determine how many and which tracks pass the track selection.
  long long currentNumberOfTracks = inputAODevent->GetNumberOfTracks(); // Number of tracks before the track selection.
  long long finalNumberOfTracks = 0.; // Number of tracks remaining after the track selection.
  Int_t *isTrackSelected = new Int_t[currentNumberOfTracks](); // Flag to indicate if the track passed the selection (1) or not (0). 
  for (long long iTrack = 0; iTrack < currentNumberOfTracks; iTrack++)
  {
  // Check if the current track belongs to the main filter bit.
    AliAODTrack *inputTrack = dynamic_cast<AliAODTrack*>(inputAODevent->GetTrack(iTrack));
    if (!inputTrack) {continue;}
    if (!inputTrack->TestFilterBit(fMainFilter)) {continue;}

  // Check if the track passes the track selection criteria.
    if (ApplyTrackSelectionAOD(inputTrack))
    {
      finalNumberOfTracks++;
      isTrackSelected[iTrack] = 1;
    }
    else {isTrackSelected[iTrack] = 0;}
  }

// Check if the event has still enough tracks to have a meaningful event weight.
  if (finalNumberOfTracks < fMultiplicityMin) {return;}

// Fill the distribution of the final number of tracks.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);

// Prepare the observables for the analysis itself.
  Double_t *currentEta = new Double_t[finalNumberOfTracks](); // Pseudorapidity for the method of the eta gaps.
  Double_t *currentPhi = new Double_t[finalNumberOfTracks](); // Azimuthal angles.
  Double_t *currentParticleWeights = new Double_t[finalNumberOfTracks](); // Particle weights.
  Int_t finalIndex = 0; // New index of the track if it passed the selection.
  for (long long iiTrack = 0; iiTrack < currentNumberOfTracks; iiTrack++)
  {
    AliAODTrack *currentTrack = dynamic_cast<AliAODTrack*>(inputAODevent->GetTrack(iiTrack));
    if (!currentTrack) {continue;}
    if (!currentTrack->TestFilterBit(fMainFilter)) {continue;}

    if (isTrackSelected[iiTrack] == 0) {continue;}  // The track has failed the track selection.

    currentEta[finalIndex] = currentTrack->Eta();
    currentPhi[finalIndex] = currentTrack->Phi();

    // Efficiency by DongJo
    Double_t pt = currentTrack->Pt();
    Double_t fCent = eventCentralityAOD;// this doesn't work eventCentrality; --> Fixed.
    Double_t effCorr = fEfficiency->GetCorrection( pt,0, fCent); //filterbitIndex 0 : TPCOnly 6: hybrid(this work for AOD86, i am not sure if it is work for new AOD)
    Double_t effInv = 1.0/effCorr;

    currentParticleWeights[finalIndex] = 1.;  // Unit particle weights are used by default.
    if (fUseParticleWeights) {currentParticleWeights[finalIndex] = effInv;}

    fHistoEfficiencyCorrection->Fill(effCorr);

    finalIndex++;
  }

// Do the analysis for the current event.
  CalculateQvectors(finalNumberOfTracks, currentPhi, currentParticleWeights);
  ComputeMultiparticleCorrelations(finalNumberOfTracks, currentPhi, currentParticleWeights);
  if (fComputeEtaGaps) {ComputeTwoParticleEtaGaps(finalNumberOfTracks, currentPhi, currentParticleWeights, currentEta);}

// Reset the variables before passing to the next event.
  currentNumberOfTracks = 0;
  finalNumberOfTracks = 0;
  delete [] isTrackSelected;
  delete [] currentEta;
  delete [] currentPhi;
  delete [] currentParticleWeights;
  finalIndex = 0;
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent(AliMCEvent *inputMCevent)
{
/* Execute the analysis for the provided MC event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent(AliMCEvent *aMCevent)";
  if (!inputMCevent) {Fatal(sMethodName.Data(), "ERROR: no MC event found.");}

// Check if the event passes the event selection criteria.
  if(!ApplyEventSelectionMC(inputMCevent)) {return;}

// Determine how many and which tracks pass the track selection.
  long long currentNumberOfTracks = inputMCevent->GetNumberOfTracks();  // Number of tracks before the track selection.
  long long finalNumberOfTracks = 0.; // Number of tracks remaining after the track selection.
  Int_t *isTrackSelected = new Int_t[currentNumberOfTracks](); // Flag to indicate if the track passed the selection (1) or not (0). 
  for (long long iTrack = 0; iTrack < currentNumberOfTracks; iTrack++)
  {
  // Check if the current track belongs to the main filter bit.
    AliAODMCParticle *inputTrack = dynamic_cast<AliAODMCParticle*>(inputMCevent->GetTrack(iTrack));
    if (!inputTrack) {continue;}

  // Check if the track passes the track selection criteria.
    if (ApplyTrackSelectionMC(inputTrack))
    {
      finalNumberOfTracks++;
      isTrackSelected[iTrack] = 1;
    }
    else {isTrackSelected[iTrack] = 0;}
  }

// Check if the event has still enough tracks to have a meaningful event weight.
  if (finalNumberOfTracks < fMultiplicityMin) {return;}

// Fill the distribution of the final number of tracks.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);

// Prepare the observables for the analysis itself.
  Double_t *currentEta = new Double_t[finalNumberOfTracks](); // Pseudorapidity for the method of the eta gaps.
  Double_t *currentPhi = new Double_t[finalNumberOfTracks](); // Azimuthal angles.
  Double_t *currentParticleWeights = new Double_t[finalNumberOfTracks](); // Particle weights.
  Int_t finalIndex = 0; // New index of the track if it passed the selection.
  for (long long iiTrack = 0; iiTrack < currentNumberOfTracks; iiTrack++)
  {
    AliAODMCParticle *currentTrack = dynamic_cast<AliAODMCParticle*>(inputMCevent->GetTrack(iiTrack));
    if (!currentTrack) {continue;}

    if (isTrackSelected[iiTrack] == 0) {continue;}  // The track has failed the track selection.

    currentEta[finalIndex] = currentTrack->Eta();
    currentPhi[finalIndex] = currentTrack->Phi();
    currentParticleWeights[finalIndex] = 1.;  // Unit particle weights are used by default.
    if (fUseParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA use of non-unit particle weights");}

    finalIndex++;
  }

// Do the analysis for the current event.
  CalculateQvectors(finalNumberOfTracks, currentPhi, currentParticleWeights);
  ComputeMultiparticleCorrelations(finalNumberOfTracks, currentPhi, currentParticleWeights);
  if (fComputeEtaGaps) {ComputeTwoParticleEtaGaps(finalNumberOfTracks, currentPhi, currentParticleWeights, currentEta);}

// Reset the variables before passing to the next event.
  currentNumberOfTracks = 0;
  finalNumberOfTracks = 0;
  delete [] isTrackSelected;
  delete [] currentEta;
  delete [] currentPhi;
  delete [] currentParticleWeights;
  finalIndex = 0;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelectionAOD(AliAODEvent *aAODevent)
{
/* Apply the event selection criteria on the given AOD event and return if the event passed it or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelectionAOD(AliAODEvent *aAODevent)";

// Set the detector for the centrality estimation as set in the task.
  TString centralityEstimator = "centralityEstimator";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Reject events not belonging to the current centrality range.
  AliMultSelection *multiplicitySelection = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if (!multiplicitySelection) {return kFALSE;}

  Double_t eventCentrality = multiplicitySelection->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));
  if ( (eventCentrality < fCentralityMin) || (eventCentrality >= fCentralityMax) ) {return kFALSE;}

// Get the observables used in the physics event selection.
  long long unfilteredInitialNumberOfTracks = aAODevent->GetNumberOfTracks(); // Unfiltered number of tracks before any selection.
  AliAODVertex *eventVertex = (AliAODVertex*)aAODevent->GetPrimaryVertex(); // 3D position of the PV.

// Fill the QA histograms before the event selection.
  fHistoInitialCentrality->Fill(eventCentrality);
  fHistoInitialNumberOfTracks->Fill(unfilteredInitialNumberOfTracks);
  fHistoInitialPVX->Fill(eventVertex->GetX());
  fHistoInitialPVY->Fill(eventVertex->GetY());
  fHistoInitialPVZ->Fill(eventVertex->GetZ());

// Apply the physics event selection criteria.
/// Cuts on the position of the PV.
  if (fCutOnPVX)
  {
    if ( (eventVertex->GetX() < fPVXMin) || (eventVertex->GetX() > fPVXMax) ) {return kFALSE;}
  }
  if (fCutOnPVY)
  {
    if ( (eventVertex->GetY() < fPVYMin) || (eventVertex->GetY() > fPVYMax) ) {return kFALSE;}
  }
  if (fCutOnPVZ)
  {
    if ( (eventVertex->GetZ() < fPVZMin) || (eventVertex->GetZ() > fPVZMax) ) {return kFALSE;}
  }

/// TBI: more physics criteria?

// Remove the high multiplicity outliers.
/// Get the number of tracks for the main and the global filter bits.
  long long numberOfTracks = aAODevent->GetNumberOfTracks();  // Total number of tracks in the given event.
  long long globalNumberOfTracks = 0; // Number of tracks in the global filter.
  long long mainNumberOfTracks = 0; // Number of tracks in the main filter bit.

  for (Int_t iTrack = 0; iTrack < numberOfTracks; iTrack++)
  {
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack));
    if (!aTrack) {continue;}

    if (aTrack->TestFilterBit(fGlobalFilter)) {globalNumberOfTracks++;}
    if (aTrack->TestFilterBit(fMainFilter)) {mainNumberOfTracks++;}
  }

/// Fill the histograms before the application of the cuts.
  fHistoInitialNumberOfTracksGlobal->Fill(globalNumberOfTracks);
  fHistoInitialNumberOfTracksMain->Fill(mainNumberOfTracks);
  if (fDoTDCorrelationHisto) {fHistoInitialFilterCorrelations->Fill(globalNumberOfTracks, mainNumberOfTracks);}

/// Apply the cuts.
  if (fCutOnTDCorrelations)
  {
    if ( (Double_t)mainNumberOfTracks < (fMultiplicityMinA*(Double_t)globalNumberOfTracks + fMultiplicityMinB) ) {return kFALSE;} // The number of tracks in the main filter is under the accepted band.
    if ( (Double_t)mainNumberOfTracks > (fMultiplicityMaxA*(Double_t)globalNumberOfTracks + fMultiplicityMaxB) ) {return kFALSE;} // The number of tracks in the main filter is above the accepted band.
  }

// Apply of the brute-force method to remove HM outliers?
  Int_t cutValueMaxNumberOfTracks = 0;  // Value of the cut on the maximum number of tracks.
  if (fCutOnTracksMax)  // If the cuts on the maximum numbers of tracks are enabled.
  {
  // Determine the value to cut depending on the centrality.
    if ( (eventCentrality >= 0.) && (eventCentrality < 5.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxZero;}
    else if ( (eventCentrality >= 5.) && (eventCentrality < 10.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFive;}
    else if ( (eventCentrality >= 10.) && (eventCentrality < 20.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTen;}
    else if ( (eventCentrality >= 20.) && (eventCentrality < 30.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTwenty;}
    else if ( (eventCentrality >= 30.) && (eventCentrality < 40.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxThirty;}
    else if ( (eventCentrality >= 40.) && (eventCentrality < 50.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxForty;}
    else if ( (eventCentrality >= 50.) && (eventCentrality < 60.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFifty;}
    else if ( (eventCentrality >= 60.) && (eventCentrality < 70.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSixty;}
    else if ( (eventCentrality >= 70.) && (eventCentrality < 80.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSeventy;}

  // Apply the cut.
    if (mainNumberOfTracks > cutValueMaxNumberOfTracks) {return kFALSE;}
  }

// Fill the histograms after the application of the cuts.
// No need to get the new numbers of tracks as the events qualified as outliers are normally removed with the last cuts.
  fHistoFinalNumberOfTracksGlobal->Fill(globalNumberOfTracks);
  fHistoFinalNumberOfTracksMain->Fill(mainNumberOfTracks);
  if (fDoTDCorrelationHisto) {fHistoFinalFilterCorrelations->Fill(globalNumberOfTracks, mainNumberOfTracks);}

// Fill the QA histograms after the event selection.
  fHistoFinalCentrality->Fill(eventCentrality);
  fHistoFinalPVX->Fill(eventVertex->GetX());
  fHistoFinalPVY->Fill(eventVertex->GetY());
  fHistoFinalPVZ->Fill(eventVertex->GetZ());

// Reset the variables.
  numberOfTracks = 0;
  globalNumberOfTracks = 0;
  mainNumberOfTracks = 0;
  cutValueMaxNumberOfTracks = 0;

  return kTRUE;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelectionMC(AliMCEvent *aMCevent)
{
/* Apply the event selection criteria on the given MC event and return if the event passed it or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelectionMC(AliMCEvent *aMCevent)";

// Set the detector for the centrality estimation as set in the task.
  TString centralityEstimator = "centralityEstimator";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

/// Reject events not belonging to the current centrality range.
  Double_t eventCentrality = 0.;  // Centrality of the given event.
  if (!fProcessOnlyMC)  // Centrality has no sense for Kine events.
  {
    AliMultSelection *multiplicitySelection = (AliMultSelection*)aMCevent->FindListObject("MultSelection");
    if (!multiplicitySelection) {return kFALSE;}

    eventCentrality = multiplicitySelection->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));
    if ( (eventCentrality < fCentralityMin) || (eventCentrality >= fCentralityMax) ) {return kFALSE;}
  }

// Get the observables used in the physics event selection.
  long long unfilteredInitialNumberOfTracks = aMCevent->GetNumberOfTracks();  // Unfiltered number of tracks before any selection.
  AliMCVertex *eventVertex = (AliMCVertex*)aMCevent->GetPrimaryVertex();  // 3D position of the PV.

// Fill the QA histograms before the event selection.
  fHistoInitialCentrality->Fill(eventCentrality);
  fHistoInitialNumberOfTracks->Fill(unfilteredInitialNumberOfTracks);
  fHistoInitialPVX->Fill(eventVertex->GetX());
  fHistoInitialPVY->Fill(eventVertex->GetY());
  fHistoInitialPVZ->Fill(eventVertex->GetZ());

// Apply the physics event selection criteria.
/// Cuts on the position of the PV.
  if (fCutOnPVX)
  {
    if ( (eventVertex->GetX() < fPVXMin) || (eventVertex->GetX() > fPVXMax) ) {return kFALSE;}
  }
  if (fCutOnPVY)
  {
    if ( (eventVertex->GetY() < fPVYMin) || (eventVertex->GetY() > fPVYMax) ) {return kFALSE;}
  }
  if (fCutOnPVZ)
  {
    if ( (eventVertex->GetZ() < fPVZMin) || (eventVertex->GetZ() > fPVZMax) ) {return kFALSE;}
  }

/// TBI: more physics criteria?

// Apply of the brute-force method to remove HM outliers?
// The 2d correlation method cannot be used with MC tracks as they have no notion of filter bits.
  long long numberOfTracks = aMCevent->GetNumberOfTracks();  // Total number of tracks in the given event.
  Int_t cutValueMaxNumberOfTracks = 0;  // Value of the cut on the maximum number of tracks.
  if (fCutOnTracksMax)  // If the cuts on the maximum numbers of tracks are enabled.
  {
  // Determine the value to cut depending on the centrality.
    if ( (eventCentrality >= 0.) && (eventCentrality < 5.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxZero;}
    else if ( (eventCentrality >= 5.) && (eventCentrality < 10.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFive;}
    else if ( (eventCentrality >= 10.) && (eventCentrality < 20.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTen;}
    else if ( (eventCentrality >= 20.) && (eventCentrality < 30.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTwenty;}
    else if ( (eventCentrality >= 30.) && (eventCentrality < 40.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxThirty;}
    else if ( (eventCentrality >= 40.) && (eventCentrality < 50.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxForty;}
    else if ( (eventCentrality >= 50.) && (eventCentrality < 60.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFifty;}
    else if ( (eventCentrality >= 60.) && (eventCentrality < 70.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSixty;}
    else if ( (eventCentrality >= 70.) && (eventCentrality < 80.) ) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSeventy;}

  // Apply the cut.
    if (numberOfTracks > cutValueMaxNumberOfTracks) {return kFALSE;}
  }

// Fill the QA histograms after the event selection.
  fHistoFinalCentrality->Fill(eventCentrality);
  fHistoFinalPVX->Fill(eventVertex->GetX());
  fHistoFinalPVY->Fill(eventVertex->GetY());
  fHistoFinalPVZ->Fill(eventVertex->GetZ());

// Reset the variables.
  eventCentrality = 0.;
  numberOfTracks = 0;
  cutValueMaxNumberOfTracks = 0;

  return kTRUE;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionAOD(AliAODTrack *aAODtrack)
{
/* Apply the track selection criteria and return if the track passed it or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)";

// Get the observables needed for the track selection.
  Double_t  pT = aAODtrack->Pt();                                               // Transverse momentum of the track.
  Double_t  eta = aAODtrack->Eta();                                             // Pseudorapidity of the track.
  Double_t  phi = aAODtrack->Phi();                                             // Azimuthal angle of the track.
  Int_t     numberOfTPCClusters = aAODtrack->GetTPCNcls();                      // Number of TPC clusters gone through by the track.
  Int_t     numberOfITSClusters = aAODtrack->GetITSNcls();                      // Number of ITS clusters gone through by the track.
  Double_t  chiSquareInTPC = (aAODtrack->GetTPCchi2())/(aAODtrack->GetNcls(1)); // Chi square in the TPC.
  Double_t  DCAx = aAODtrack->XAtDCA();                                         // x-coordinate of the DCA.
  Double_t  DCAy = aAODtrack->YAtDCA();                                         // y-coordinate of the DCA.
  Double_t  DCAz = aAODtrack->ZAtDCA();                                         // z-coordinate of the DCA.
  Double_t  DCAxy = TMath::Sqrt( (DCAx*DCAx) + (DCAy*DCAy) );                   // xy-coordinate of the DCA.
  Int_t     charge = aAODtrack->Charge();                                       // Electric charge.

// Fill the initial histograms before the track selection.
  fHistoInitialPt->Fill(pT);
  fHistoInitialEta->Fill(eta);
  fHistoInitialPhi->Fill(phi);
  fHistoInitialNumberOfTPC->Fill(numberOfTPCClusters);
  fHistoInitialChiSquare->Fill(chiSquareInTPC);
  fHistoInitialDCAxy->Fill(DCAxy);
  fHistoInitialDCAz->Fill(DCAz);
  fHistoInitialCharge->Fill(charge);
  fHistoInitialNumberOfITS->Fill(numberOfITSClusters);

// Apply the selection to the track.
  if (fCutOnPt)
  {
    if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;}
  }
  if (fCutOnEta)
  {
    if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;}
  }
  if (fCutOnNumberOfTPC)
  {
    if (numberOfTPCClusters < fNumberOfTPCMin) {return kFALSE;}
  }
  if (fCutOnChiSquarePInTPC)
  {
    if ( (chiSquareInTPC < fChiSquarePInTPCMin) || (chiSquareInTPC > fChiSquarePInTPCMax) ) {return kFALSE;}
  }
  if (fCutOnDCAxy)
  {
    if (DCAxy > fDCAxyMax) {return kFALSE;}
  }
  if (fCutOnDCAz)
  {
    if (DCAz > fDCAzMax) {return kFALSE;}
  }
  if (fCutOnCharge)
  {
    if (charge != fCharge) {return kFALSE;}
  }
  if (fCutOnNumberOfITS)
  {
    if (numberOfITSClusters < fNumberOfITSMin) {return kFALSE;}
  }

// Fill the histograms after the track selection.
  fHistoFinalPt->Fill(pT);
  fHistoFinalEta->Fill(eta);
  fHistoFinalPhi->Fill(phi);
  fHistoFinalNumberOfTPC->Fill(numberOfTPCClusters);
  fHistoFinalChiSquare->Fill(chiSquareInTPC);
  fHistoFinalDCAxy->Fill(DCAxy);
  fHistoFinalDCAz->Fill(DCAz);
  fHistoFinalCharge->Fill(charge);
  fHistoFinalNumberOfITS->Fill(numberOfITSClusters);

// Reset the variables.
  pT = 0.;
  eta = 0.;
  phi = 0.;
  numberOfTPCClusters = 0;
  numberOfITSClusters = 0;
  chiSquareInTPC = 0.;
  DCAx = 0.;
  DCAy = 0.;
  DCAz = 0.;
  DCAxy = 0.;
  charge = 0;

  return kTRUE;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionMC(AliAODMCParticle *aMCtrack)
{
/* Apply the track selection criteria to the given MC track and return if it passed or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionMC(AliAODMCParticle *aMCtrack)";

// Get the observables needed for the track selection.
  Double_t  pT = aMCtrack->Pt();          // Transverse momentum of the track.
  Double_t  eta = aMCtrack->Eta();        // Pseudorapidity of the track.
  Double_t  phi = aMCtrack->Phi();        // Azimuthal angle of the track.
  Int_t     charge = aMCtrack->Charge();  // Electric charge.

// Fill the initial histograms before the track selection.
  fHistoInitialPt->Fill(pT);
  fHistoInitialEta->Fill(eta);
  fHistoInitialPhi->Fill(phi);
  fHistoInitialCharge->Fill(charge);

// Apply the selection to the track.
  if (fCutOnPt)
  {
    if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;}
  }
  if (fCutOnEta)
  {
    if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;}
  }
  if (fCutOnCharge)
  {
    if (charge != fCharge) {return kFALSE;}
  }

// Fill the histograms after the track selection.
  fHistoFinalPt->Fill(pT);
  fHistoFinalEta->Fill(eta);
  fHistoFinalPhi->Fill(phi);
  fHistoFinalCharge->Fill(charge);

// Reset the variables.
  pT = 0.;
  eta = 0.;
  phi = 0.;
  charge = 0;

  return kTRUE;
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(long long numberOfParticles, Double_t angles[], Double_t pWeights[])
{
/* Calculate all the Q-vectors for the given arrays of azimuthal angles and particle weights. */
  Double_t iAngle = 0.; // Azimuthal angle of the current particle.
  Double_t iWeight = 0.;  // Particle weight of the current particle.
  Double_t iWeightToPowerP = 0.;  // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    }
  }

// Compute the Q-vectors.
  for (long long iTrack = 0; iTrack < numberOfParticles; iTrack++)
  {
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < 9; iPower++)
      {
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle), iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }
    }
  }
}

//======================================================================================//
TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
/* Simplify the use of the Q-vectors. */
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]); // Use that Q*(n,p) = Q(-n,p).
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::ComputeMultiparticleCorrelations(long long numberOfParticles, Double_t angles[], Double_t pWeights[])
{
/* Compute the considered 2-, 4- and 6-particle correlations for the harmonics v_1 to v_6. */
// Compute the 2-particle correlations.
  Int_t twoZerosArray[2] = {0}; // Zero array for the denominator.
  Double_t twoParticleDenominator = CalculateRecursion(2, twoZerosArray).Re();  // 2-particle event weight.

  Int_t twoParticleArray[2] = {0};  // Array to save the harmonics with the format: cos(nphi1 - nphi2).
  TComplex twoParticleComplex = TComplex(0., 0.); // Complex value for the 2-particle correlation.
  Double_t twoParticleValue = 0.; // Real part of the 2-particle correlation.
  Double_t iBinTwoMiddle = 0.;  // Index of the corresponding bin in the TProfile.
  for (Int_t n = 1; n < 7; n++)
  {
  // Compute the correlation.
    twoParticleArray[0] = n;
    twoParticleArray[1] = -1*n;
    twoParticleComplex = (CalculateRecursion(2, twoParticleArray))/twoParticleDenominator;
    twoParticleValue = twoParticleComplex.Re();

  // Fill the corresponding bin in the TProfile.
    iBinTwoMiddle = (1.*n) - 0.5;
    fProfileTwoParticleCorrelations->Fill(iBinTwoMiddle, twoParticleValue, twoParticleDenominator);
    fProfileTwoParticleCorrelations->GetXaxis()->SetBinLabel(n, Form("%d", n));

  // Cross-check with nested loops if needed.
    if (fCrosscheckWithNestedLoops)
    {
    // Fill the corresponding bin in the TProfile.
      ComputeTwoNestedLoops(numberOfParticles, twoParticleArray, angles, pWeights, fProfileTwoParticleCorrelationsNestedLoops, iBinTwoMiddle);
      fProfileTwoParticleCorrelationsNestedLoops->GetXaxis()->SetBinLabel(n, Form("%d", n));
    }

  // Reset of the variables for the next harmonic.
    twoParticleComplex = TComplex(0., 0.);
    twoParticleValue = 0.;
    iBinTwoMiddle = 0.;
  }

// Compute the 4-particle correlations.
  Int_t fourZerosArray[4] = {0};  // Zero array for the denominator.
  Double_t fourParticleDenominator = CalculateRecursion(4, fourZerosArray).Re();  // 4-particle event weight.

  Int_t fourParticleArray[4] = {0}; // Array to save the harmonics with the format: cos(mphi1 + nphi2 - mphi3 - nphi4).
  TComplex fourParticleComplex = TComplex(0., 0.);  // Complex value for the 4-particle correlation.
  Double_t fourParticleValue = 0.;  // Real part of the 4-particle correlation.
  Int_t iBinFour = 1; // Index of the corresponding bin in the TProfile.
  Int_t iBinFourMiddle = 0.;  // Middle of the corresponding bin in the TProfile.
  for(Int_t m = 1; m < 7; m++)
  {
    for (Int_t n = m; n < 7; n++)
    {
    // Compute the correlation.
      fourParticleArray[0] = m;
      fourParticleArray[1] = n;
      fourParticleArray[2] = -1*m;
      fourParticleArray[3] = -1*n;
      fourParticleComplex = (CalculateRecursion(4, fourParticleArray))/fourParticleDenominator;
      fourParticleValue = fourParticleComplex.Re();

    // Fill the corresponding bin in the TProfile.
      iBinFourMiddle = (1.*iBinFour) - 0.5;
      fProfileFourParticleCorrelations->Fill(iBinFourMiddle, fourParticleValue, fourParticleDenominator);
      fProfileFourParticleCorrelations->GetXaxis()->SetBinLabel(iBinFour, Form("(%d,%d)", m, n));

    // Cross-check with nested loops if needed.
      if (fCrosscheckWithNestedLoops)
      {
      // Fill the corresponding bin in the TProfile.
        ComputeFourNestedLoops(numberOfParticles, fourParticleArray, angles, pWeights, fProfileFourParticleCorrelationsNestedLoops, iBinFourMiddle);
        fProfileFourParticleCorrelationsNestedLoops->GetXaxis()->SetBinLabel(iBinFour, Form("(%d,%d)", m, n));
      }

    // Reset of the variables for the next pair of harmonics.
      iBinFour++;
      fourParticleComplex = TComplex(0., 0.);
      fourParticleValue = 0.;
      iBinFourMiddle = 0.;
    }
  }

// Compute the 4-particle correlations for the cross-check if needed.
  Int_t fourParticleCrossCheckArray[4] = {0}; // Array to save the harmonics with the format: cos(mphi1 + nphi2 - mphi3 - nphi4).
  TComplex fourParticleCrossCheckComplex = TComplex(0., 0.);  // Complex value for the 4-particle correlation.
  Double_t fourParticleCrossCheckValue = 0.;  // Real part of the 4-particle correlation.
  Int_t iBinFourCrossCheck = 1; // Index of the corresponding bin in the TProfile.
  Int_t iBinFourCrossCheckMiddle = 0.;  // Middle of the corresponding bin in the TProfile.
  for(Int_t m = 2; m < 7; m++)
  {
    for (Int_t n = 1; n < m; n++)
    {
    // Compute the correlation.
      fourParticleCrossCheckArray[0] = m;
      fourParticleCrossCheckArray[1] = n;
      fourParticleCrossCheckArray[2] = -1*m;
      fourParticleCrossCheckArray[3] = -1*n;
      fourParticleCrossCheckComplex = (CalculateRecursion(4, fourParticleCrossCheckArray))/fourParticleDenominator;
      fourParticleCrossCheckValue = fourParticleCrossCheckComplex.Re();

    // Fill the corresponding bin in the TProfile.
      iBinFourCrossCheckMiddle = (1.*iBinFourCrossCheck) - 0.5;
      fProfileFourParticleCorrelationsCrossCheck->Fill(iBinFourCrossCheckMiddle, fourParticleCrossCheckValue, fourParticleDenominator);
      fProfileFourParticleCorrelationsCrossCheck->GetXaxis()->SetBinLabel(iBinFourCrossCheck, Form("(%d,%d)", m, n));

    // Reset of the variables for the next pair of harmonics.
     iBinFourCrossCheck++;
     fourParticleCrossCheckComplex = TComplex(0., 0.);
     fourParticleCrossCheckValue = 0.;
     iBinFourCrossCheckMiddle = 0.;
    }
  }

// Compute the 6-particle correlations.
  Int_t sixZerosArray[6] = {0}; // Zero array for the denominator.
  Double_t sixParticleDenominator = CalculateRecursion(6, sixZerosArray).Re();  // 6-particle event weight.

  Int_t sixParticleArray[6] = {0};  // Array to save the harmonics with the format: cos(lphi1 + mphi2 + nphi3 - lphi4 - mphi5 - nphi6).
  TComplex sixParticleComplex = TComplex(0., 0.); // Complex value for the 6-particle correlation.
  Double_t sixParticleValue = 0.; // Real part of the 6-particle correlation.
  Int_t iBinSix = 1;  // Index of the corresponding bin in the TProfile.
  Int_t iBinSixMiddle = 0.; // Middle of the corresponding bin in the TProfile.
  for (Int_t l = 1; l < 5; l++)
  {
    for (Int_t m = 2; m < 6; m++)
    {
      if (l >= m) {continue;}
      for (Int_t n = 3; n < 7; n++)
      {
        if ((l >= n) || (m >= n)) {continue;}

      // Compute the correlation.
        sixParticleArray[0] = l;
        sixParticleArray[1] = m;
        sixParticleArray[2] = n;
        sixParticleArray[3] = -1*l;
        sixParticleArray[4] = -1*m;
        sixParticleArray[5] = -1*n;
        sixParticleComplex = (CalculateRecursion(6, sixParticleArray))/sixParticleDenominator;
        sixParticleValue = sixParticleComplex.Re();

      // Fill the corresponding bin in the TProfile.
        iBinSixMiddle = (1.*iBinSix) - 0.5;
        fProfileSixParticleCorrelations->Fill(iBinSixMiddle, sixParticleValue, sixParticleDenominator);
        fProfileSixParticleCorrelations->GetXaxis()->SetBinLabel(iBinSix, Form("(%d,%d,%d)", l, m, n));

      // Reset of the variables for the next pair of harmonics.
        iBinSix++;
        sixParticleComplex = TComplex(0., 0.);
        sixParticleValue = 0.;
        iBinSixMiddle = 0.; 
      }
    }
  }

// Reset the variables for the next event.
  twoParticleDenominator = 0.;
  fourParticleDenominator = 0.;
  sixParticleDenominator = 0.;
  iBinFour = 1;
  iBinFourCrossCheck =1;
  iBinSix = 1;
}

//======================================================================================//
TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
/* Calculate the multi-particle correlators by using the recursion method (an improved faster version) originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). */
  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= CalculateRecursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(CalculateRecursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += CalculateRecursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::ComputeTwoNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin)
{
/* Compute the 2-particle correlations using nested loops. */
  Double_t twoParticleCosine = 0.;  // cos(kphi_1 + lphi_2)).
  Double_t twoParticleWeight = 1.;  // Total particle weights.

  for (long long m = 0; m < nParticles; m++)
  {
    for (long long n = 0; n < nParticles; n++)
    {
      if (m == n) {continue;} // Remove the autocorrelations.

      twoParticleCosine = TMath::Cos(harmonic[0]*aAngles[m] + harmonic[1]*aAngles[n]);
      twoParticleWeight = weights[m] + weights[n];
      profile->Fill(middleBin, twoParticleCosine, twoParticleWeight);

    // Reset the variables.
      twoParticleCosine = 0.;
      twoParticleWeight = 1.;
    }
  }
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::ComputeFourNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin)
{
/* Compute the 4-particle correlations using nested loops. */
  Double_t fourParticleCosine = 0.; // cos(kphi_1 +l phi_2 + mphi_3 + nphi_4).
  Double_t fourParticleWeight = 1.; // Total particle weights.

  for (long long k = 0; k < nParticles; k++)
  {
    for (long long l = 0; l < nParticles; l++)
    {
      if (k == l) {continue;} // Remove the autocorrelations.
      for (long long m = 0; m < nParticles; m++)
      {
        if ((k == m) || (l == m)) {continue;}
        for (long long n = 0; n < nParticles; n++)
        {
          if ((k == n) || (l == n) || (m == n)) {continue;}

          fourParticleCosine = TMath::Cos(harmonic[0]*aAngles[k] + harmonic[1]*aAngles[l] + harmonic[2]*aAngles[m] + harmonic[3]*aAngles[n]);
          fourParticleWeight = weights[k] + weights[l] + weights[m] + weights[n];
          profile->Fill(middleBin, fourParticleCosine, fourParticleWeight);

        // Reset the variables.
          fourParticleCosine = 0.;
          fourParticleWeight = 1.;
        }
      }
    }
  }
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::ComputeTwoParticleEtaGaps(long long nParticles, Double_t angles[], Double_t pWeights[], Double_t pseudorapidity[])
{
/* Compute the 2-particle correlations with eta gaps. */
  TComplex Qminus[6][11] = {{TComplex(0., 0.)}};  // Q-vectors for the negative subset of the eta range, for v_1 to v_6.
  TComplex Qplus[6][11] = {{TComplex(0., 0.)}}; // Q-vectors for the positive subset of the eta range, for v_1 to v_6.
  Double_t Mminus[6][11] = {{0.}};  // Number of particles in the negative subset of the eta range.
  Double_t Mplus[6][11] = {{0.}}; // Number of particle in the positive subset of the eta range.
  Double_t EtaGaps[11] = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.}; // Choice of values of the eta gap.
  Double_t iAngle = 0.; // Azimuthal angle of the current particle.
  Double_t iWeight = 0.;  // Particle weight of the current particle.
  Double_t iEta = 0.; // Pseudorapidity of the current particle.
  Double_t iWeightToPowerP = 1.;  // Particle weight rised to the power p. TBA: use of non-unit particle weights.
  TComplex twoParticleComplex = TComplex(0., 0.); // Complex value of the 2-p correlations.
  Double_t twoParticleValue = 0.; // Real value of the 2-p correlations.

// Compute the Q-vectors for the negative and positive subsets of the eta range.
  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
    iAngle = angles[iParticle];
    iWeight = pWeights[iParticle];
    iEta = pseudorapidity[iParticle];

    if (iEta < 0.)  // Negative subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 6; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta < ((-0.5)*EtaGaps[iGap]))
          {
            Qminus[iHarmo][iGap] += TComplex(iWeightToPowerP*TMath::Cos((iHarmo+1)*iAngle), iWeightToPowerP*TMath::Sin((iHarmo+1)*iAngle));
            Mminus[iHarmo][iGap] += iWeightToPowerP;
          }
        }
      }
    }
    else if (iEta > 0.) // Positive subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 6; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta > (0.5*EtaGaps[iGap]))
          {
            Qplus[iHarmo][iGap] += TComplex(iWeightToPowerP*TMath::Cos((iHarmo+1)*iAngle), iWeightToPowerP*TMath::Sin((iHarmo+1)*iAngle));
            Mplus[iHarmo][iGap] += iWeightToPowerP;
          }
        }
      }
    }
    else {continue;}  // Particle with iEta = 0.
  }

// Compute the 2-p correlations using Qminus and Qplus.
  for (Int_t iHarmo = 0; iHarmo < 6; iHarmo++)
  {
    for (Int_t iGap = 0; iGap < 11; iGap++)
    {
      if (!( (Qminus[iHarmo][iGap].TComplex::Rho() > 0.) && (Qplus[iHarmo][iGap].TComplex::Rho() > 0.) )) {continue;}
      if (!( (Mminus[iHarmo][iGap] > 0.) && (Mplus[iHarmo][iGap] > 0.) )) {continue;}

      twoParticleComplex = Qminus[iHarmo][iGap]*TComplex::Conjugate(Qplus[iHarmo][iGap]);
      twoParticleValue = (twoParticleComplex.Re())/(Mminus[iHarmo][iGap]*Mplus[iHarmo][iGap]);
      fProfileTwoParticleCorrelationsWithEtaGaps[iHarmo]->Fill(iGap + 0.5, twoParticleValue, Mminus[iHarmo][iGap]*Mplus[iHarmo][iGap]);

    // Reset of the value of the 2-p correlations.
      twoParticleComplex = TComplex(0.,0.);
      twoParticleValue = 0.;
    }
  }

  iWeight++;  // Dummy operation to use iWeight.

}

//######################################################################################//
// Methods called in 'Terminate'.
//======================================================================================//


