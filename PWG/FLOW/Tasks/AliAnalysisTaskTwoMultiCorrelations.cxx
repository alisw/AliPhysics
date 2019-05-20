
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
  fHistoCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoIntermediateNumberOfTracks(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHistoFilterCorrelations(NULL),
  fHistoFinalFilterCorrelations(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
  fNumberOfBinsHNOT(30000),
  fMaxBinHNOT(30000.),
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
  fHistoIntermediatePt(NULL),
  fHistoFinalPt(NULL),
  fHistoIntermediateEta(NULL),
  fHistoFinalEta(NULL),
  fHistoIntermediatePhi(NULL),
  fHistoFinalPhi(NULL),
  fHistoIntermediateNumberOfTPC(NULL),
  fHistoFinalNumberOfTPC(NULL),
  fHistoIntermediateChiSquare(NULL),
  fHistoFinalChiSquare(NULL),
  fHistoIntermediateDCAxy(NULL),
  fHistoFinalDCAxy(NULL),
  fHistoIntermediateDCAz(NULL),
  fHistoFinalDCAz(NULL),
  fHistoIntermediateCharge(NULL),
  fHistoFinalCharge(NULL),
  fHistoIntermediateNumberOfITS(NULL),
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
  fHistoCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoIntermediateNumberOfTracks(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHistoFilterCorrelations(NULL),
  fHistoFinalFilterCorrelations(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
  fNumberOfBinsHNOT(30000),
  fMaxBinHNOT(30000.),
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
  fHistoIntermediatePt(NULL),
  fHistoFinalPt(NULL),
  fHistoIntermediateEta(NULL),
  fHistoFinalEta(NULL),
  fHistoIntermediatePhi(NULL),
  fHistoFinalPhi(NULL),
  fHistoIntermediateNumberOfTPC(NULL),
  fHistoFinalNumberOfTPC(NULL),
  fHistoIntermediateChiSquare(NULL),
  fHistoFinalChiSquare(NULL),
  fHistoIntermediateDCAxy(NULL),
  fHistoFinalDCAxy(NULL),
  fHistoIntermediateDCAz(NULL),
  fHistoFinalDCAz(NULL),
  fHistoIntermediateCharge(NULL),
  fHistoFinalCharge(NULL),
  fHistoIntermediateNumberOfITS(NULL),
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

// Create the mother list.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE); // Gives ownership of the elements inside the TList to the TList itself.

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
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the analysis for each event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";

  AliAODEvent *currentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());  // Pointer to an AOD event.
  AliMCEvent *currentMCEvent = MCEvent(); // Pointer to a MC event.

// Select the analysis method according to the type of files.
  if ((Int_t)fProcessOnlyAOD + (Int_t)fProcessOnlyMC + (Int_t)fProcessBothMCandAOD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one fProcess must be set to kTRUE.");
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
/* Book all the lists in the output file */
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

// Daughter list with the 2-p correlations
  fTwoParticleCorrelationsWithEtaGapsList = new TList();
  fTwoParticleCorrelationsWithEtaGapsList->SetName("fTwoParticleCorrelationsWithEtaGapsList");
  fTwoParticleCorrelationsWithEtaGapsList->SetOwner(kTRUE);
  fMainList->Add(fTwoParticleCorrelationsWithEtaGapsList);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookMultiplicityList()
{
/* Book the histograms containing the number of tracks. */
// Distribution of the centrality of the events.
  fHistoCentrality = new TH1D("fHistoCentrality", "Distribution of the centrality", 100, 0., 100.);
  fHistoCentrality->SetStats(kTRUE);
  fHistoCentrality->GetXaxis()->SetTitle("Centrality percentile");
  fHistoCentrality->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoCentrality);

// Distribution of the initial number of tracks.
  fHistoInitialNumberOfTracks = new TH1I("fHistoInitialNumberOfTracks", "Distribution of the initial number of tracks", fNumberOfBinsHNOT, 0., fMaxBinHNOT);
  fHistoInitialNumberOfTracks->SetStats(kTRUE);
  fHistoInitialNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fHistoInitialNumberOfTracks->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoInitialNumberOfTracks);

// Distribution of the number of tracks before the track selection.
  fHistoIntermediateNumberOfTracks = new TH1I("fHistoIntermediateNumberOfTracks", "Distribution of the number of tracks before the track selection", fNumberOfBinsHNOT, 0., fMaxBinHNOT);
  fHistoIntermediateNumberOfTracks->SetStats(kTRUE);
  fHistoIntermediateNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fHistoIntermediateNumberOfTracks->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoIntermediateNumberOfTracks);

// Distribution of the final number of tracks.
  fHistoFinalNumberOfTracks = new TH1I("fHistoFinalNumberOfTracks", "Distribution of the final number of tracks", fNumberOfBinsHNOT, 0., fMaxBinHNOT);
  fHistoFinalNumberOfTracks->SetStats(kTRUE);
  fHistoFinalNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fHistoFinalNumberOfTracks->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoFinalNumberOfTracks);

// 2D correlation histograms between two filters.
  if (fDoTDCorrelationHisto) // Manually enabled only to decide the cuts, as it is memory-greedy.
  {
    fHistoFilterCorrelations = new TH2I("fHistoFilterCorrelations", "2D multiplicity correlations between filters", fNumberOfBinsHFC, 0., fMaxBinHFC, fNumberOfBinsHFC, 0., fMaxBinHFC);
    fHistoFilterCorrelations->SetStats(kTRUE);
    fHistoFilterCorrelations->GetXaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fGlobalFilter));
    fHistoFilterCorrelations->GetYaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fMainFilter));
    fHistoFilterCorrelations->SetMarkerSize(0.5);
    fHistoFilterCorrelations->SetMarkerColor(kBlue);
    fHistoFilterCorrelations->SetMarkerStyle(kFullCircle);
    fMultiplicityList->Add(fHistoFilterCorrelations);

    fHistoFinalFilterCorrelations = new TH2I("fHistoFinalFilterCorrelations", "2D multiplicity correlations between filters after the track selection", fNumberOfBinsHFC, 0., fMaxBinHFC, fNumberOfBinsHFC, 0., fMaxBinHFC);
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
// Initial distribution of the x-position of the PV.
  fHistoInitialPVX = new TH1D("fHistoInitialPVX", "Initial distribution of PV_{x}", 1000, -20., 20.);
  fHistoInitialPVX->SetStats(kTRUE);
  fHistoInitialPVX->GetXaxis()->SetTitle("PV_{x}");
  fHistoInitialPVX->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVX);

// Final distribution of the x-position of the PV.
  fHistoFinalPVX = new TH1D("fHistoFinalPVX", "Final distribution of PV_{x}", 1000, -20., 20.);
  fHistoFinalPVX->SetStats(kTRUE);
  fHistoFinalPVX->GetXaxis()->SetTitle("PV_{x}");
  fHistoFinalPVX->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVX);

// Initial distribution of the y-position of the PV.
  fHistoInitialPVY = new TH1D("fHistoInitialPVY", "Initial distribution of PV_{y}", 1000, -20., 20.);
  fHistoInitialPVY->SetStats(kTRUE);
  fHistoInitialPVY->GetXaxis()->SetTitle("PV_{y}");
  fHistoInitialPVY->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVY);

// Final distribution of the y-position of the PV.
  fHistoFinalPVY = new TH1D("fHistoFinalPVY", "Final distribution of PV_{y}", 1000, -20., 20.);
  fHistoFinalPVY->SetStats(kTRUE);
  fHistoFinalPVY->GetXaxis()->SetTitle("PV_{y}");
  fHistoFinalPVY->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVY);

// Initial distribution of the z-position of the PV.
  fHistoInitialPVZ = new TH1D("fHistoInitialPVZ", "Initial distribution of PV_{z}", 1000, -20., 20.);
  fHistoInitialPVZ->SetStats(kTRUE);
  fHistoInitialPVZ->GetXaxis()->SetTitle("PV_{z}");
  fHistoInitialPVZ->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoInitialPVZ);

// Final distribution of the z-position of the PV.
  fHistoFinalPVZ = new TH1D("fHistoFinalPVZ", "Final distribution of PV_{z}", 1000, -20., 20.);
  fHistoFinalPVZ->SetStats(kTRUE);
  fHistoFinalPVZ->GetXaxis()->SetTitle("PV_{z}");
  fHistoFinalPVZ->GetYaxis()->SetTitle("Number of events");
  fEventSelectionList->Add(fHistoFinalPVZ);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookTrackSelectionList()
{
/* Book the control histograms for the track selection criteria. */
// Distributions of the transverse momentum.
  fHistoIntermediatePt = new TH1D("fHistoIntermediatePt", "Distribution of p_{T} before the track selection", 1000, 0., 20.);
  fHistoIntermediatePt->SetStats(kTRUE);
  fHistoIntermediatePt->GetXaxis()->SetTitle("p_{T}");
  fHistoIntermediatePt->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediatePt);

  fHistoFinalPt = new TH1D("fHistoFinalPt", "Distribution of p_{T} before the track selection", 1000, 0., 20.);
  fHistoFinalPt->SetStats(kTRUE);
  fHistoFinalPt->GetXaxis()->SetTitle("p_{T}");
  fHistoFinalPt->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalPt);

// Distributions of the pseudorapidity.
  fHistoIntermediateEta = new TH1D("fHistoIntermediateEta", "Distribution of #eta before the track selection", 1000, -5.5, 5.5);
  fHistoIntermediateEta->SetStats(kTRUE);
  fHistoIntermediateEta->GetXaxis()->SetTitle("#eta");
  fHistoIntermediateEta->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateEta);

  fHistoFinalEta = new TH1D("fHistoFinalEta", "Distribution of #eta after the track selection", 1000, -5.5, 5.5);
  fHistoFinalEta->SetStats(kTRUE);
  fHistoFinalEta->GetXaxis()->SetTitle("#eta");
  fHistoFinalEta->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalEta);

// Distributions of the azimuthal angles.
  fHistoIntermediatePhi = new TH1D("fHistoIntermediatePhi", "Distribution of #phi before the track selection", 1000, 0., 6.3);
  fHistoIntermediatePhi->SetStats(kTRUE);
  fHistoIntermediatePhi->GetXaxis()->SetTitle("#phi");
  fHistoIntermediatePhi->GetXaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediatePhi);

  fHistoFinalPhi = new TH1D("fHistoFinalPhi", "Distribution of #phi after the track selection", 1000, 0., 6.3);
  fHistoFinalPhi->SetStats(kTRUE);
  fHistoFinalPhi->GetXaxis()->SetTitle("#phi");
  fHistoFinalPhi->GetXaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalPhi);

// Distributions of the number of TPC clusters before the track selection.
  fHistoIntermediateNumberOfTPC = new TH1I("fHistoIntermediateNumberOfTPC", "Distribution of the number of TPC clusters before the track selection", 1000, 0., 170.);
  fHistoIntermediateNumberOfTPC->SetStats(kTRUE);
  fHistoIntermediateNumberOfTPC->GetXaxis()->SetTitle("Number of TPC clusters");
  fHistoIntermediateNumberOfTPC->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateNumberOfTPC);

  fHistoFinalNumberOfTPC = new TH1I("fHistoFinalNumberOfTPC", "Distribution of the number of TPC clusters after the track selection", 1000, 0., 170.);
  fHistoFinalNumberOfTPC->SetStats(kTRUE);
  fHistoFinalNumberOfTPC->GetXaxis()->SetTitle("Number of TPC clusters");
  fHistoFinalNumberOfTPC->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalNumberOfTPC);

// Distributions of the chi^2 of the track momentum in TPC.
  fHistoIntermediateChiSquare = new TH1D("fHistoIntermediateChiSquare", "Distribution of the #chi^{2} in the TPC before the track selection", 1000, 0., 20.);
  fHistoIntermediateChiSquare->SetStats(kTRUE);
  fHistoIntermediateChiSquare->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoIntermediateChiSquare->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateChiSquare);

  fHistoFinalChiSquare = new TH1D("fHistoFinalChiSquare", "Distribution of the #chi^{2} in the TPC after the track selection", 1000, 0., 20.);
  fHistoFinalChiSquare->SetStats(kTRUE);
  fHistoFinalChiSquare->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoFinalChiSquare->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalChiSquare);

// Distributions of the xy-coordinate of the DCA.
  fHistoIntermediateDCAxy = new TH1D("fHistoIntermediateDCAxy", "Distribution of DCA_{xy} before the track selection", 1000, 0., 10.);
  fHistoIntermediateDCAxy->SetStats(kTRUE);
  fHistoIntermediateDCAxy->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoIntermediateDCAxy->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateDCAxy);

  fHistoFinalDCAxy = new TH1D("fHistoFinalDCAxy", "Distribution of DCA_{xy} after the track selection", 1000, 0., 10.);
  fHistoFinalDCAxy->SetStats(kTRUE);
  fHistoFinalDCAxy->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fHistoFinalDCAxy->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalDCAxy);

// Distributions of the z-coordinate of the DCA.
  fHistoIntermediateDCAz = new TH1D("fHistoIntermediateDCAz", "Distribution of DCA_{z} before the track selection", 1000, 0., 10.);
  fHistoIntermediateDCAz->SetStats(kTRUE);
  fHistoIntermediateDCAz->GetXaxis()->SetTitle("DCA_{z}");
  fHistoIntermediateDCAz->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateDCAz);

  fHistoFinalDCAz = new TH1D("fHistoFinalDCAz", "Distribution of DCA_{z} after the track selection", 1000, 0., 10.);
  fHistoFinalDCAz->SetStats(kTRUE);
  fHistoFinalDCAz->GetXaxis()->SetTitle("DCA_{z}");
  fHistoFinalDCAz->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalDCAz);

// Distributions of the electric charge of the tracks.
  fHistoIntermediateCharge = new TH1I("fHistoIntermediateCharge", "Distribution of the electric charge before the track selection", 2, -2, 2);
  fHistoIntermediateCharge->SetStats(kTRUE);
  fHistoIntermediateCharge->GetXaxis()->SetTitle("Charge");
  fHistoIntermediateCharge->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateCharge);

  fHistoFinalCharge = new TH1I("fHistoFinalCharge", "Distribution of the electric charge after the track selection", 2, -2, 2);
  fHistoFinalCharge->SetStats(kTRUE);
  fHistoFinalCharge->GetXaxis()->SetTitle("Charge");
  fHistoFinalCharge->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoFinalCharge);

// Distributions of the number of clusters in the ITS.
  fHistoIntermediateNumberOfITS = new TH1I("fHistoIntermediateNumberOfITS", "Distribution of the number of ITS clusters before the track selection", 1000, 0., 170.);
  fHistoIntermediateNumberOfITS->SetStats(kTRUE);
  fHistoIntermediateNumberOfITS->GetXaxis()->SetTitle("Number of ITS clusters");
  fHistoIntermediateNumberOfITS->GetYaxis()->SetTitle("Number of tracks");
  fTrackSelectionList->Add(fHistoIntermediateNumberOfITS);

  fHistoFinalNumberOfITS = new TH1I("fHistoFinalNumberOfITS", "Distribution of the number of ITS clusters after the track selection", 1000, 0., 170.);
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
  fProfileTwoParticleCorrelations = new TProfile("fProfileTwoParticleCorrelations", "2-particle correlations", 6, 0., 6.);
  fProfileTwoParticleCorrelations->SetStats(kTRUE);
  fProfileTwoParticleCorrelations->Sumw2();
  fProfileTwoParticleCorrelations->GetXaxis()->SetTitle("n");
  fProfileTwoParticleCorrelations->GetYaxis()->SetTitle("#LT#LT2#GT#GT_{n,-n}");
  fMultiParticleCorrelationsList->Add(fProfileTwoParticleCorrelations);

  if (fCrosscheckWithNestedLoops)
  {
    fProfileTwoParticleCorrelationsNestedLoops = new TProfile("fProfileTwoParticleCorrelationsNestedLoops", "2-particle correlations with nested loops", 6, 0., 6.);
    fProfileTwoParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileTwoParticleCorrelationsNestedLoops->Sumw2();
    fProfileTwoParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("n");
    fProfileTwoParticleCorrelationsNestedLoops->GetYaxis()->SetTitle("#LT#LT2#GT#GT_{n,-n}");
    fMultiParticleCorrelationsList->Add(fProfileTwoParticleCorrelationsNestedLoops);
  }

// 4-particle correlations.
  fProfileFourParticleCorrelations = new TProfile("fProfileFourParticleCorrelations", "4-particle correlations", 21, 0., 21.);
  fProfileFourParticleCorrelations->SetStats(kTRUE);
  fProfileFourParticleCorrelations->Sumw2();
  fProfileFourParticleCorrelations->GetXaxis()->SetTitle("(m,n)");
  fProfileFourParticleCorrelations->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
  fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelations);

  fProfileFourParticleCorrelationsCrossCheck = new TProfile("fProfileFourParticleCorrelationsCrossCheck", "4-particle correlations for cross-check", 15, 0., 15.);
  fProfileFourParticleCorrelationsCrossCheck->SetStats(kTRUE);
  fProfileFourParticleCorrelationsCrossCheck->Sumw2();
  fProfileFourParticleCorrelationsCrossCheck->GetXaxis()->SetTitle("(m,n)");
  fProfileFourParticleCorrelationsCrossCheck->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
  fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelationsCrossCheck);

  if (fCrosscheckWithNestedLoops)
  {
    fProfileFourParticleCorrelationsNestedLoops = new TProfile("fProfileFourParticleCorrelationsNestedLoops", "4-particle correlations with nested loops", 21, 0., 21.);
    fProfileFourParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileFourParticleCorrelationsNestedLoops->Sumw2();
    fProfileFourParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("(m,n)");
    fProfileFourParticleCorrelationsNestedLoops->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
    fMultiParticleCorrelationsList->Add(fProfileFourParticleCorrelationsNestedLoops);
  }

// 6-particle correlations.
  fProfileSixParticleCorrelations = new TProfile("fProfileSixParticleCorrelations", "6-particle correlations", 20, 0., 20.);
  fProfileSixParticleCorrelations->SetStats(kTRUE);
  fProfileSixParticleCorrelations->Sumw2();
  fProfileSixParticleCorrelations->GetXaxis()->SetTitle("(l,m,n)");
  fProfileSixParticleCorrelations->GetYaxis()->SetTitle("#LT#LT6#GT#GT_{l,m,n,-l,-m,-n}");
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
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->SetStats(kTRUE);
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->Sumw2();
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->GetXaxis()->SetTitle("#eta gap");
      fProfileTwoParticleCorrelationsWithEtaGaps[i]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,-%d}", i+1, i+1));
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
void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent(AliMCEvent *aMCevent)
{
/* Execute the analysis for the provided MC event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent(AliMCEvent *aMCevent)";

// Check if there is an event or not.
  if (!aMCevent) {Fatal(sMethodName.Data(), "ERROR: no MC event found.");}

// Select the detector to use for the estimation of the centrality.
  TString centralityEstimator = "centralityEstimator";  // Name of the detector used for the centrality estimation.
  if ((Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be selected in 'SetCentralityEstimation'.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Determine if the event belongs to this centrality range (for reconstructed particles only).
  Double_t aCentrality = 0.;  // Centrality of the given event.
  if (!fProcessOnlyMC)
  {
    AliMultSelection *ams = (AliMultSelection*)aMCevent->FindListObject("MultSelection");
    if (!ams) {return;} // Protection against NULL pointer.
    aCentrality = ams->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));
    if ((aCentrality >= fCentralityMin) && (aCentrality < fCentralityMax))
    {
      fHistoCentrality->Fill(aCentrality);
    }
    else {return;}  // This event does not belong to this centrality range.
  }

// Get the number of tracks before the event selection.
  long long initialNumberOfTracks = aMCevent->GetNumberOfTracks();
  fHistoInitialNumberOfTracks->Fill(initialNumberOfTracks);

// Cuts on the position of the Primary Vertex.
  AliMCVertex *avtx = (AliMCVertex*)aMCevent->GetPrimaryVertex(); // 3d position of the PV.
  fHistoInitialPVX->Fill(avtx->GetX());
  fHistoInitialPVY->Fill(avtx->GetY());
  fHistoInitialPVZ->Fill(avtx->GetZ());

  if (fCutOnPVX) {if ((avtx->GetX() < fPVXMin) || (avtx->GetX() > fPVXMax)) {return;}}
  if (fCutOnPVY) {if ((avtx->GetY() < fPVYMin) || (avtx->GetY() > fPVYMax)) {return;}}
  if (fCutOnPVZ) {if ((avtx->GetZ() < fPVZMin) || (avtx->GetZ() > fPVZMax)) {return;}}

/// TBA: more event cuts?

// Preparations for the track selection.
  long long numberOfTracksBeforeTrackSelection = aMCevent->GetNumberOfTracks();  // Number of tracks before the track selection.
  fHistoIntermediateNumberOfTracks->Fill(numberOfTracksBeforeTrackSelection);
  long long finalNumberOfTracks = 0;  // Number of tracks after the full selection.
  Int_t *IsTrackSelected = new Int_t[numberOfTracksBeforeTrackSelection](); // Flag to indicate a track passed the track selection (1) or not (0).

  Double_t pT = 0.; // Transverse momentum.
  Double_t eta = 0.;  // Pseudorapidity.
  Double_t phi = 0.;  // Azimuthal angle.
  Int_t charge = 0; // Electric charge.

// Look at each track in the event to mark them as selected or not.
  for (long long iTrack = 0; iTrack < numberOfTracksBeforeTrackSelection; iTrack++)
  {
    AliAODMCParticle *currentTrack = dynamic_cast<AliAODMCParticle*>(aMCevent->GetTrack(iTrack));  // Pointer to the MC track.
    if (!currentTrack) {continue;}  // Protection against NULL pointer.

  // Get all the observables for the track selection.
    pT = currentTrack->Pt();
    eta = currentTrack->Eta();
    phi = currentTrack->Phi();
    charge = currentTrack->Charge();

  // Fill the histograms before the track selection.
    fHistoIntermediatePt->Fill(pT);
    fHistoIntermediateEta->Fill(eta);
    fHistoIntermediatePhi->Fill(phi);
    fHistoIntermediateCharge->Fill(charge);

  // Apply the track selection to the provided track.
    Bool_t cutOnCharge = kTRUE; // Set to kTRUE by default in case no selection over the charge is done.
    if (fCutOnCharge) // Check if the track passes the cut in case it is applied.
    {
      if (charge != fCharge) {cutOnCharge = kFALSE;}
    }

    if ((fPtMin <= pT) && (pT <= fPtMax) && (fEtaMin <= eta) && (eta <= fEtaMax) && (cutOnCharge))  // Apply the cuts to the track.
    {
      IsTrackSelected[iTrack] = 1;
      finalNumberOfTracks++;
    }
    else {IsTrackSelected[iTrack] = 0;}  // The track failed the selection.
  }

// Remove the events with too few or too many tracks.
  Int_t cutValueMaxNumberOfTracks = 0;  // Value of the cut on the maximum number of tracks.
  if (finalNumberOfTracks <= fMultiplicityMin) {return;}
  if (fCutOnTracksMax)  // If the cuts on the maximum numbers of tracks are enabled.
  {
  // Determine the value to cut depending on the centrality.
    if ((aCentrality >= 0.) && (aCentrality < 5.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxZero;}
    else if ((aCentrality >= 5.) && (aCentrality < 10.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFive;}
    else if ((aCentrality >= 10.) && (aCentrality < 20.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTen;}
    else if ((aCentrality >= 20.) && (aCentrality < 30.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTwenty;}
    else if ((aCentrality >= 30.) && (aCentrality < 40.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxThirty;}
    else if ((aCentrality >= 40.) && (aCentrality < 50.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxForty;}
    else if ((aCentrality >= 50.) && (aCentrality < 60.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFifty;}
    else if ((aCentrality >= 60.) && (aCentrality < 70.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSixty;}
    else if ((aCentrality >= 70.) && (aCentrality < 80.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSeventy;}

  // Apply the cut.
    if (finalNumberOfTracks >= cutValueMaxNumberOfTracks) {return;}
  }

// Fill all the event histograms after the full selection.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);
  fHistoFinalPVX->Fill(avtx->GetX());
  fHistoFinalPVY->Fill(avtx->GetY());
  fHistoFinalPVZ->Fill(avtx->GetZ());

// Define the arrays for the azimuthal angles and particle weights to use in the analysis.
  Double_t *etaArray = new Double_t[finalNumberOfTracks](); // Pseudorapidity.
  Double_t *phiArray = new Double_t[finalNumberOfTracks](); // Azimuthal angles.
  Double_t *particleWeightArray = new Double_t[finalNumberOfTracks](); // Particle weights.
  Int_t indexInNewArrays = 0; // New index of the track if it passed the selection.

// Loop over the tracks to keep only the selected ones.
  for (long long iTrack = 0; iTrack < numberOfTracksBeforeTrackSelection; iTrack++)
  {
    AliAODMCParticle *aTrack = dynamic_cast<AliAODMCParticle*>(aMCevent->GetTrack(iTrack)); // Pointer to the MC track.
    if (!aTrack) {continue;}  // Protection against NULL pointer.

    if (IsTrackSelected[iTrack] == 1) // The particle passed the selection.
    {
    // Get all the observables used in the track selection.
      pT = aTrack->Pt();
      etaArray[indexInNewArrays] = aTrack->Eta();
      phiArray[indexInNewArrays] = aTrack->Phi();
      charge = aTrack->Charge();

      if (fUseParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA.");}
      else {particleWeightArray[indexInNewArrays] = 1.;}

    // Fill all the track histograms after the full selection.
      fHistoFinalPt->Fill(pT);
      fHistoFinalEta->Fill(etaArray[indexInNewArrays]);
      fHistoFinalPhi->Fill(phiArray[indexInNewArrays]);
      fHistoFinalCharge->Fill(charge);

    // Increase the value of 'indexInNewArrays' by one.
      indexInNewArrays++;
    }
    else {continue;}
  }

// Calculate the Q-vectors for the current event.
  CalculateQvectors(finalNumberOfTracks, phiArray, particleWeightArray);

// Compute all the multiparticle correlations for the current event.
  ComputeMultiparticleCorrelations(finalNumberOfTracks, phiArray, particleWeightArray);

// Calculate the 2-particle correlations with eta gaps if selected.
  if (fComputeEtaGaps) {ComputeTwoParticleEtaGaps(finalNumberOfTracks, phiArray, particleWeightArray, etaArray);}

// Reset everything to zero for the next event.
  aCentrality = 0.;
  numberOfTracksBeforeTrackSelection = 0;
  cutValueMaxNumberOfTracks = 0;
  finalNumberOfTracks = 0;
  delete [] IsTrackSelected;
  pT = 0.;
  eta = 0.;
  phi = 0.;
  charge = 0;
  delete [] etaArray;
  delete [] phiArray;
  delete [] particleWeightArray;
  indexInNewArrays = 0;
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *aAODevent)
{
/* Execute the analysis for the provided AOD event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *aAODevent)";
  if (!aAODevent) {Fatal(sMethodName.Data(), "ERROR: no AOD event found.");}

// Select the detector used to determine the centrality.
  TString centralityEstimator = "centralityEstimator";  // Name of the selected detector.
  if ((Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Determine if the event belongs to the current centrality range.
  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if (!ams) {return;} // Protection against NULL pointer.

  Double_t aCentrality = ams->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));  // Centrality of the given event.
  if ((aCentrality >= fCentralityMin) && (aCentrality < fCentralityMax))
  {
    fHistoCentrality->Fill(aCentrality);
  }
  else {return;}  // This event does not belong to this centrality range.

// Get the initial number of tracks.
  long long initialNumberOfTracks = aAODevent->GetNumberOfTracks();
  fHistoInitialNumberOfTracks->Fill(initialNumberOfTracks);

// Apply the event selection criteria.
/// Cuts on the position of the PV.
  AliAODVertex *aVertex = (AliAODVertex*)aAODevent->GetPrimaryVertex(); // 3D position of the PV.
  fHistoInitialPVX->Fill(aVertex->GetX());
  fHistoInitialPVY->Fill(aVertex->GetY());
  fHistoInitialPVZ->Fill(aVertex->GetZ());

  if (fCutOnPVX) {if ((aVertex->GetX() < fPVXMin) || (aVertex->GetX() > fPVXMax)) {return;}}
  if (fCutOnPVY) {if ((aVertex->GetY() < fPVYMin) || (aVertex->GetY() > fPVYMax)) {return;}}
  if (fCutOnPVZ) {if ((aVertex->GetZ() < fPVZMin) || (aVertex->GetZ() > fPVZMax)) {return;}}

/// TBA: Add more event cuts?

// Get the number of tracks after the event selection.
  long long intermediateNumberOfTracks = aAODevent->GetNumberOfTracks();  // Number of tracks before the track selection.
  fHistoIntermediateNumberOfTracks->Fill(intermediateNumberOfTracks);

// Do the 2D correlation histogram?
  if (fDoTDCorrelationHisto)
  {
    Int_t multiplicityGlobalFilter = 0; // Number of tracks in the global filter.
    Int_t multiplicityMainFilter = 0; // Number of tracks in the main filter bit.
    for (Int_t inTrack = 0; inTrack < intermediateNumberOfTracks; inTrack++)
    {
      AliAODTrack *aliTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(inTrack));  // Pointer to an AOD track.
      if (!aliTrack) {continue;}  // Protection against NULL pointer.

      //if ((aliTrack->TestFilterBit(fMainFilter)) && (!aliTrack->TestFilterBit(fGlobalFilter))) {multiplicityMainFilter++;}
      if (aliTrack->TestFilterBit(fGlobalFilter)) {multiplicityGlobalFilter++;}
      if (aliTrack->TestFilterBit(fMainFilter)) {multiplicityMainFilter++;}
    }

  // Fill the histogram.
    fHistoFilterCorrelations->Fill(multiplicityGlobalFilter, multiplicityMainFilter);

  // Reset the variables.
    multiplicityMainFilter = 0;
    multiplicityGlobalFilter = 0;
  }

// Prepare the variables for the track selection.
  long long finalNumberOfTracks = 0;  // Final number of tracks for the main filter.
  long long globalNumberOfTracks = 0; // Final number of tracks for the global filter.
  Int_t *IsTrackSelected = new Int_t[intermediateNumberOfTracks](); // Flag to indicate if a track has passed the track selection (1) or not (0).
  Double_t pT = 0.; // Transverse momentum of the track.
  Double_t eta = 0.;  // Pseudorapidity of the track.
  Double_t phi = 0.;  // Azimuthal angle of the track.
  Int_t numberOfTPCClusters = 0;  // Number of TPC clusters gone through by the track.
  Int_t numberOfITSClusters = 0;  // Number of ITS clusters gone through by the track.
  Double_t chiSquareInTPC = 0.; // Chi square in the TPC.
  Double_t DCAx = 0.; // x-coordinate of the DCA.
  Double_t DCAy = 0.; // y-coordinate of the DCA.
  Double_t DCAz = 0.;  // z-coordinate of the DCA.
  Double_t DCAxy = 0.;  // xy-coordinate of the DCA.
  Int_t charge = 0; // Electric charge.

// Apply the track selection to each track and mark them according to if they passed or not.
  for (long long iTrack = 0; iTrack < intermediateNumberOfTracks; iTrack++)
  {
  // Get a pointer to the current track.
    AliAODTrack *currentTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack));
    if (!currentTrack) {continue;}  // Protection against NULL pointer.
    if ((!currentTrack->TestFilterBit(fMainFilter)) && (!currentTrack->TestFilterBit(fGlobalFilter))) {continue;} // The track belongs neither to the main nor the global filter.

  // Get all the observables for the current track.
    pT = currentTrack->Pt();
    eta = currentTrack->Eta();
    phi = currentTrack->Phi();
    numberOfTPCClusters = currentTrack->GetTPCNcls();
    numberOfITSClusters = currentTrack->GetITSNcls();
    chiSquareInTPC = (currentTrack->GetTPCchi2())/(currentTrack->GetNcls(1));
    DCAx = currentTrack->XAtDCA();
    DCAy = currentTrack->YAtDCA();
    DCAz = currentTrack->ZAtDCA();
    charge = currentTrack->Charge();

    DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));

  // Fill the histograms before the application of the track selection criteria.
    fHistoIntermediatePt->Fill(pT);
    fHistoIntermediateEta->Fill(eta);
    fHistoIntermediatePhi->Fill(phi);
    fHistoIntermediateNumberOfTPC->Fill(numberOfTPCClusters);
    fHistoIntermediateChiSquare->Fill(chiSquareInTPC);
    fHistoIntermediateDCAxy->Fill(DCAxy);
    fHistoIntermediateDCAz->Fill(DCAz);
    fHistoIntermediateCharge->Fill(charge);
    fHistoIntermediateNumberOfITS->Fill(numberOfITSClusters);

  // Apply the track selection criteria.
    if (ApplyTrackSelection(currentTrack))  // The track passed the track selection.
    {
      IsTrackSelected[iTrack] = 1;
      if (currentTrack->TestFilterBit(fMainFilter)) {finalNumberOfTracks++;}
      if (currentTrack->TestFilterBit(fGlobalFilter)) {globalNumberOfTracks++;}
    }
    else {IsTrackSelected[iTrack] = 0;} // The track failed the selection.
  }

// Do the final 2D correlation histograms for filters?
  if (fDoTDCorrelationHisto)
  {
    fHistoFinalFilterCorrelations->Fill(globalNumberOfTracks, finalNumberOfTracks);
  }

// Remove the events with too few or too many tracks.
  if (finalNumberOfTracks <= fMultiplicityMin) {return;}  // The event does not have enough tracks.
  if (fCutOnTDCorrelations) // Apply the multiplicity cut from the 2D correlation plot.
  {
    if (finalNumberOfTracks <= (fMultiplicityMinA*globalNumberOfTracks + fMultiplicityMinB)) {return;}  // The number of tracks is under the accepted band.
    if (finalNumberOfTracks >= (fMultiplicityMaxA*globalNumberOfTracks + fMultiplicityMaxB)) {return;}  // The number of tracks is above the accepted band.
  }

/// Application of the brute-force method to remove HM outliers.
  Int_t cutValueMaxNumberOfTracks = 0;  // Value of the cut on the maximum number of tracks.
  if (fCutOnTracksMax)  // If the cuts on the maximum numbers of tracks are enabled.
  {
  // Determine the value to cut depending on the centrality.
    if ((aCentrality >= 0.) && (aCentrality < 5.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxZero;}
    else if ((aCentrality >= 5.) && (aCentrality < 10.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFive;}
    else if ((aCentrality >= 10.) && (aCentrality < 20.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTen;}
    else if ((aCentrality >= 20.) && (aCentrality < 30.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxTwenty;}
    else if ((aCentrality >= 30.) && (aCentrality < 40.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxThirty;}
    else if ((aCentrality >= 40.) && (aCentrality < 50.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxForty;}
    else if ((aCentrality >= 50.) && (aCentrality < 60.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxFifty;}
    else if ((aCentrality >= 60.) && (aCentrality < 70.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSixty;}
    else if ((aCentrality >= 70.) && (aCentrality < 80.)) {cutValueMaxNumberOfTracks = fNumberOfTracksMaxSeventy;}

  // Apply the cut.
    if (finalNumberOfTracks >= cutValueMaxNumberOfTracks) {return;}
  }

// Fill all the histograms after the full selection.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);
  fHistoFinalPVX->Fill(aVertex->GetX());
  fHistoFinalPVY->Fill(aVertex->GetY());
  fHistoFinalPVZ->Fill(aVertex->GetZ());

// Define the arrays to use in the reminding part of the analysis.
  Double_t *etaArray = new Double_t[finalNumberOfTracks](); // Pseudorapidity.
  Double_t *phiArray = new Double_t[finalNumberOfTracks](); // Azimuthal angles.
  Double_t *particleWeightArray = new Double_t[finalNumberOfTracks]();  // Particle weights.
  Int_t indexInNewArrays = 0; // New index of the track if it passed the selection.

// Loop over the tracks to keep only the selected ones.
  for (long long iiTrack = 0; iiTrack < intermediateNumberOfTracks; iiTrack++)
  {
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iiTrack));  // Pointer to the AOD track.
    if (!aTrack) {continue;}  // Protection against NULL pointer.
    if (!aTrack->TestFilterBit(fMainFilter)) {continue;}

    if (IsTrackSelected[iiTrack] == 1)  // The particle passed the selection.
    {
    // Get all the observables used in the track selection.
      pT = aTrack->Pt();
      etaArray[indexInNewArrays] = aTrack->Eta();
      phiArray[indexInNewArrays] = aTrack->Phi();
      numberOfTPCClusters = aTrack->GetTPCNcls();
      numberOfITSClusters = aTrack->GetITSNcls();
      chiSquareInTPC = (aTrack->GetTPCchi2())/(aTrack->GetNcls(1));
      DCAx = aTrack->XAtDCA();
      DCAy = aTrack->YAtDCA();
      DCAz = aTrack->ZAtDCA();
      DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));
      charge = aTrack->Charge();

      if (fUseParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA use of non-unit particle weights");}
      else {particleWeightArray[indexInNewArrays] = 1.;}

    // Fill all the track histograms after the full selection.
      fHistoFinalPt->Fill(pT);
      fHistoFinalEta->Fill(etaArray[indexInNewArrays]);
      fHistoFinalPhi->Fill(phiArray[indexInNewArrays]);
      fHistoFinalNumberOfTPC->Fill(numberOfTPCClusters);
      fHistoFinalChiSquare->Fill(chiSquareInTPC);
      fHistoFinalDCAxy->Fill(DCAxy);
      fHistoFinalDCAz->Fill(DCAz);
      fHistoFinalCharge->Fill(charge);
      fHistoFinalNumberOfITS->Fill(numberOfITSClusters);

    // Increase the value of 'indexInNewArrays' by one.
      indexInNewArrays++;
    }
    else {continue;}
  }

/* All the selection criteria have been applied. Start of the analysis itself. */
// Calculate the Q-vectors for the current event.
  CalculateQvectors(finalNumberOfTracks, phiArray, particleWeightArray);

// Compute all the multiparticle correlations for the current event.
  ComputeMultiparticleCorrelations(finalNumberOfTracks, phiArray, particleWeightArray);

// Calculate the 2-particle correlations with eta gaps if selected.
  if (fComputeEtaGaps) {ComputeTwoParticleEtaGaps(finalNumberOfTracks, phiArray, particleWeightArray, etaArray);}

// Reset everything to zero for the next event.
  intermediateNumberOfTracks = 0;
  cutValueMaxNumberOfTracks = 0;
  finalNumberOfTracks = 0;
  delete [] IsTrackSelected;
  pT = 0.;
  eta = 0.;
  phi = 0.;
  numberOfTPCClusters = 0;
  numberOfITSClusters = 0;
  chiSquareInTPC = 0.;
  DCAx = 0.;
  DCAy = 0.;
  DCAz = 0.;
  charge = 0;
  delete [] etaArray;
  delete [] phiArray;
  delete [] particleWeightArray;
  indexInNewArrays = 0;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)
{
/* Apply the track selection criteria and return if the track passed it or not. */
  Bool_t testOfPt = kTRUE;  // Cut on the transverse momentum.
  Bool_t testOfEta = kTRUE; // Cut on the pseudorapidity.
  Bool_t testOfNumberOfTPC = kTRUE; // Cut on the number of TPC clusters.
  Bool_t testOfChiSquareTPC = kTRUE;  // Cut on the chi^2 of the momentum in the TPC.
  Bool_t testOfDCAxy = kTRUE; // Cut on the DCA of the track in the xy-plane.
  Bool_t testOfDCAz = kTRUE;  // Cut on the DCA of the track along the z-direction.
  Bool_t testOfCharge = kTRUE;  // Cut on the electric charge.
  Bool_t testOfNumberOfITS = kTRUE; // Cut on the number of ITS clusters.

// Get the observables for the given track.
    Double_t pT = aAODtrack->Pt();
    Double_t eta = aAODtrack->Eta();
    Int_t numberOfTPCClusters = aAODtrack->GetTPCNcls();
    Int_t numberOfITSClusters = aAODtrack->GetITSNcls();
    Double_t chiSquareInTPC = (aAODtrack->GetTPCchi2())/(aAODtrack->GetNcls(1));
    Double_t DCAx = aAODtrack->XAtDCA();
    Double_t DCAy = aAODtrack->YAtDCA();
    Double_t DCAz = aAODtrack->ZAtDCA();
    Int_t charge = aAODtrack->Charge();

    Double_t DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));

// Apply the selection to the track.
  if (fCutOnPt) {testOfPt = (fPtMin <= pT) && (pT <= fPtMax);}
  if (fCutOnEta) {testOfEta = (fEtaMin <= eta) && (eta <= fEtaMax);}
  if (fCutOnNumberOfTPC) {testOfNumberOfTPC = (fNumberOfTPCMin < numberOfTPCClusters);}
  if (fCutOnChiSquarePInTPC) {testOfChiSquareTPC = (fChiSquarePInTPCMin <= chiSquareInTPC) && (chiSquareInTPC <= fChiSquarePInTPCMax);}
  if (fCutOnDCAxy) {testOfDCAxy = (DCAxy < fDCAxyMax);}
  if (fCutOnDCAz) {testOfDCAz = (DCAz < fDCAzMax);}
  if (fCutOnCharge) {testOfCharge = (charge == fCharge);}
  if (fCutOnNumberOfITS) {testOfNumberOfITS = (fNumberOfITSMin < numberOfITSClusters);}

  return testOfPt && testOfEta && testOfNumberOfTPC && testOfChiSquareTPC && testOfDCAxy && testOfDCAz && testOfCharge && testOfNumberOfITS;
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


