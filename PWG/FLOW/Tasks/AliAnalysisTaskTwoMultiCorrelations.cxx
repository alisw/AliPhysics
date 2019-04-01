
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
// Structure of the output file.
  fMainList(NULL),
  fQAListBeforeSelection(NULL),
  fQAListAfterSelection(NULL),
  fListCorrelations(NULL),
// General parameters.
  fMaxNumberOfParticlesInCorrelations(8),
  fHighestFlowHarmonic(6),
  fUseParticleWeights(kFALSE),
  fCrossCheckFourParticleCorrelations(kFALSE),
  fCrossCheckWithNestedLoops(kFALSE),
// Type of files used in the analysis.
  fProcessOnlyAOD(kFALSE),
  fProcessOnlyMC(kFALSE),
  fProcessBothMCandAOD(kFALSE),
// Determination of the centrality.
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
// Event selection.
  fCutOnVertexX(kFALSE),
  fVertexMinX(-44.),
  fVertexMaxX(-44.),
  fCutOnVertexY(kFALSE),
  fVertexMinY(-44.),
  fVertexMaxY(-44.),
  fCutOnVertexZ(kFALSE),
  fVertexMinZ(-10.),
  fVertexMaxZ(10.),
  fNumberOfTracksMin(6),
// Track selection.
  fCutOnPt(kFALSE),
  fPtMin(0.2),
  fPtMax(5.),
  fCutOnEta(kFALSE),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fFilter(128),
  fCutOnNumberOfTPC(kFALSE),
  fNumberOfTPCMin(70),
  fCutOnChiSquarePInTPC(kFALSE),
  fChiSquarePInTPCMin(0.1),
  fChiSquarePInTPCMax(4.),
  fCutOnDCA(kFALSE),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
  fCutOnCharge(kFALSE),
  fCharge(0),
// TH1D with the observables for the event selection.
  fHistoCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoNumberOfTracksBeforeTrackSelection(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHNOTNumberOfBins(30000),
  fHNOTMax(30000.),
  fHistoVertexXBeforeSelection(NULL),
  fHistoVertexXAfterSelection(NULL),
  fHistoVertexYBeforeSelection(NULL),
  fHistoVertexYAfterSelection(NULL),
  fHistoVertexZBeforeSelection(NULL),
  fHistoVertexZAfterSelection(NULL),
// TH1D with the observables for the track selection.
  fHistoPtBeforeSelection(NULL),
  fHistoPtAfterSelection(NULL),
  fHistoEtaBeforeSelection(NULL),
  fHistoEtaAfterSelection(NULL),
  fHistoPhiBeforeSelection(NULL),
  fHistoPhiAfterSelection(NULL),
  fHistoTPCClustersBeforeSelection(NULL),
  fHistoTPCClustersAfterSelection(NULL),
  fHistoTPCChiSquareBeforeSelection(NULL),
  fHistoTPCChiSquareAfterSelection(NULL),
  fHistoDCAXYBeforeSelection(NULL),
  fHistoDCAXYAfterSelection(NULL),
  fHistoDCAZBeforeSelection(NULL),
  fHistoDCAZAfterSelection(NULL),
  fHistoChargeBeforeSelection(NULL),
  fHistoChargeAfterSelection(NULL),
// TProfiles with the final multiparticle correlations.
  fProfileTwoParticleCorrelations(NULL),
  fProfileFourParticleCorrelations(NULL),
  fProfileFourParticleCorrelationsCrossCheck(NULL),
  fProfileSixParticleCorrelations(NULL),
  fProfileTwoParticleCorrelationsNestedLoops(NULL),
  fProfileFourParticleCorrelationsNestedLoops(NULL)
{
/* Dummy constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Initialise 'fQvectors' to zero.
  InitialiseArraysOfQvectors();
}

//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
  AliAnalysisTaskSE(name),
// Structure of the output file.
  fMainList(NULL),
  fQAListBeforeSelection(NULL),
  fQAListAfterSelection(NULL),
  fListCorrelations(NULL),
// General parameters.
  fMaxNumberOfParticlesInCorrelations(8),
  fHighestFlowHarmonic(6),
  fUseParticleWeights(kFALSE),
  fCrossCheckFourParticleCorrelations(kFALSE),
  fCrossCheckWithNestedLoops(kFALSE),
// Type of files used in the analysis.
  fProcessOnlyAOD(kFALSE),
  fProcessOnlyMC(kFALSE),
  fProcessBothMCandAOD(kFALSE),
// Determination of the centrality.
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE),
  fCentralityMin(0.),
  fCentralityMax(100.),
// Event selection.
  fCutOnVertexX(kFALSE),
  fVertexMinX(-44.),
  fVertexMaxX(-44.),
  fCutOnVertexY(kFALSE),
  fVertexMinY(-44.),
  fVertexMaxY(-44.),
  fCutOnVertexZ(kFALSE),
  fVertexMinZ(-10.),
  fVertexMaxZ(10.),
  fNumberOfTracksMin(6),
// Track selection.
  fCutOnPt(kFALSE),
  fPtMin(0.2),
  fPtMax(5.),
  fCutOnEta(kFALSE),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fFilter(128),
  fCutOnNumberOfTPC(kFALSE),
  fNumberOfTPCMin(70),
  fCutOnChiSquarePInTPC(kFALSE),
  fChiSquarePInTPCMin(0.1),
  fChiSquarePInTPCMax(4.),
  fCutOnDCA(kFALSE),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
  fCutOnCharge(kFALSE),
  fCharge(0),
// TH1D with the observables for the event selection.
  fHistoCentrality(NULL),
  fHistoInitialNumberOfTracks(NULL),
  fHistoNumberOfTracksBeforeTrackSelection(NULL),
  fHistoFinalNumberOfTracks(NULL),
  fHNOTNumberOfBins(30000),
  fHNOTMax(30000.),
  fHistoVertexXBeforeSelection(NULL),
  fHistoVertexXAfterSelection(NULL),
  fHistoVertexYBeforeSelection(NULL),
  fHistoVertexYAfterSelection(NULL),
  fHistoVertexZBeforeSelection(NULL),
  fHistoVertexZAfterSelection(NULL),
// TH1D with the observables for the track selection.
  fHistoPtBeforeSelection(NULL),
  fHistoPtAfterSelection(NULL),
  fHistoEtaBeforeSelection(NULL),
  fHistoEtaAfterSelection(NULL),
  fHistoPhiBeforeSelection(NULL),
  fHistoPhiAfterSelection(NULL),
  fHistoTPCClustersBeforeSelection(NULL),
  fHistoTPCClustersAfterSelection(NULL),
  fHistoTPCChiSquareBeforeSelection(NULL),
  fHistoTPCChiSquareAfterSelection(NULL),
  fHistoDCAXYBeforeSelection(NULL),
  fHistoDCAXYAfterSelection(NULL),
  fHistoDCAZBeforeSelection(NULL),
  fHistoDCAZAfterSelection(NULL),
  fHistoChargeBeforeSelection(NULL),
  fHistoChargeAfterSelection(NULL),
// TProfiles with the final multiparticle correlations.
  fProfileTwoParticleCorrelations(NULL),
  fProfileFourParticleCorrelations(NULL),
  fProfileFourParticleCorrelationsCrossCheck(NULL),
  fProfileSixParticleCorrelations(NULL),
  fProfileTwoParticleCorrelationsNestedLoops(NULL),
  fProfileFourParticleCorrelationsNestedLoops(NULL)
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

// Define the procedure to follow if non-unit particle weights are used.
  /*if (fUseParticleWeights)
  {
    // TBA, needed for data periods after 2010.
  }*/
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

// Book the histograms in all the daughter lists.
  this->BookQAListBeforeSelection();
  this->BookQAListAfterSelection();
  this->BookListCorrelations();

// Continue to avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1, fMainList);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the analysis for each event. */
/// TBA: non-unit particle weight for non-uniform acceptance.
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";

// Select the type of file for the analysis (MC/AOD) from TaskSE.
  AliAODEvent *currentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent()); // Pointer to an AOD event.
  AliMCEvent *currentMCEvent = MCEvent(); // Pointer to a Monte Carlo event.

  if ((Int_t)fProcessOnlyAOD + (Int_t)fProcessOnlyMC + (Int_t)fProcessBothMCandAOD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one fProcess must be kTRUE in 'SetAnalysisType'.");
  }
  else if (fProcessOnlyAOD) {AnalyseAODevent(currentAODEvent);}
  else if (fProcessOnlyMC) {AnalyseMCevent(currentMCEvent);}
  else if (fProcessBothMCandAOD) {Fatal(sMethodName.Data(),"ERROR: TBA.");}

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
/* Initialise all the elements in 'fQvectors' to zero. */
  for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    }
  }
}

//######################################################################################//
// Methods called in 'UserCreateOutputObjects'.
//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists in the output file */
// Check if the mother list exists.
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()";
  if (!fMainList) {Fatal(sMethodName.Data(), "Error: 'fMainList' is NULL.");}

// Daughter list with the observables involved in the event selection.
  fQAListBeforeSelection = new TList();
  fQAListBeforeSelection->SetName("fQAListBeforeSelection");
  fQAListBeforeSelection->SetOwner(kTRUE);
  fMainList->Add(fQAListBeforeSelection);

// Daughter list with the observables involved in the track selection.
  fQAListAfterSelection = new TList();
  fQAListAfterSelection->SetName("fQAListAfterSelection");
  fQAListAfterSelection->SetOwner(kTRUE);
  fMainList->Add(fQAListAfterSelection);

// Daughter list with the multiparticle correlations.
  fListCorrelations = new TList();
  fListCorrelations->SetName("fListCorrelations");
  fListCorrelations->SetOwner(kTRUE);
  fMainList->Add(fListCorrelations);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookQAListBeforeSelection()
{
/* Book the TH1D with the observables involved in the event selection. */
// Distribution of the centrality.
  fHistoCentrality = new TH1D("fHistoCentrality", "Distribution of the centrality before the event selection", 100, 0., 100.);
  fHistoCentrality->SetStats(kTRUE);
  fHistoCentrality->GetXaxis()->SetTitle("Centrality percentile");
  fQAListBeforeSelection->Add(fHistoCentrality);

// Distribution of the initial number of tracks.
  fHistoInitialNumberOfTracks = new TH1I("fHistoInitialNumberOfTracks", "Initial distribution of the number of tracks", fHNOTNumberOfBins, 0., fHNOTMax);
  fHistoInitialNumberOfTracks->SetStats(kTRUE);
  fHistoInitialNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fQAListBeforeSelection->Add(fHistoInitialNumberOfTracks);

// Distribution of the number of tracks before the track selection.
  fHistoNumberOfTracksBeforeTrackSelection = new TH1I("fHistoNumberOfTracksBeforeTrackSelection", "Distribution of the number of tracks before the track selection", fHNOTNumberOfBins, 0., fHNOTMax);
  fHistoNumberOfTracksBeforeTrackSelection->SetStats(kTRUE);
  fHistoNumberOfTracksBeforeTrackSelection->GetXaxis()->SetTitle("Number of tracks");
  fQAListBeforeSelection->Add(fHistoNumberOfTracksBeforeTrackSelection);

// Distribution of the x-position of the PV.
  fHistoVertexXBeforeSelection = new TH1D("fHistoVertexXBeforeSelection", "Distribution of PV_{x} before the selection", 1000, -20., 20.);
  fHistoVertexXBeforeSelection->SetStats(kTRUE);
  fHistoVertexXBeforeSelection->GetXaxis()->SetTitle("PV_{x}");
  fQAListBeforeSelection->Add(fHistoVertexXBeforeSelection);

// Distribution of the y-position of the PV.
  fHistoVertexYBeforeSelection = new TH1D("fHistoVertexYBeforeSelection", "Distribution of PV_{y} before the selection", 1000, -20., 20.);
  fHistoVertexYBeforeSelection->SetStats(kTRUE);
  fHistoVertexYBeforeSelection->GetXaxis()->SetTitle("PV_{y}");
  fQAListBeforeSelection->Add(fHistoVertexYBeforeSelection);

// Distribution of the z-position of the PV.
  fHistoVertexZBeforeSelection = new TH1D("fHistoVertexZBeforeSelection", "Distribution of PV_{z} before the selection", 1000, -20., 20.);
  fHistoVertexZBeforeSelection->SetStats(kTRUE);
  fHistoVertexZBeforeSelection->GetXaxis()->SetTitle("PV_{z}");
  fQAListBeforeSelection->Add(fHistoVertexZBeforeSelection);

// Distribution of the transverse momentum.
  fHistoPtBeforeSelection = new TH1D("fHistoPtBeforeSelection", "Distribution of p_{T} before the track selection", 1000, 0., 20.);
  fHistoPtBeforeSelection->SetStats(kTRUE);
  fHistoPtBeforeSelection->GetXaxis()->SetTitle("p_{T}");
  fQAListBeforeSelection->Add(fHistoPtBeforeSelection);

// Distribution of the pseudorapidity.
  fHistoEtaBeforeSelection = new TH1D("fHistoEtaBeforeSelection", "Distribution of #eta before the track selection", 1000, -1., 1.);
  fHistoEtaBeforeSelection->SetStats(kTRUE);
  fHistoEtaBeforeSelection->GetXaxis()->SetTitle("#eta");
  fQAListBeforeSelection->Add(fHistoEtaBeforeSelection);

// Distribution of the azimuthal angles.
  fHistoPhiBeforeSelection = new TH1D("fHistoPhiBeforeSelection", "Distributiion of #phi before the track selection", 1000, 0., 6.3);
  fHistoPhiBeforeSelection->SetStats(kTRUE);
  fHistoPhiBeforeSelection->GetXaxis()->SetTitle("#phi");
  fQAListBeforeSelection->Add(fHistoPhiBeforeSelection);

// Distribution of the number of TPC clusters.
  fHistoTPCClustersBeforeSelection = new TH1I("fHistoTPCClustersBeforeSelection", "Distribution of the number of TPC clusters before the track selection", 1000, 0., 170.);
  fHistoTPCClustersBeforeSelection->SetStats(kTRUE);
  fHistoTPCClustersBeforeSelection->GetXaxis()->SetTitle("Number of TPC clusters");
  fQAListBeforeSelection->Add(fHistoTPCClustersBeforeSelection);

// Distribution of the chi square of the track momentum in the TPC.
  fHistoTPCChiSquareBeforeSelection = new TH1D("fHistoTPCChiSquareBeforeSelection", "Distribution of the #chi^{2} of the track momentum in the TPC before the track selection", 1000, 0., 20.);
  fHistoTPCChiSquareBeforeSelection->SetStats(kTRUE);
  fHistoTPCChiSquareBeforeSelection->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fQAListBeforeSelection->Add(fHistoTPCChiSquareBeforeSelection);

// Distribution of the xy-plane of the DCA.
  fHistoDCAXYBeforeSelection = new TH1D("fHistoDCAXYBeforeSelection", "Distribution of DCA_{xy} before the track selection", 1000, 0., 10.);
  fHistoDCAXYBeforeSelection->SetStats(kTRUE);
  fHistoDCAXYBeforeSelection->GetXaxis()->SetTitle("DCA_{xy}");
  fQAListBeforeSelection->Add(fHistoDCAXYBeforeSelection);

// Distribution of the z-coordinate of the DCA.
  fHistoDCAZBeforeSelection = new TH1D("fHistoDCAZBeforeSelection", "Distribution of DCA_{z} before the track selection", 1000, 0., 10.);
  fHistoDCAZBeforeSelection->SetStats(kTRUE);
  fHistoDCAZBeforeSelection->GetXaxis()->SetTitle("DCA_{z}");
  fQAListBeforeSelection->Add(fHistoDCAZBeforeSelection);

// Distribution of the electric charge.
  fHistoChargeBeforeSelection = new TH1I("fHistoChargeBeforeSelection", "Distribution of the electric charge before the track selection", 2, -2, 2);
  fHistoChargeBeforeSelection->SetStats(kTRUE);
  fHistoChargeBeforeSelection->GetXaxis()->SetTitle("Charge");
  fQAListBeforeSelection->Add(fHistoChargeBeforeSelection);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookQAListAfterSelection()
{
/* Book the TH1D with the observables involved in the track selection. */
// Final number of tracks after both the event and track selection.
  fHistoFinalNumberOfTracks = new TH1I("fHistoFinalNumberOfTracks", "Final distribution of the number of tracks", fHNOTNumberOfBins, 0., fHNOTMax);
  fHistoFinalNumberOfTracks->SetStats(kTRUE);
  fHistoFinalNumberOfTracks->GetXaxis()->SetTitle("Number of tracks");
  fQAListAfterSelection->Add(fHistoFinalNumberOfTracks);

// Distribution of the x-position of the PV.
  fHistoVertexXAfterSelection = new TH1D("fHistoVertexXAfterSelection", "Distribution of PV_{x} after the full selection", 1000, -20., 20.);
  fHistoVertexXAfterSelection->SetStats(kTRUE);
  fHistoVertexXAfterSelection->GetXaxis()->SetTitle("PV_{x}");
  fQAListAfterSelection->Add(fHistoVertexXAfterSelection);

// Distribution of the y-position of the PV.
  fHistoVertexYAfterSelection = new TH1D("fHistoVertexYAfterSelection", "Distribution of PV_{y} after the full selection", 1000, -20., 20.);
  fHistoVertexYAfterSelection->SetStats(kTRUE);
  fHistoVertexYAfterSelection->GetXaxis()->SetTitle("PV_{y}");
  fQAListAfterSelection->Add(fHistoVertexYAfterSelection);

// Distribution of the z-position of the PV.
  fHistoVertexZAfterSelection = new TH1D("fHistoVertexZAfterSelection", "Distribution of PV_{z} after the full selection", 1000, -20., 20.);
  fHistoVertexZAfterSelection->SetStats(kTRUE);
  fHistoVertexZAfterSelection->GetXaxis()->SetTitle("PV_{z}");
  fQAListAfterSelection->Add(fHistoVertexZAfterSelection);

// Distribution of the transverse momentum.
  fHistoPtAfterSelection = new TH1D("fHistoPtAfterSelection", "Distribution of p_{T} after the full selection", 1000, 0., 20.);
  fHistoPtAfterSelection->SetStats(kTRUE);
  fHistoPtAfterSelection->GetXaxis()->SetTitle("p_{T}");
  fQAListAfterSelection->Add(fHistoPtAfterSelection);

// Distribution of the pseudorapidity.
  fHistoEtaAfterSelection = new TH1D("fHistoEtaAfterSelection", "Distribution of #eta after the full selection", 1000, -1., 1.);
  fHistoEtaAfterSelection->SetStats(kTRUE);
  fHistoEtaAfterSelection->GetXaxis()->SetTitle("#eta");
  fQAListAfterSelection->Add(fHistoEtaAfterSelection);

// Distribution of the azimuthal angles.
  fHistoPhiAfterSelection = new TH1D("fHistoPhiAfterSelection", "Distribution of #phi after the full selection", 1000, 0., 6.3);
  fHistoPhiAfterSelection->SetStats(kTRUE);
  fHistoPhiAfterSelection->GetXaxis()->SetTitle("#phi");
  fQAListAfterSelection->Add(fHistoPhiAfterSelection);

// Distribution of the number of TPC clusters.
  fHistoTPCClustersAfterSelection = new TH1I("fHistoTPCClustersAfterSelection", "Distribution of the number of TPC clusters after the full selection", 1000, 0., 170.);
  fHistoTPCClustersAfterSelection->SetStats(kTRUE);
  fHistoTPCClustersAfterSelection->GetXaxis()->SetTitle("Number of TPC clusters");
  fQAListAfterSelection->Add(fHistoTPCClustersAfterSelection);

// Distribution of the chi square of the track momentum in the TPC.
  fHistoTPCChiSquareAfterSelection = new TH1D("fHistoTPCChiSquareAfterSelection", "Distribution of the #chi^{2} of the track momentum in the TPC after the full selection", 1000, 0., 20.);
  fHistoTPCChiSquareAfterSelection->SetStats(kTRUE);
  fHistoTPCChiSquareAfterSelection->GetXaxis()->SetTitle("#chi^{2}/NDF in TPC");
  fQAListAfterSelection->Add(fHistoTPCChiSquareAfterSelection);

// Distribution of the xy-plane of the DCA.
  fHistoDCAXYAfterSelection = new TH1D("fHistoDCAXYAfterSelection", "Distribution of DCA_{xy} after the full selection", 1000, 0., 10.);
  fHistoDCAXYAfterSelection->SetStats(kTRUE);
  fHistoDCAXYAfterSelection->GetXaxis()->SetTitle("DCA_{xy}");
  fQAListAfterSelection->Add(fHistoDCAXYAfterSelection);

// Distribution of the z-coordinate of the DCA.
  fHistoDCAZAfterSelection = new TH1D("fHistoDCAZAfterSelection", "Distribution of DCA_{z} after the full selection", 1000, 0., 10.);
  fHistoDCAZAfterSelection->SetStats(kTRUE);
  fHistoDCAZAfterSelection->GetXaxis()->SetTitle("DCA_{z}");
  fQAListAfterSelection->Add(fHistoDCAZAfterSelection);

// Distribution of the electric charge.
  fHistoChargeAfterSelection = new TH1I("fHistoChargeAfterSelection", "Distribution of the electric charge after the full selection", 2, -2, 2);
  fHistoChargeAfterSelection->SetStats(kTRUE);
  fHistoChargeAfterSelection->GetXaxis()->SetTitle("Charge");
  fQAListAfterSelection->Add(fHistoChargeAfterSelection);
}

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookListCorrelations()
{
/* Book the TProfiles with the multiparticle correlations. */
// 2-particle correlations.
  fProfileTwoParticleCorrelations = new TProfile("fProfileTwoParticleCorrelations", "2-particle correlations", 6, 0., 6.);
  fProfileTwoParticleCorrelations->SetStats(kTRUE);
  fProfileTwoParticleCorrelations->Sumw2();
  fProfileTwoParticleCorrelations->GetXaxis()->SetTitle("n");
  fProfileTwoParticleCorrelations->GetYaxis()->SetTitle("#LT#LT2#GT#GT_{n,-n}");
  fListCorrelations->Add(fProfileTwoParticleCorrelations);

  if (fCrossCheckWithNestedLoops)
  {
    fProfileTwoParticleCorrelationsNestedLoops = new TProfile("fProfileTwoParticleCorrelationsNestedLoops", "2-particle correlations with nested loops", 6, 0., 6.);
    fProfileTwoParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileTwoParticleCorrelationsNestedLoops->Sumw2();
    fProfileTwoParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("n");
    fProfileTwoParticleCorrelationsNestedLoops->GetYaxis()->SetTitle("#LT#LT2#GT#GT_{n,-n}");
    fListCorrelations->Add(fProfileTwoParticleCorrelationsNestedLoops);
  }

// 4-particle correlations.
  fProfileFourParticleCorrelations = new TProfile("fProfileFourParticleCorrelations", "4-particle correlations", 21, 0., 21.);
  fProfileFourParticleCorrelations->SetStats(kTRUE);
  fProfileFourParticleCorrelations->Sumw2();
  fProfileFourParticleCorrelations->GetXaxis()->SetTitle("(m,n)");
  fProfileFourParticleCorrelations->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
  fListCorrelations->Add(fProfileFourParticleCorrelations);

  if (fCrossCheckFourParticleCorrelations)
  {
    fProfileFourParticleCorrelationsCrossCheck = new TProfile("fProfileFourParticleCorrelationsCrossCheck", "4-particle correlations for cross-check", 15, 0., 15.);
    fProfileFourParticleCorrelationsCrossCheck->SetStats(kTRUE);
    fProfileFourParticleCorrelationsCrossCheck->Sumw2();
    fProfileFourParticleCorrelationsCrossCheck->GetXaxis()->SetTitle("(m,n)");
    fProfileFourParticleCorrelationsCrossCheck->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
    fListCorrelations->Add(fProfileFourParticleCorrelationsCrossCheck);
  }

  if (fCrossCheckWithNestedLoops)
  {
    fProfileFourParticleCorrelationsNestedLoops = new TProfile("fProfileFourParticleCorrelationsNestedLoops", "4-particle correlations with nested loops", 21, 0., 21.);
    fProfileFourParticleCorrelationsNestedLoops->SetStats(kTRUE);
    fProfileFourParticleCorrelationsNestedLoops->Sumw2();
    fProfileFourParticleCorrelationsNestedLoops->GetXaxis()->SetTitle("(m,n)");
    fProfileFourParticleCorrelationsNestedLoops->GetYaxis()->SetTitle("#LT#LT4#GT#GT_{m,n,-m,-n}");
    fListCorrelations->Add(fProfileFourParticleCorrelationsNestedLoops);
  }

// 6-particle correlations.
  fProfileSixParticleCorrelations = new TProfile("fProfileSixParticleCorrelations", "6-particle correlations", 20, 0., 20.);
  fProfileSixParticleCorrelations->SetStats(kTRUE);
  fProfileSixParticleCorrelations->Sumw2();
  fProfileSixParticleCorrelations->GetXaxis()->SetTitle("(l,m,n)");
  fProfileSixParticleCorrelations->GetYaxis()->SetTitle("#LT#LT6#GT#GT_{l,m,n,-l,-m,-n}");
  fListCorrelations->Add(fProfileSixParticleCorrelations);
}

//######################################################################################//
// Methods called in 'UserExec'.
//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *aAODevent)
{
/* Execute the analysis for the provided AOD event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent(AliAODEvent *aAODevent)";

// Check if there is an event or not.
  if (!aAODevent) {Fatal(sMethodName.Data(), "ERROR: no AOD event found.");}

// Select the detector to use for the estimation of the centrality.
  TString centralityEstimator = "centralityEstimator";  // Name of the detector used for the centrality estimation.
  if ((Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1)
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be selected in 'SetCentralityEstimation'.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Determine if the event belongs to this centrality range.
  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if (!ams) {return;} // Protection against NULL pointer.
  Double_t aCentrality = ams->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));  // Centrality of the given event.
  if ((aCentrality >= fCentralityMin) && (aCentrality < fCentralityMax))
  {
    fHistoCentrality->Fill(aCentrality);
  }
  else {return;}  // This event does not belong to this centrality range.

// Get the number of tracks before the event selection.
  long long initialNumberOfTracks = aAODevent->GetNumberOfTracks();
  fHistoInitialNumberOfTracks->Fill(initialNumberOfTracks);

// Cuts on the position of the Primary Vertex.
  AliAODVertex *avtx = (AliAODVertex*)aAODevent->GetPrimaryVertex();  // 3d position of the PV.
  fHistoVertexXBeforeSelection->Fill(avtx->GetX());
  fHistoVertexYBeforeSelection->Fill(avtx->GetY());
  fHistoVertexZBeforeSelection->Fill(avtx->GetZ());

  if (fCutOnVertexX)
  {
    if ((avtx->GetX() < fVertexMinX) || (avtx->GetX() > fVertexMaxX)) {return;}
  }
  if (fCutOnVertexY)
  {
    if ((avtx->GetY() < fVertexMinY) || (avtx->GetY() > fVertexMaxY)) {return;}
  }
  if (fCutOnVertexZ)
  {
    if ((avtx->GetZ() < fVertexMinZ) || (avtx->GetZ() > fVertexMaxZ)) {return;}
  }

/// TBA: more event cuts?

// Preparations for the track selection.
  long long numberOfTracksBeforeTrackSelection = aAODevent->GetNumberOfTracks();  // Number of tracks before the track selection.
  fHistoNumberOfTracksBeforeTrackSelection->Fill(numberOfTracksBeforeTrackSelection);
  long long finalNumberOfTracks = 0;  // Number of tracks after the full selection.
  Int_t *IsTrackSelected = new Int_t[numberOfTracksBeforeTrackSelection](); // Flag to indicate a track passed the track selection (1) or not (0).

  Double_t pT = 0.; // Transverse momentum.
  Double_t eta = 0.;  // Pseudorapidity.
  Double_t phi = 0.;  // Azimuthal angle.
  Int_t numberOfTPCClusters = 0;  // Number of TPC clusters.
  Double_t chiSquareInTPC = 0.; // Chi square of the track momentum in the TPC.
  Double_t DCAx = 0.; // x-value of the DCA.
  Double_t DCAy = 0.; // y-value of the DCA.
  Double_t DCAz = 0.; // z-value of the DCA.
  Double_t DCAxy = 0.;  // xy-value of the DCA.
  Int_t charge = 0; // Electric charge.

// Look at each track in the event to mark them as selected or not.
  for (long long iTrack = 0; iTrack < numberOfTracksBeforeTrackSelection; iTrack++)
  {
    AliAODTrack *currentTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack));  // Pointer to the AOD track.
    if (!currentTrack) {continue;}  // Protection against NULL pointer.
    if (!currentTrack->TestFilterBit(fFilter)) {continue;}  // Filter bit 128 denotes TPC-only tracks.

  // Get all the observables for the track selection.
    pT = currentTrack->Pt();
    eta = currentTrack->Eta();
    phi = currentTrack->Phi();
    numberOfTPCClusters = currentTrack->GetTPCNcls();
    chiSquareInTPC = currentTrack->Chi2perNDF();
    DCAx = currentTrack->XAtDCA();
    DCAy = currentTrack->YAtDCA();
    DCAz = currentTrack->ZAtDCA();
    charge = currentTrack->Charge();

    DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));

  // Fill the histograms before the track selection.
    fHistoPtBeforeSelection->Fill(pT);
    fHistoEtaBeforeSelection->Fill(eta);
    fHistoPhiBeforeSelection->Fill(phi);
    fHistoTPCClustersBeforeSelection->Fill(numberOfTPCClusters);
    fHistoTPCChiSquareBeforeSelection->Fill(chiSquareInTPC);
    fHistoDCAXYBeforeSelection->Fill(DCAxy);
    fHistoDCAZBeforeSelection->Fill(DCAz);
    fHistoChargeBeforeSelection->Fill(charge);

  // Apply the track selection to the provided track.
    if (ApplyTrackSelection(pT, eta, numberOfTPCClusters, chiSquareInTPC, DCAxy, DCAz, charge))  // The track passed the selection.
    {
      IsTrackSelected[iTrack] = 1;
      finalNumberOfTracks++;
    }
    else {IsTrackSelected[iTrack] = 0;}  // The track failed the selection.
  }

// Remove the events with not enough tracks for the event weight.
  if (finalNumberOfTracks <= fNumberOfTracksMin) {return;}

// Fill all the event histograms after the full selection.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);
  fHistoVertexXAfterSelection->Fill(avtx->GetX());
  fHistoVertexYAfterSelection->Fill(avtx->GetY());
  fHistoVertexZAfterSelection->Fill(avtx->GetZ());

// Define the arrays for the azimuthal angles and particle weights to use in the analysis.
  Double_t *phiArray = new Double_t[finalNumberOfTracks](); // Azimuthal angles.
  Double_t *particleWeightArray = new Double_t[finalNumberOfTracks](); // Particle weights.
  Int_t indexInNewArrays = 0; // New index of the track if it passed the selection.

// Loop over the tracks to keep only the selected ones.
  for (long long iTrack = 0; iTrack < numberOfTracksBeforeTrackSelection; iTrack++)
  {
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack));  // Pointer to the AOD track.
    if (!aTrack) {continue;}  // Protection against NULL pointer.
    if (!aTrack->TestFilterBit(fFilter)) {continue;}  // Filter bit 128 denotes TPC-only tracks.

    if (IsTrackSelected[iTrack] == 1) // The particle passed the selection.
    {
    // Get all the observables used in the track selection.
      pT = aTrack->Pt();
      eta = aTrack->Eta();
      phiArray[indexInNewArrays] = aTrack->Phi();
      numberOfTPCClusters = aTrack->GetTPCNcls();
      chiSquareInTPC = aTrack->Chi2perNDF();
      DCAx = aTrack->XAtDCA();
      DCAy = aTrack->YAtDCA();
      DCAz = aTrack->ZAtDCA();
      DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));
      charge = aTrack->Charge();

      if (fUseParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA.");}
      else {particleWeightArray[indexInNewArrays] = 1.;}

    // Fill all the track histograms after the full selection.
      fHistoPtAfterSelection->Fill(pT);
      fHistoEtaAfterSelection->Fill(eta);
      fHistoPhiAfterSelection->Fill(phiArray[indexInNewArrays]);
      fHistoTPCClustersAfterSelection->Fill(numberOfTPCClusters);
      fHistoTPCChiSquareAfterSelection->Fill(chiSquareInTPC);
      fHistoDCAXYAfterSelection->Fill(DCAxy);
      fHistoDCAZAfterSelection->Fill(DCAz);
      fHistoChargeAfterSelection->Fill(charge);

    // Increase the value of 'indexInNewArrays' by one.
      indexInNewArrays++;
    }
    else {continue;}
  }

// Calculate the Q-vectors for the current event.
  CalculateQvectors(finalNumberOfTracks, phiArray, particleWeightArray);

// Compute all the multiparticle correlations for the current event.
  ComputeMultiparticleCorrelations(finalNumberOfTracks, phiArray, particleWeightArray);

// Reset everything to zero for the next event.
  numberOfTracksBeforeTrackSelection = 0;
  finalNumberOfTracks = 0;
  delete [] IsTrackSelected;
  pT = 0.;
  eta = 0.;
  phi = 0.;
  numberOfTPCClusters = 0;
  chiSquareInTPC = 0.;
  DCAx = 0.;
  DCAy = 0.;
  DCAz = 0.;
  charge = 0;
  delete [] phiArray;
  delete [] particleWeightArray;
  indexInNewArrays = 0;
}

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
  if (!fProcessOnlyMC)
  {
    AliMultSelection *ams = (AliMultSelection*)aMCevent->FindListObject("MultSelection");
    if (!ams) {return;} // Protection against NULL pointer.
    Double_t aCentrality = ams->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));  // Centrality of the given event.
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
  fHistoVertexXBeforeSelection->Fill(avtx->GetX());
  fHistoVertexYBeforeSelection->Fill(avtx->GetY());
  fHistoVertexZBeforeSelection->Fill(avtx->GetZ());

  if (fCutOnVertexX)
  {
    if ((avtx->GetX() < fVertexMinX) || (avtx->GetX() > fVertexMaxX)) {return;}
  }
  if (fCutOnVertexY)
  {
    if ((avtx->GetY() < fVertexMinY) || (avtx->GetY() > fVertexMaxY)) {return;}
  }
  if (fCutOnVertexZ)
  {
    if ((avtx->GetZ() < fVertexMinZ) || (avtx->GetZ() > fVertexMaxZ)) {return;}
  }

/// TBA: more event cuts?

// Preparations for the track selection.
  long long numberOfTracksBeforeTrackSelection = aMCevent->GetNumberOfTracks();  // Number of tracks before the track selection.
  fHistoNumberOfTracksBeforeTrackSelection->Fill(numberOfTracksBeforeTrackSelection);
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
    fHistoPtBeforeSelection->Fill(pT);
    fHistoEtaBeforeSelection->Fill(eta);
    fHistoPhiBeforeSelection->Fill(phi);
    fHistoChargeBeforeSelection->Fill(charge);

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

// Remove the events with not enough tracks for the event weight.
  if (finalNumberOfTracks <= fNumberOfTracksMin) {return;}

// Fill all the event histograms after the full selection.
  fHistoFinalNumberOfTracks->Fill(finalNumberOfTracks);
  fHistoVertexXAfterSelection->Fill(avtx->GetX());
  fHistoVertexYAfterSelection->Fill(avtx->GetY());
  fHistoVertexZAfterSelection->Fill(avtx->GetZ());

// Define the arrays for the azimuthal angles and particle weights to use in the analysis.
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
      eta = aTrack->Eta();
      phiArray[indexInNewArrays] = aTrack->Phi();
      charge = aTrack->Charge();

      if (fUseParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA.");}
      else {particleWeightArray[indexInNewArrays] = 1.;}

    // Fill all the track histograms after the full selection.
      fHistoPtAfterSelection->Fill(pT);
      fHistoEtaAfterSelection->Fill(eta);
      fHistoPhiAfterSelection->Fill(phiArray[indexInNewArrays]);
      fHistoChargeAfterSelection->Fill(charge);

    // Increase the value of 'indexInNewArrays' by one.
      indexInNewArrays++;
    }
    else {continue;}
  }

// Calculate the Q-vectors for the current event.
  CalculateQvectors(finalNumberOfTracks, phiArray, particleWeightArray);

// Compute all the multiparticle correlations for the current event.
  ComputeMultiparticleCorrelations(finalNumberOfTracks, phiArray, particleWeightArray);

// Reset everything to zero for the next event.
  numberOfTracksBeforeTrackSelection = 0;
  finalNumberOfTracks = 0;
  delete [] IsTrackSelected;
  pT = 0.;
  eta = 0.;
  phi = 0.;
  charge = 0;
  delete [] phiArray;
  delete [] particleWeightArray;
  indexInNewArrays = 0;
}

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(Double_t momentum, Double_t pseudorapidity, Int_t NclustersInTPC, Double_t TPCchiSquare, Double_t xyDCA, Double_t zDCA, Int_t eCharge)
{
/* Apply the track selection to the arguments and return if it is passed or not. */
  Bool_t testOfPt = kTRUE;  // Cut on the transverse momentum.
  Bool_t testOfEta = kTRUE; // Cut on the pseudorapidity.
  Bool_t testOfNumberOfTPC = kTRUE; // Cut on the number of TPC clusters.
  Bool_t testOfChiSquareTPC = kTRUE;  // Cut on the chi^2 of the momentum in the TPC.
  Bool_t testOfDCAxy = kTRUE; // Cut on the DCA of the track in the xy-plane.
  Bool_t testOfDCAz = kTRUE;  // Cut on the DCA of the track along the z-direction.
  Bool_t testOfCharge = kTRUE;  // Cut on the electric charge.

  if (fCutOnPt) {testOfPt = (fPtMin <= momentum) && (momentum <= fPtMax);}
  if (fCutOnEta) {testOfEta = (fEtaMin <= pseudorapidity) && (pseudorapidity <= fEtaMax);}
  if (fCutOnNumberOfTPC) {testOfNumberOfTPC = (fNumberOfTPCMin < NclustersInTPC);}
  if (fCutOnChiSquarePInTPC) {testOfChiSquareTPC = (fChiSquarePInTPCMin <= TPCchiSquare) && (TPCchiSquare <= fChiSquarePInTPCMax);}
  if (fCutOnDCA) {testOfDCAxy = (xyDCA < fDCAxyMax); testOfDCAz = (zDCA < fDCAzMax);}
  if (fCutOnCharge) {testOfCharge = (eCharge == fCharge);}

  return testOfPt && testOfEta && testOfNumberOfTPC && testOfChiSquareTPC && testOfDCAxy && testOfDCAz && testOfCharge;
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
    if (fCrossCheckWithNestedLoops)
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
      if (fCrossCheckWithNestedLoops)
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
  if (fCrossCheckFourParticleCorrelations)
  {
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
  // Reset the variables for the next event.
    iBinFourCrossCheck =1;
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

//######################################################################################//
// Methods called in 'Terminate'.
//======================================================================================//


