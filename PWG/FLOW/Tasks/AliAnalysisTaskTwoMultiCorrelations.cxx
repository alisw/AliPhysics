
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
// Analysis task for the computation of the multiTrack correlations with different   //
// flow harmonics v_n. This script can takes the Monte Carlo simulations data (e.g.     //
// HIJING), as well as the experimental Pb-Pb data taken by the ALICE experiment.       //
// The current script computes the multiTrack correlators using the method of the    //
// Q-vectors for a maximum of 6 different harmonics and 8 particles).                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
// Version: 12.02.2019                                                                  //
//--------------------------------------------------------------------------------------//

#include "Riostream.h"
#include <vector>
#include "AliAnalysisTaskTwoMultiCorrelations.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODVertex.h"
#include "AliMCVertex.h"
#include "TList.h"
#include "TH1D.h"
#include "TFile.h"
#include "TComplex.h"
#include "TMath.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskTwoMultiCorrelations)

//======================================================================================//
// Main methods for AliAnalysisTaskSE.                                                  //
//--------------------------------------------------------------------------------------//
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights):
  AliAnalysisTaskSE(name),
// General parameters.
  fMaxNumberCorrelations(8),
  fMaxFlowHarmonic(6),
  fNumberHarmonicsInSC(3),
  fUseParticleWeights(kFALSE),
  fComputeNestedLoops(kFALSE),
// Structure of the output file.
  fMainList(NULL),
  fControlListEventCuts(NULL),
  fControlListTrackCuts(NULL),
  fListCorrelations(NULL),
  fListTwoParticles(NULL),
  fListThreeParticles(NULL),
  fListFourParticles(NULL),
  fListSixParticles(NULL),
// Control histograms for the distribution of observables for the event selection.
  fHistoCentrality(NULL),
  fHistoNumberOfTracksBeforeCuts(NULL),
  fHistoNumberOfTracksAfterEventCuts(NULL),
  fHistoVertexXBeforeCuts(NULL),
  fHistoVertexYBeforeCuts(NULL),
  fHistoVertexZBeforeCuts(NULL),
  fHistoVertexXAfterEventCuts(NULL),
  fHistoVertexYAfterEventCuts(NULL),
  fHistoVertexZAfterEventCuts(NULL),
// Control histrograms for the distribution of observables for the track selection.
  fHistoNumberOfTracksAfterAllCuts(NULL),
  fHistoPtBeforeCuts(NULL),
  fHistoPtAfterCuts(NULL),
  fHistoEtaBeforeCuts(NULL),
  fHistoEtaAfterCuts(NULL),
  fHistoPhiBeforeCuts(NULL),
  fHistoPhiAfterCuts(NULL),
  fHistoTPCClustersBeforeCuts(NULL),
  fHistoTPCClustersAfterCuts(NULL),
  fHistoDCAXYBeforeCuts(NULL),
  fHistoDCAZBeforeCuts(NULL),
  fHistoDCAXYAfterCuts(NULL),
  fHistoDCAZAfterCuts(NULL),
// TProfiles with the final multiTrack correlations.
  fProfileCosineEightParticles(NULL),
// Type of analysis.
  fAnalysisType(NULL),
  fProcessBothKineAndReco(kFALSE),
  fProcessOnlyKine(kFALSE),
  fProcessOnlyReco(kFALSE),
// Determination of the centrality.
  fCentralitySelection(NULL),
  fUseSPDForCentrality(kFALSE),
  fUseVZeroForCentrality(kFALSE),
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
// Track selection.
  fPtMin(0.2),
  fPtMax(5.),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fNumberOfTPCMin(70),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
// Harmonics.
  fHarmonicOne(2),
  fHarmonicTwo(-2),
  fHarmonicThree(3),
  fHarmonicFour(-3),
  fHarmonicFive(4),
  fHarmonicSix(-4),
  fHarmonicSeven(0),
  fHarmonicEight(0)
{
/* Constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Create the main list.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE);

// Define the input and output slots.
  DefineOutput(1, TList::Class());

// Initialise the fQvectors to zero.
  InitialiseArraysOfQvectors();

// Initialise the pointers of the TProfiles to NULL.
  InitialiseArraysOfTProfiles();

// Use the particle weights?
  if(useParticleWeights)
  {
    /* TBA, but not needed for LHC10h */
  }
} // End: AliAnalysisTaskTwoMultiCorrelations().

//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations():
  AliAnalysisTaskSE(),
// General parameters.
  fMaxNumberCorrelations(8),
  fMaxFlowHarmonic(6),
  fNumberHarmonicsInSC(3),
  fUseParticleWeights(kFALSE),
  fComputeNestedLoops(kFALSE),
// Structure of the output file.
  fMainList(NULL),
  fControlListEventCuts(NULL),
  fControlListTrackCuts(NULL),
  fListCorrelations(NULL),
  fListTwoParticles(NULL),
  fListThreeParticles(NULL),
  fListFourParticles(NULL),
  fListSixParticles(NULL),
// Control histograms for the distribution of observables for the event selection.
  fHistoCentrality(NULL),
  fHistoNumberOfTracksBeforeCuts(NULL),
  fHistoNumberOfTracksAfterEventCuts(NULL),
  fHistoVertexXBeforeCuts(NULL),
  fHistoVertexYBeforeCuts(NULL),
  fHistoVertexZBeforeCuts(NULL),
  fHistoVertexXAfterEventCuts(NULL),
  fHistoVertexYAfterEventCuts(NULL),
  fHistoVertexZAfterEventCuts(NULL),
// Control histrograms for the distribution of observables for the track selection.
  fHistoNumberOfTracksAfterAllCuts(NULL),
  fHistoPtBeforeCuts(NULL),
  fHistoPtAfterCuts(NULL),
  fHistoEtaBeforeCuts(NULL),
  fHistoEtaAfterCuts(NULL),
  fHistoPhiBeforeCuts(NULL),
  fHistoPhiAfterCuts(NULL),
  fHistoTPCClustersBeforeCuts(NULL),
  fHistoTPCClustersAfterCuts(NULL),
  fHistoDCAXYBeforeCuts(NULL),
  fHistoDCAZBeforeCuts(NULL),
  fHistoDCAXYAfterCuts(NULL),
  fHistoDCAZAfterCuts(NULL),
// TProfiles with the final multiTrack correlations.
  fProfileCosineEightParticles(NULL),
// Type of analysis.
  fAnalysisType(NULL),
  fProcessBothKineAndReco(kFALSE),
  fProcessOnlyKine(kFALSE),
  fProcessOnlyReco(kFALSE),
// Determination of the centrality.
  fCentralitySelection(NULL),
  fUseSPDForCentrality(kFALSE),
  fUseVZeroForCentrality(kFALSE),
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
// Track selection.
  fPtMin(0.2),
  fPtMax(5.),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fNumberOfTPCMin(70),
  fDCAxyMax(3.2),
  fDCAzMax(2.4),
// Harmonics.
  fHarmonicOne(2),
  fHarmonicTwo(-2),
  fHarmonicThree(3),
  fHarmonicFour(-3),
  fHarmonicFive(4),
  fHarmonicSix(-4),
  fHarmonicSeven(0),
  fHarmonicEight(0)
{
/* Dummy constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Initialise the fQvectors to zero.
  InitialiseArraysOfQvectors();

// Initialise the pointers of the TProfiles to NULL.
  InitialiseArraysOfTProfiles();

} // End: AliAnalysisTaskTwoMultiCorrelations().

//======================================================================================//
AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
/* Destructor of the class. */
  if (fMainList) {delete fMainList;}

} // End: ~AliAnalysisTaskTwoMultiCorrelations().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Steps to do at the beginning of the analysis. */
// Avoid name clashes.
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(kFALSE);

// Book all the lists.
  this->BookAllLists();

// Book the graphs in all the secondary lists.
  this->BookControlListEventCuts();
  this->BookControlListTrackCuts();
  this->BookListCorrelations();

// Initialise the pointers for the TStrings.
  fAnalysisType = new TString();
  fCentralitySelection = new TString();

// Still avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1,fMainList);

} // End: UserCreateOutputObjects().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Deal with the main analysis for each event. */
// TBA: non-unit particle weights for non-uniform acceptance.
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";

// Select the type of analysis (MC/AOD) from TaskSE.
  AliMCEvent *currentMCEvent = MCEvent(); // Pointer to a Monte Carlo event.
  AliAODEvent *currentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());  // Pointer to an AOD event.
  if (1 != (Int_t)fProcessBothKineAndReco + (Int_t)fProcessOnlyKine + (Int_t)fProcessOnlyReco)
  {
    Fatal(sMethodName.Data(), "ERROR: only one fProcess must be kTRUE in 'SetAnalysisType'.");
  }
  else if (fProcessBothKineAndReco) {*fAnalysisType = "MC_AOD";}
  else if (fProcessOnlyKine) {*fAnalysisType = "MC";}
  else if (fProcessOnlyReco) {*fAnalysisType = "AOD";}

// Call the specific analysis.
  if (fAnalysisType->EqualTo("AOD")) {AODanalysis(currentAODEvent);}
  else if (fAnalysisType->EqualTo("MC")) {MCanalysis(currentMCEvent);}
  else if (fAnalysisType->EqualTo("MC_AOD")) {Fatal(sMethodName.Data(),"ERROR: TBA.");}
  else {Fatal(sMethodName.Data(), "ERROR: fAnalysisType not defined.");}

// PostData.
  PostData(1, fMainList);

} // End: UserExec().

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
/* Post-analysis steps done once the run over the events is done. */
// Access the main list.
  fMainList = (TList*)GetOutputData(1);
  if (!fMainList) {exit(1);}

// Create the output file and save the main list inside.
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fMainList->Write(fMainList->GetName(),TObject::kSingleKey);
  delete outputFile;

} // End: Terminate().

//======================================================================================//
// Methods called in the constructor.                                                   //
//--------------------------------------------------------------------------------------//
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfQvectors()
{
/* Initialisation of the Q-vectors to zero. */
  for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    }
  }

} // End: InitialiseArraysOfQvectors().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfTProfiles()
{
/* Initialisation of the pointers to the TProfiles. */
  for (Int_t i = 0; i < 6; i++)
  {
    fProfileCosineTwoParticles[i] = NULL;
    fProfileCosineTwoNestedLoops[i] = NULL;
    fProfileCosineFourParticles[i] = NULL;
    fProfileCosineFourNestedLoops[i] = NULL;
  }

  for (Int_t j = 0; j < 5; j++)
  {
    fProfileTwoCosine[j] = NULL;
    fProfileTwoCosineNestedLoops[j] = NULL;
  }

  for (Int_t k = 0; k < 4; k++)
  {
    fProfileCosineThreeParticles[k] = NULL;
    fProfileCosineThreeNestedLoops[k] = NULL;
    fProfileCosineSixParticles[k] = NULL;
  }

} // End: InitialiseArraysOfTProfiles().

//======================================================================================//
// Methods called in UserCreateOutputObjects.                                           //
//--------------------------------------------------------------------------------------//
void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()";
  if (!fMainList) {Fatal(sMethodName.Data(), "Error: 'fMainList' is NULL.");}

// List with the observables used to select the events.
  fControlListEventCuts = new TList();
  fControlListEventCuts->SetName("ControlListEventCuts");
  fControlListEventCuts->SetOwner(kTRUE);
  fMainList->Add(fControlListEventCuts);

// List with the observables used to select the tracks.
  fControlListTrackCuts = new TList();
  fControlListTrackCuts->SetName("ControlListTrackCuts");
  fControlListTrackCuts->SetOwner(kTRUE);
  fMainList->Add(fControlListTrackCuts);

// List with the final multiTrack correlations, separated into sublists according to the number of particles involved in the correlations.
  fListCorrelations = new TList();
  fListCorrelations->SetName("ListCorrelations");
  fListCorrelations->SetOwner(kTRUE);
  fMainList->Add(fListCorrelations);

  fListTwoParticles = new TList();
  fListTwoParticles->SetName("ListTwoParticles");
  fListTwoParticles->SetOwner(kTRUE);
  fListCorrelations->Add(fListTwoParticles);

  fListThreeParticles = new TList();
  fListThreeParticles->SetName("ListThreeParticles");
  fListThreeParticles->SetOwner(kTRUE);
  fListCorrelations->Add(fListThreeParticles);

  fListFourParticles = new TList();
  fListFourParticles->SetName("ListFourParticles");
  fListFourParticles->SetOwner(kTRUE);
  fListCorrelations->Add(fListFourParticles);

  fListSixParticles = new TList();
  fListSixParticles->SetName("ListSixParticles");
  fListSixParticles->SetOwner(kTRUE);
  fListCorrelations->Add(fListSixParticles);

} // End: BookAllLists().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookControlListEventCuts()
{
/* Book the histograms with the observables for the event selection. */
// Centrality distribution.
  fHistoCentrality = new TH1D("fHistoCentrality", "Distribution of the centrality", 100, 0., 100.);
  fHistoCentrality->SetStats(kTRUE);
  fHistoCentrality->GetXaxis()->SetTitle("Centrality percentile");
  fControlListEventCuts->Add(fHistoCentrality);

// Multiplicity distribution.
  fHistoNumberOfTracksBeforeCuts = new TH1D("fHistoNumberOfTracksBeforeCuts", "Distribution of the number of tracks before any selection", 50000, 0., 50000.);
  fHistoNumberOfTracksBeforeCuts->SetStats(kTRUE);
  fHistoNumberOfTracksBeforeCuts->GetXaxis()->SetTitle("Number of tracks");
  fControlListEventCuts->Add(fHistoNumberOfTracksBeforeCuts);

  fHistoNumberOfTracksAfterEventCuts = new TH1D("fHistoNumberOfTracksAfterEventCuts", "Distribution of the number of tracks after the event selection", 50000, 0., 50000.);
  fHistoNumberOfTracksAfterEventCuts->SetStats(kTRUE);
  fHistoNumberOfTracksAfterEventCuts->GetXaxis()->SetTitle("Number of tracks");
  fControlListEventCuts->Add(fHistoNumberOfTracksAfterEventCuts);

// x-position of the PV distribution.
  fHistoVertexXBeforeCuts = new TH1D("fHistoVertexXBeforeCuts", "Distribution of PV_{x} before any selection", 1000, -20., 20.);
  fHistoVertexXBeforeCuts->SetStats(kTRUE);
  fHistoVertexXBeforeCuts->GetXaxis()->SetTitle("PV_{x}");
  fControlListEventCuts->Add(fHistoVertexXBeforeCuts);

  fHistoVertexXAfterEventCuts = new TH1D("fHistoVertexXAfterEventCuts", "Distribution of PV_{x} after the event selection", 1000, -20., 20.);
  fHistoVertexXAfterEventCuts->SetStats(kTRUE);
  fHistoVertexXAfterEventCuts->GetXaxis()->SetTitle("PV_{x}");
  fControlListEventCuts->Add(fHistoVertexXAfterEventCuts);

// y-position of the PV distribution.
  fHistoVertexYBeforeCuts = new TH1D("fHistoVertexYBeforeCuts", "Distribution of PV_{y} before any selection", 1000, -20., 20.);
  fHistoVertexYBeforeCuts->SetStats(kTRUE);
  fHistoVertexYBeforeCuts->GetXaxis()->SetTitle("PV_{y}");
  fControlListEventCuts->Add(fHistoVertexYBeforeCuts);

  fHistoVertexYAfterEventCuts = new TH1D("fHistoVertexYAfterEventCuts", "Distribution of PV_{y} after the event selection", 1000, -20., 20.);
  fHistoVertexYAfterEventCuts->SetStats(kTRUE);
  fHistoVertexYAfterEventCuts->GetXaxis()->SetTitle("PV_{y}");
  fControlListEventCuts->Add(fHistoVertexYAfterEventCuts);

// z-position of the PV distribution.
  fHistoVertexZBeforeCuts = new TH1D("fHistoVertexZBeforeCuts", "Distribution of PV_{z} before any selection", 1000, -20., 20.);
  fHistoVertexZBeforeCuts->SetStats(kTRUE);
  fHistoVertexZBeforeCuts->GetXaxis()->SetTitle("PV_{z}");
  fControlListEventCuts->Add(fHistoVertexZBeforeCuts);

  fHistoVertexZAfterEventCuts = new TH1D("fHistoVertexZAfterEventCuts", "Distribution of PV_{z} after the event selection", 1000, -20., 20.);
  fHistoVertexZAfterEventCuts->SetStats(kTRUE);
  fHistoVertexZAfterEventCuts->GetXaxis()->SetTitle("PV_{z}");
  fControlListEventCuts->Add(fHistoVertexZAfterEventCuts);

} // End: BookControlListEventCuts().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookControlListTrackCuts()
{
/* Book the histograms with the observables for the track selection. */
// Multiplicity distribution.
  fHistoNumberOfTracksAfterAllCuts = new TH1D("fHistoNumberOfTracksAfterAllCuts", "Distribution of the number of tracks after the track selection", 50000, 0., 50000.);
  fHistoNumberOfTracksAfterAllCuts->SetStats(kTRUE);
  fHistoNumberOfTracksAfterAllCuts->GetXaxis()->SetTitle("Number of tracks");
  fControlListTrackCuts->Add(fHistoNumberOfTracksAfterAllCuts);

// Transverse momentum distribution.
  fHistoPtBeforeCuts = new TH1D("fHistoPtBeforeCuts", "Distribution of p_{T} before the track selection", 1000, 0., 20.);
  fHistoPtBeforeCuts->SetStats(kTRUE);
  fHistoPtBeforeCuts->GetXaxis()->SetTitle("p_{T}");
  fControlListTrackCuts->Add(fHistoPtBeforeCuts);

  fHistoPtAfterCuts = new TH1D("fHistoPtAfterCuts", "Distribution of p_{T} after the track selection", 1000, 0., 20.);
  fHistoPtAfterCuts->SetStats(kTRUE);
  fHistoPtAfterCuts->GetXaxis()->SetTitle("p_{T}");
  fControlListTrackCuts->Add(fHistoPtAfterCuts);

// Pseudorapidity distribution.
  fHistoEtaBeforeCuts = new TH1D("fHistoEtaBeforeCuts", "Distribution of #eta before the track selection", 1000, -1., 1.);
  fHistoEtaBeforeCuts->SetStats(kTRUE);
  fHistoEtaBeforeCuts->GetXaxis()->SetTitle("#eta");
  fControlListTrackCuts->Add(fHistoEtaBeforeCuts);

  fHistoEtaAfterCuts = new TH1D("fHistoEtaAfterCuts", "Distribution of #eta after the track selection", 1000, -1., 1.);
  fHistoEtaAfterCuts->SetStats(kTRUE);
  fHistoEtaAfterCuts->GetXaxis()->SetTitle("#eta");
  fControlListTrackCuts->Add(fHistoEtaAfterCuts);

// Azimuthal angles distribution.
  fHistoPhiBeforeCuts = new TH1D("fHistoPhiBeforeCuts", "Distributiion of #phi before the track selection", 1000, 0., 6.3);
  fHistoPhiBeforeCuts->SetStats(kTRUE);
  fHistoPhiBeforeCuts->GetXaxis()->SetTitle("#phi");
  fControlListTrackCuts->Add(fHistoPhiBeforeCuts);

  fHistoPhiAfterCuts = new TH1D("fHistoPhiAfterCuts", "Distribution of #phi after the track selection", 1000, 0., 6.3);
  fHistoPhiAfterCuts->SetStats(kTRUE);
  fHistoPhiAfterCuts->GetXaxis()->SetTitle("#phi");
  fControlListTrackCuts->Add(fHistoPhiAfterCuts);

// TPC cluster distribution.
  fHistoTPCClustersBeforeCuts = new TH1D("fHistoTPCClustersBeforeCuts", "Distribution of the number of TPC clusters before the track selection", 1000, 0., 160.);
  fHistoTPCClustersBeforeCuts->SetStats(kTRUE);
  fHistoTPCClustersBeforeCuts->GetXaxis()->SetTitle("TPC clusters");
  fControlListTrackCuts->Add(fHistoTPCClustersBeforeCuts);

  fHistoTPCClustersAfterCuts = new TH1D("fHistoTPCClustersAfterCuts", "Distribution of the number of TPC clusters after the track selection", 1000, 0., 160.);
  fHistoTPCClustersAfterCuts->SetStats(kTRUE);
  fHistoTPCClustersAfterCuts->GetXaxis()->SetTitle("TPC clusters");
  fControlListTrackCuts->Add(fHistoTPCClustersAfterCuts);

// xy-plane of the DCA distribution.
  fHistoDCAXYBeforeCuts = new TH1D("fHistoDCAXYBeforeCuts", "Distribution of DCA_{xy} before the track selection", 1000, 0., 10.);
  fHistoDCAXYBeforeCuts->SetStats(kTRUE);
  fHistoDCAXYBeforeCuts->GetXaxis()->SetTitle("DCA_{xy}");
  fControlListTrackCuts->Add(fHistoDCAXYBeforeCuts);

  fHistoDCAXYAfterCuts = new TH1D("fHistoDCAXYAfterCuts", "Distribution of DCA_{xy} after the track selection", 1000, 0., 10.);
  fHistoDCAXYAfterCuts->SetStats(kTRUE);
  fHistoDCAXYAfterCuts->GetXaxis()->SetTitle("DCA_{xy}");
  fControlListTrackCuts->Add(fHistoDCAXYAfterCuts);

// z-coordinate of the DCA distribution.
  fHistoDCAZBeforeCuts = new TH1D("fHistoDCAZBeforeCuts", "Distribution of DCA_{z} before the track selection", 1000, 0., 10.);
  fHistoDCAZBeforeCuts->SetStats(kTRUE);
  fHistoDCAZBeforeCuts->GetXaxis()->SetTitle("DCA_{z}");
  fControlListTrackCuts->Add(fHistoDCAZBeforeCuts);

  fHistoDCAZAfterCuts = new TH1D("fHistoDCAZAfterCuts", "Distribution of DCA_{z} after the track selection", 1000, 0., 10.);
  fHistoDCAZAfterCuts->SetStats(kTRUE);
  fHistoDCAZAfterCuts->GetXaxis()->SetTitle("DCA_{z}");
  fControlListTrackCuts->Add(fHistoDCAZAfterCuts);

} // End: BookControlListTrackCuts().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::BookListCorrelations()
{
/* Book the TProfiles with the results for the multiTrack correlations. */
// For 2-particle correlations.
  for (Int_t t = 0; t < 6; t++)
  {
    fProfileCosineTwoParticles[t] = new TProfile("", "", 1, 0., 1.);
    fProfileCosineTwoParticles[t]->SetStats(kTRUE);
    fProfileCosineTwoParticles[t]->Sumw2();
    fProfileCosineTwoParticles[t]->SetName(Form("fProfileCosineTwoParticles_%d", t));
    fListTwoParticles->Add(fProfileCosineTwoParticles[t]);

    if (fComputeNestedLoops)
    {
      fProfileCosineTwoNestedLoops[t] = new TProfile("", "", 1, 0., 1.);
      fProfileCosineTwoNestedLoops[t]->SetStats(kTRUE);
      fProfileCosineTwoNestedLoops[t]->Sumw2();
      fProfileCosineTwoNestedLoops[t]->SetName(Form("fProfileCosineTwoNestedLoops_%d", t));
      fListTwoParticles->Add(fProfileCosineTwoNestedLoops[t]);
    }
  }

  for (Int_t c = 0; c < 5; c++)
  {
    fProfileTwoCosine[c] = new TProfile("", "", 1, 0., 1.);
    fProfileTwoCosine[c]->SetStats(kTRUE);
    fProfileTwoCosine[c]->Sumw2();
    fProfileTwoCosine[c]->SetName(Form("fProfileTwoCosine_%d", c));
    fListTwoParticles->Add(fProfileTwoCosine[c]);

    if (fComputeNestedLoops)
    {
      fProfileTwoCosineNestedLoops[c] = new TProfile("", "", 1, 0., 1.);
      fProfileTwoCosineNestedLoops[c]->SetStats(kTRUE);
      fProfileTwoCosineNestedLoops[c]->Sumw2();
      fProfileTwoCosineNestedLoops[c]->SetName(Form("fProfileTwoCosineNestedLoops_%d", c));
      fListTwoParticles->Add(fProfileTwoCosineNestedLoops[c]);
    }
  }

// For 3-particle correlations.
  for (Int_t th = 0; th < 4; th++)
  {
    fProfileCosineThreeParticles[th] = new TProfile("","", 1, 0., 1.);
    fProfileCosineThreeParticles[th]->SetStats(kTRUE);
    fProfileCosineThreeParticles[th]->Sumw2();
    fProfileCosineThreeParticles[th]->SetName(Form("fProfileCosineThreeParticles_%d", th));
    fListThreeParticles->Add(fProfileCosineThreeParticles[th]);

    if (fComputeNestedLoops)
    {
      fProfileCosineThreeNestedLoops[th] = new TProfile("","", 1, 0., 1.);
      fProfileCosineThreeNestedLoops[th]->SetStats(kTRUE);
      fProfileCosineThreeNestedLoops[th]->Sumw2();
      fProfileCosineThreeNestedLoops[th]->SetName(Form("fProfileCosineThreeNestedLoops_%d", th));
      fListThreeParticles->Add(fProfileCosineThreeNestedLoops[th]);
    }
  }

// For 4-particle correlations.
  for (Int_t f = 0; f < 6; f++)
  {
    fProfileCosineFourParticles[f] = new TProfile("", "", 1, 0., 1.);
    fProfileCosineFourParticles[f]->SetStats(kTRUE);
    fProfileCosineFourParticles[f]->Sumw2();
    fProfileCosineFourParticles[f]->SetName(Form("fProfileCosineFourParticles_%d", f));
    fListFourParticles->Add(fProfileCosineFourParticles[f]);

    if (fComputeNestedLoops)
    {
      fProfileCosineFourNestedLoops[f] = new TProfile("", "", 1, 0., 1.);
      fProfileCosineFourNestedLoops[f]->SetStats(kTRUE);
      fProfileCosineFourNestedLoops[f]->Sumw2();
      fProfileCosineFourNestedLoops[f]->SetName(Form("fProfileCosineFourNestedLoops_%d", f));
      fListFourParticles->Add(fProfileCosineFourNestedLoops[f]);
    }
  }

// For 6- and 8-particle correlations.
  for (Int_t s = 0; s < 4; s++)
  {
    fProfileCosineSixParticles[s] = new TProfile("", "", 1, 0., 1.);
    fProfileCosineSixParticles[s]->SetStats(kTRUE);
    fProfileCosineSixParticles[s]->Sumw2();
    fProfileCosineSixParticles[s]->SetName(Form("fProfileCosineSixParticles_%d", s));
    fListSixParticles->Add(fProfileCosineSixParticles[s]);
  }

  fProfileCosineEightParticles = new TProfile("", "", 1, 0., 1.);
  fProfileCosineEightParticles->SetStats(kTRUE);
  fProfileCosineEightParticles->Sumw2();
  fProfileCosineEightParticles->SetName("fProfileCosineEightParticles");
  fListSixParticles->Add(fProfileCosineEightParticles);

} // End: BookListCorrelations().

//======================================================================================//
// Methods called in UserExec(Option_t *).
//--------------------------------------------------------------------------------------//
void AliAnalysisTaskTwoMultiCorrelations::AODanalysis(AliAODEvent *aAODevent)
{
/* Compute the multiTrack correlations for an AOD event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AODanalysis(AliAODEvent *aAODevent)";

// Determine the detector to use for the centrality estimations.
  if (1 != (Int_t)fUseSPDForCentrality + (Int_t)fUseVZeroForCentrality)
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be selected in 'SetCentralityEstimation'.");
  }
  else if (fUseSPDForCentrality) {*fCentralitySelection = "CL1";}
  else if (fUseVZeroForCentrality) {*fCentralitySelection = "V0M";}

// Determine if the event belongs to this centrality range.
  AliMultSelection *ams = (AliMultSelection*)aAODevent->FindListObject("MultSelection");
  if (!ams) {return;} // Protection against NULL pointer.
  if (ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())) >= fCentralityMin && ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())) < fCentralityMax)
  {
    fHistoCentrality->Fill(ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())));
  }
  else {return;}  // This event does not belong to this centrality range.

// Fill the histogram for the distribution of the number of tracks before selection.
  long long numberOfTracksBeforeCuts = aAODevent->GetNumberOfTracks();  // Number of tracks before any selection.
  fHistoNumberOfTracksBeforeCuts->Fill(numberOfTracksBeforeCuts);

//......................................................................................//
/* Application of the event selection. */
// Cuts on the PV position.
  AliAODVertex *avtx = (AliAODVertex*)aAODevent->GetPrimaryVertex();  // Position of the PV.
  fHistoVertexXBeforeCuts->Fill(avtx->GetX());
  fHistoVertexYBeforeCuts->Fill(avtx->GetY());
  fHistoVertexZBeforeCuts->Fill(avtx->GetZ());

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

// Fill the histograms after the event selection.
  fHistoVertexXAfterEventCuts->Fill(avtx->GetX());
  fHistoVertexYAfterEventCuts->Fill(avtx->GetY());
  fHistoVertexZAfterEventCuts->Fill(avtx->GetZ());

//......................................................................................//
/* Application of the track selection. */
  long long numberOfTracksAfterEventCuts = aAODevent->GetNumberOfTracks();  // Number of tracks after the event selection and before the track selection.
  long long nParticlesAfterCuts = 0;  // Number of tracks after the track selection (= multiplicity).
  Int_t *hasTrackPassedTheCuts = new Int_t[numberOfTracksAfterEventCuts](); // Flag if a track passed the selection (1) or not (0).

  fHistoNumberOfTracksAfterEventCuts->Fill(numberOfTracksAfterEventCuts);

  for (long long iTrack = 0; iTrack < numberOfTracksAfterEventCuts; iTrack++)
  {
  // Get a pointer to the current track.
    AliAODTrack *currentTrack = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iTrack));
    if (!currentTrack) {continue;}  // Protection against NULL pointer.
    if (!currentTrack->TestFilterBit(128)){continue;} // Filter bit 128 denotes TPC-only tracks.

  // Get all the observables for the selection.
    Double_t ptBeforeCuts = currentTrack->Pt(); // Transverse momentum.    
    Double_t etaBeforeCuts = currentTrack->Eta(); // Pseudorapidity.
    Double_t phiBeforeCuts = currentTrack->Phi(); // Azimuthal angle.
    Int_t numberOfTPCclustersBeforeCuts = currentTrack->GetTPCncls(); // Number of TPC clusters.
    Double_t DCAxBeforeCuts = currentTrack->XAtDCA(); // x-coordinate of the DCA.
    Double_t DCAyBeforeCuts = currentTrack->YAtDCA(); // y-coordinate of the DCA.
    Double_t DCAzBeforeCuts = currentTrack->ZAtDCA(); // z-coordinate of the DCA.

    Double_t DCAxyBeforeCuts = TMath::Sqrt((DCAxBeforeCuts*DCAxBeforeCuts) + (DCAyBeforeCuts*DCAyBeforeCuts));

  // Fill the control histograms of some observables before the cuts.
    fHistoPtBeforeCuts->Fill(ptBeforeCuts);
    fHistoEtaBeforeCuts->Fill(etaBeforeCuts);
    fHistoPhiBeforeCuts->Fill(phiBeforeCuts);
    fHistoTPCClustersBeforeCuts->Fill(numberOfTPCclustersBeforeCuts);
    fHistoDCAXYBeforeCuts->Fill(DCAxyBeforeCuts);
    fHistoDCAZBeforeCuts->Fill(DCAzBeforeCuts);

  // Apply the track selection to the current track.
    if (CreateTrackSelection(ptBeforeCuts, etaBeforeCuts, numberOfTPCclustersBeforeCuts, DCAxyBeforeCuts, DCAzBeforeCuts))
    {
      nParticlesAfterCuts++;
      hasTrackPassedTheCuts[iTrack] = 1;
    }
    else {hasTrackPassedTheCuts[iTrack] = 0;}
  }

  fHistoNumberOfTracksAfterAllCuts->Fill(nParticlesAfterCuts);

// Define the observables after the track selection for the AOD analysis.
  Double_t *phi = new Double_t[nParticlesAfterCuts]();  // Azimuthal angle.
  Double_t *particleWeights = new Double_t[nParticlesAfterCuts]();  // Particle weight.
  Int_t indexNewArray = 0;  // New index of the selected particle after the selection.
  Double_t pt = 0.; // Transverse momentum.
  Double_t eta = 0.;  // Pseudorapidity.
  Int_t numberOfTPCclusters = 0;  // Number of TPC clusters.
  Double_t DCAx = 0.; // x-coordinate of the DCA.
  Double_t DCAy = 0.; // y-coordinate of the DCA.
  Double_t DCAz = 0.; // z-coordinate of the DCA.
  Double_t DCAxy = 0.;  // xy-plane of the DCA.

  if (nParticlesAfterCuts <= 2*fNumberHarmonicsInSC) {return;}  // Do the analysis only if the number of particles after the cuts is at least equal to the maximum number of particles needed for the studied GSC.

  for (long long iiTrack = 0; iiTrack < numberOfTracksAfterEventCuts; iiTrack++)
  {
  // Pointer to the current particle.
    AliAODTrack *keptParticle = dynamic_cast<AliAODTrack*>(aAODevent->GetTrack(iiTrack));
    if (!keptParticle) {continue;}  // Protection against NULL pointer.
    if(!keptParticle->TestFilterBit(128)) {continue;} // Filter bit 128 denotes TPC-only tracks.

    if (hasTrackPassedTheCuts[iiTrack] == 1)  // The particle passed the selection.
    {
    // Get the values of the observables for the current particle.
      phi[indexNewArray] = keptParticle->Phi();
      pt = keptParticle->Pt();
      eta = keptParticle->Eta();
      numberOfTPCclusters = keptParticle->GetTPCncls();
      DCAx = keptParticle->XAtDCA();
      DCAy = keptParticle->YAtDCA();
      DCAz = keptParticle->ZAtDCA();

      DCAxy = TMath::Sqrt((DCAx*DCAx) + (DCAy*DCAy));

      particleWeights[indexNewArray] = 1.;
      //if (fUseParticleWeights) {continue;}  // TBA
      //else {particleWeights[indexNewArray] = 1.;}

    // Fill the control histograms after the track selection.
      fHistoPhiAfterCuts->Fill(phi[indexNewArray]);
      fHistoPtAfterCuts->Fill(pt);
      fHistoEtaAfterCuts->Fill(eta);
      fHistoTPCClustersAfterCuts->Fill(numberOfTPCclusters);
      fHistoDCAXYAfterCuts->Fill(DCAxy);
      fHistoDCAZAfterCuts->Fill(DCAz);

      indexNewArray++;
    }
    else {continue;}
  }

// Calculate all the Q-vectors for the current event.
  CalculateQvectors(nParticlesAfterCuts, phi, particleWeights);

// Compute all the multiTrack correlations for the current event.
  GSCfullAnalysis(nParticlesAfterCuts, phi, particleWeights);

// Reset everything to zero for the next event.
  delete [] phi;
  delete [] particleWeights;
  indexNewArray = 0;
  pt = 0.;
  eta = 0.;
  numberOfTPCclusters = 0;
  DCAx = 0.;
  DCAy = 0.;
  DCAz = 0.;
  DCAxy = 0.;

} // End: AODanalysis().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::MCanalysis(AliMCEvent *aMCevent)
{
/* Compute the multiTrack correlations for a MC event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::MCanalysis(AliMCEvent *aMCevent)";

// Determine the detector to use for the centrality estimations.
  if (1 != (Int_t)fUseSPDForCentrality + (Int_t)fUseVZeroForCentrality)
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be selected in 'SetCentralityEstimation'.");
  }
  else if (fUseSPDForCentrality) {*fCentralitySelection = "CL1";}
  else if (fUseVZeroForCentrality) {*fCentralitySelection = "V0M";}

// Determine if the event belongs to this centrality range (for reconstructed particles only).
  if (!fProcessOnlyKine)
  {
    AliMultSelection *ams = (AliMultSelection*)aMCevent->FindListObject("MultSelection");
    if (!ams) {return;} // Protection against NULL pointer.
    if (ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())) >= fCentralityMin && ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())) < fCentralityMax)
    {
      fHistoCentrality->Fill(ams->GetMultiplicityPercentile(Form("%s", fCentralitySelection->Data())));
    }
    else {return;}  // This event does not belong to this centrality range.
  }

// Fill the histogram for the distribution of the number of tracks before selection.
  long long numberOfTracksBeforeCuts = aMCevent->GetNumberOfTracks(); // Number of tracks before any selection.
  fHistoNumberOfTracksBeforeCuts->Fill(numberOfTracksBeforeCuts);

//......................................................................................//
/* Application of the event selection. */
// Cuts on the PV position.
  AliMCVertex *avtx = (AliMCVertex*)aMCevent->GetPrimaryVertex();  // Position of the PV.
  fHistoVertexXBeforeCuts->Fill(avtx->GetX());
  fHistoVertexYBeforeCuts->Fill(avtx->GetY());
  fHistoVertexZBeforeCuts->Fill(avtx->GetZ());

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

// Fill the histograms after the event selection.
  fHistoVertexXAfterEventCuts->Fill(avtx->GetX());
  fHistoVertexYAfterEventCuts->Fill(avtx->GetY());
  fHistoVertexZAfterEventCuts->Fill(avtx->GetZ());

//......................................................................................//
/* Application of the track selection. */
  long long numberOfTracksAfterEventCuts = aMCevent->GetNumberOfTracks(); // Number of tracks after the event selection and before the track selection.
  long long nParticlesAfterCuts = 0;  // Number of tracks after the track selection (= multiplicity).
  Int_t *hasTrackPassedTheCuts = new Int_t[numberOfTracksAfterEventCuts](); // Flag if a track passed the selection (1) or not (0).

  fHistoNumberOfTracksAfterEventCuts->Fill(numberOfTracksAfterEventCuts);

  for (long long iTrack = 0; iTrack < numberOfTracksAfterEventCuts; iTrack++)
  {
  // Get a pointer to the current track.
    AliAODMCParticle *currentTrack = dynamic_cast<AliAODMCParticle*>(aMCevent->GetTrack(iTrack));
    if (!currentTrack) {continue;}  // Protection against NULL pointer.

  // Get all the observables for the selection.
    Double_t ptBeforeCuts = currentTrack->Pt(); // Transverse momentum.    
    Double_t etaBeforeCuts = currentTrack->Eta(); // Pseudorapidity.
    Double_t phiBeforeCuts = currentTrack->Phi(); // Azimuthal angle.

  // Fill the control histograms of some observables before the cuts.
    fHistoPtBeforeCuts->Fill(ptBeforeCuts);
    fHistoEtaBeforeCuts->Fill(etaBeforeCuts);
    fHistoPhiBeforeCuts->Fill(phiBeforeCuts);

  // Apply the track selection to the current track.
    if ((fPtMin <= ptBeforeCuts) && (ptBeforeCuts <= fPtMax) && (fEtaMin <= etaBeforeCuts) && (etaBeforeCuts <= fEtaMax))
    {
      nParticlesAfterCuts++;
      hasTrackPassedTheCuts[iTrack] = 1;
    }
    else {hasTrackPassedTheCuts[iTrack] = 0;}
  }

  fHistoNumberOfTracksAfterAllCuts->Fill(nParticlesAfterCuts);

// Define the observables after the track selection for the AOD analysis.
  Double_t *phi = new Double_t[nParticlesAfterCuts]();  // Azimuthal angle.
  Double_t *particleWeights = new Double_t[nParticlesAfterCuts]();  // Particle weight.
  Int_t indexNewArray = 0;  // New index of the selected particle after the selection.
  Double_t pt = 0.; // Transverse momentum.
  Double_t eta = 0.;  // Pseudorapidity.

  if (nParticlesAfterCuts <= 2*fNumberHarmonicsInSC) {return;}  // Do the analysis only if the number of particles after the cuts is at least equal to the maximum number of particles needed for the studied GSC.

  for (long long iiTrack = 0; iiTrack < numberOfTracksAfterEventCuts; iiTrack++)
  {
  // Pointer to the current particle.
    AliAODMCParticle *keptParticle = dynamic_cast<AliAODMCParticle*>(aMCevent->GetTrack(iiTrack));
    if (!keptParticle) {continue;}  // Protection against NULL pointer.

    if (hasTrackPassedTheCuts[iiTrack] == 1)  // The particle passed the selection.
    {
    // Get the values of the observables for the current particle.
      phi[indexNewArray] = keptParticle->Phi();
      pt = keptParticle->Pt();
      eta = keptParticle->Eta();

      particleWeights[indexNewArray] = 1.;
      //if (fUseParticleWeights) {continue;}  // TBA
      //else {particleWeights[indexNewArray] = 1.;}

    // Fill the control histograms after the track selection.
      fHistoPhiAfterCuts->Fill(phi[indexNewArray]);
      fHistoPtAfterCuts->Fill(pt);
      fHistoEtaAfterCuts->Fill(eta);

      indexNewArray++;
    }
    else {continue;}
  }

// Calculate all the Q-vectors for the current event.
  CalculateQvectors(nParticlesAfterCuts, phi, particleWeights);

// Compute all the multiTrack correlations for the current event.
  GSCfullAnalysis(nParticlesAfterCuts, phi, particleWeights);

// Reset everything to zero for the next event.
  delete [] phi;
  delete [] particleWeights;
  indexNewArray = 0;
  pt = 0.;
  eta = 0.;

} // End: MCDanalysis().

//======================================================================================//
Bool_t AliAnalysisTaskTwoMultiCorrelations::CreateTrackSelection(Double_t currentPt, Double_t currentEta, Int_t currentNumberOfTPC, Double_t currentDCAXY, Double_t currentDCAZ)
{
/* List the cuts to apply to the particles in a given event. */
  Bool_t cutOnPt = (fPtMin <= currentPt) && (currentPt <= fPtMax);
  Bool_t cutOnEta = (fEtaMin <= currentEta) && (currentEta <= fEtaMax);
  Bool_t cutOnNumberOfTPC = (fNumberOfTPCMin < currentNumberOfTPC);
  Bool_t cutOnDCAxy = (currentDCAXY < fDCAxyMax);
  Bool_t cutOnDCAz = (currentDCAZ < fDCAzMax);

  return cutOnPt && cutOnEta && cutOnNumberOfTPC && cutOnDCAxy && cutOnDCAz;

} // End: Bool_t CreateTrackSelection().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(long long nParticles, Double_t angles[], Double_t weights[])
{
/* Calculate all the Q-vectors for a given set of azimuthal angles and particles weights. */
  Double_t iAngle = 0.; // Angle of the current particle.
  Double_t iWeight = 1.;  // Particle weight of the current particle.
  Double_t weightToPowerP = 1.; // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    }
  }

// Compute the Q-vectors.
  for (long long iTrack = 0; iTrack < nParticles; iTrack++)
  {
    iAngle = angles[iTrack];
    if (fUseParticleWeights) {iWeight = weights[iTrack];}
    for (Int_t iHarmo = 0; iHarmo < 49; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < 9; iPower++)
      {
        if (fUseParticleWeights) {weightToPowerP = TMath::Power(iWeight, iPower);}
        fQvectors[iHarmo][iPower] += TComplex(weightToPowerP*TMath::Cos(iHarmo*iAngle), weightToPowerP*TMath::Sin(iHarmo*iAngle));
      }
    }
  }
} // End: CalculateQvectors().

//======================================================================================//
void AliAnalysisTaskTwoMultiCorrelations::GSCfullAnalysis(long long nParticles, Double_t angles[], Double_t weights[])
{
/* Compute all the needed multiTrack correlators (with possible cross-check from nested loops) for the chosen harmonics. */
  Int_t harmonicsArray[8] = {fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, fHarmonicFive, fHarmonicSix, fHarmonicSeven, fHarmonicEight}; // Harmonics (n_1,... n_8).

// Compute the denominator of each case of correlators.
  Int_t twoZerosArray[2] = {0};  // For 2p correlations.
  Double_t twoParticleDenominator = CalculateRecursion(2, twoZerosArray).Re();
  Int_t threeZerosArray[3] = {0}; // For 3p correlations.
  Double_t threeParticleDenominator = CalculateRecursion(3, threeZerosArray).Re();
  Int_t fourZerosArray[4] = {0}; // For 4p correlations.
  Double_t fourParticleDenominator = CalculateRecursion(4, fourZerosArray).Re();
  Int_t sixZerosArray[6] = {0};  // For 6p correlations.
  Double_t sixParticleDenominator = CalculateRecursion(6, sixZerosArray).Re();
  Int_t eightZerosArray[8] = {0}; // For 8p correlations.
  Double_t eightParticleDenominator = CalculateRecursion(8, eightZerosArray).Re();
  eightParticleDenominator++;  // Dummy line since this variable is not useful yet.

//**************************************************************************************//
// Compute the 2-particle correlations.
  Int_t firstHarmonicArray[2] = {harmonicsArray[0], harmonicsArray[1]}; // Array with (k,-k).
  Int_t secondHarmonicArray[2] = {harmonicsArray[2], harmonicsArray[3]};  // Array with (l,-l).
  Int_t thirdHarmonicArray[2] = {harmonicsArray[4], harmonicsArray[5]}; // Array with (m,-m).
  Int_t fourthHarmonicArray[2] = {harmonicsArray[6], harmonicsArray[7]};  // Array with (n,-n).

  Double_t tempoThirdHarmoNested = 0.;

/// 1st harmonic: <2>_{k,-k}.
  TComplex firstHarmonicComplex = (CalculateRecursion(2, firstHarmonicArray))/twoParticleDenominator; // Complex value of <2>_{k,-k}.
  Double_t firstHarmonicValue = firstHarmonicComplex.Re();  // Value of <2>_{k,-k}.
  fProfileCosineTwoParticles[0]->Fill(0.5, firstHarmonicValue, twoParticleDenominator);
  fProfileCosineTwoParticles[0]->SetTitle(Form("#LT2#GT_{k,-k}, k = %d", firstHarmonicArray[0]));    
  fProfileCosineTwoParticles[0]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", firstHarmonicArray[0], firstHarmonicArray[1]));

/// 2nd harmonic: <2>_{l,-l}.
  TComplex secondHarmonicComplex = (CalculateRecursion(2, secondHarmonicArray))/twoParticleDenominator; // Complex value of <2>_{l,-l}.
  Double_t secondHarmonicValue = secondHarmonicComplex.Re();  // Value of <2>_{l,-l}.
  fProfileCosineTwoParticles[1]->Fill(0.5, secondHarmonicValue, twoParticleDenominator);
  fProfileCosineTwoParticles[1]->SetTitle(Form("#LT2#GT_{l,-l}, l = %d", secondHarmonicArray[0]));
  fProfileCosineTwoParticles[1]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", secondHarmonicArray[0], secondHarmonicArray[1]));

/// 3rd harmonic: <2>_{m,-m}.
  TComplex thirdHarmonicComplex = TComplex(0.,0.);  // Complex value of <2>_{m,-m}.
  Double_t thirdHarmonicValue = 0.; // Value of <2>_{m,-m}.
  if (fNumberHarmonicsInSC >= 3)
  {
    thirdHarmonicComplex = (CalculateRecursion(2, thirdHarmonicArray))/twoParticleDenominator;
    thirdHarmonicValue = thirdHarmonicComplex.Re();
    fProfileCosineTwoParticles[2]->Fill(0.5, thirdHarmonicValue, twoParticleDenominator);
    fProfileCosineTwoParticles[2]->SetTitle(Form("#LT2#GT_{m,-m}, m = %d", thirdHarmonicArray[0]));
    fProfileCosineTwoParticles[2]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", thirdHarmonicArray[0], thirdHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC >= 3)
    
/// 4th harmonic: <2>_{n,-n}.
  TComplex fourthHarmonicComplex = TComplex(0.,0.); // Complex value of <2>_{n,-n}.
  Double_t fourthHarmonicValue = 0.;  // Value of <2>_{n,-n}.
  if (fNumberHarmonicsInSC == 4)
  {
    fourthHarmonicComplex = (CalculateRecursion(2, fourthHarmonicArray))/twoParticleDenominator;
    fourthHarmonicValue = fourthHarmonicComplex.Re();
    fProfileCosineTwoParticles[3]->Fill(0.5, fourthHarmonicValue, twoParticleDenominator);
    fProfileCosineTwoParticles[3]->SetTitle(Form("#LT2#GT_{n,-n}, n = %d", fourthHarmonicArray[0]));
    fProfileCosineTwoParticles[3]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", fourthHarmonicArray[0], fourthHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC == 4)

/// Sum first two harmonics: <2>_{k+l,-k-l}.
  Int_t sumTwoHarmonicsArray[2] = {firstHarmonicArray[0]+secondHarmonicArray[0], firstHarmonicArray[1]+secondHarmonicArray[1]}; // Array with (k+l,-k-l).
  TComplex sumTwoHarmonicsComplex = TComplex(0.,0.);  // Complex value of <2>_{k+l,-k-l}.
  Double_t sumTwoHarmonicsValue = 0.; // <2>_{k+l,-k-l}.
  if (fNumberHarmonicsInSC == 2)
  {
    sumTwoHarmonicsComplex = (CalculateRecursion(2, sumTwoHarmonicsArray))/twoParticleDenominator;
    sumTwoHarmonicsValue = sumTwoHarmonicsComplex.Re();
    fProfileCosineTwoParticles[4]->Fill(0.5, sumTwoHarmonicsValue, twoParticleDenominator);
    fProfileCosineTwoParticles[4]->SetTitle(Form("#LT2#GT_{k+l,-k-l}, k = %d, l = %d", firstHarmonicArray[0], secondHarmonicArray[0]));
    fProfileCosineTwoParticles[4]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", sumTwoHarmonicsArray[0], sumTwoHarmonicsArray[1]));
  } // End: if (fNumberHarmonicsInSC == 2)

/// Difference first two harmonics: <2>_{k-l,-k+l}.
  Int_t diffTwoHarmonicsArray[2] = {firstHarmonicArray[0]+secondHarmonicArray[1], firstHarmonicArray[1]+secondHarmonicArray[0]};  // Array with (k-l,-k+l).
  TComplex diffTwoHarmonicsComplex = TComplex(0.,0.); // Complex value of <2>_{k-l,-k+l}.
  Double_t diffTwoHarmonicsValue = 0.;  // <2>_{k-l,-k+l}.
  if (fNumberHarmonicsInSC == 2)
  {
    diffTwoHarmonicsComplex = (CalculateRecursion(2, diffTwoHarmonicsArray))/twoParticleDenominator;
    diffTwoHarmonicsValue = diffTwoHarmonicsComplex.Re();
    fProfileCosineTwoParticles[5]->Fill(0.5, diffTwoHarmonicsValue, twoParticleDenominator);
    fProfileCosineTwoParticles[5]->SetTitle(Form("#LT2#GT_{k-l,-k+l}, k = %d, l = %d", firstHarmonicArray[0], secondHarmonicArray[0]));
    fProfileCosineTwoParticles[5]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", diffTwoHarmonicsArray[0], diffTwoHarmonicsArray[1]));
  } // End: if (fNumberHarmonicsInSC == 2)

/// <2>_{k,-k}<2>_{l,-l}.
  Double_t doubleCosineValueFirstCase = firstHarmonicValue*secondHarmonicValue;
  fProfileTwoCosine[0]->Fill(0.5, doubleCosineValueFirstCase);
  fProfileTwoCosine[0]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{l,-l}, k = %d, l = %d", firstHarmonicArray[0], secondHarmonicArray[0]));
  fProfileTwoCosine[0]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], secondHarmonicArray[0], secondHarmonicArray[1]));

/// <2>_{k,-k}<2>_{m,-m}.
  Double_t doubleCosineValueSecondCase = 0.;
  if (fNumberHarmonicsInSC >= 3)
  {
    doubleCosineValueSecondCase = firstHarmonicValue*thirdHarmonicValue;
    fProfileTwoCosine[1]->Fill(0.5, doubleCosineValueSecondCase);
    fProfileTwoCosine[1]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{m,-m}, k = %d, m = %d", firstHarmonicArray[0], thirdHarmonicArray[0]));
    fProfileTwoCosine[1]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// <2>_{l,-l}<2>_{m,-m}.
  Double_t doubleCosineValueThirdCase = 0.;
  if (fNumberHarmonicsInSC >= 3)
  {
    doubleCosineValueThirdCase = secondHarmonicValue*thirdHarmonicValue;
    fProfileTwoCosine[2]->Fill(0.5, doubleCosineValueThirdCase);
    fProfileTwoCosine[2]->SetTitle(Form("#LT2#GT_{l,-l}#LT2#GT_{m,-m}, l = %d, m = %d", secondHarmonicArray[0], thirdHarmonicArray[0]));
    fProfileTwoCosine[2]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", secondHarmonicArray[0], secondHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// <2>_{k,-k}<2>_{l,-l}<2>_{m,-m}.
  Double_t doubleCosineValueFourthCase = 0.;
  if (fNumberHarmonicsInSC >= 3)
  {
    doubleCosineValueFourthCase = firstHarmonicValue*secondHarmonicValue*thirdHarmonicValue;
    fProfileTwoCosine[3]->Fill(0.5, doubleCosineValueFourthCase);
    fProfileTwoCosine[3]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{l,-l}#LT2#GT_{m,-m}, k = %d, l = %d, m = %d", firstHarmonicArray[0], secondHarmonicArray[0], thirdHarmonicArray[0]));
    fProfileTwoCosine[3]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], secondHarmonicArray[0], secondHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// Cross-check with nested loops for two-particle correlations.
  if (fComputeNestedLoops)
  {
  // 1st harmonic: <2>_{k,-k}.
    Double_t firstHarmonicValueNested = ComputeTwoNestedLoops(nParticles, firstHarmonicArray, angles, weights, fProfileCosineTwoNestedLoops[0]);
    fProfileCosineTwoNestedLoops[0]->SetTitle(Form("#LT2#GT_{k,-k}, k = %d", firstHarmonicArray[0]));
    fProfileCosineTwoNestedLoops[0]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", firstHarmonicArray[0], firstHarmonicArray[1]));

  // 2nd harmonic: <2>_{l,-l}.
    Double_t secondHarmonicValueNested = ComputeTwoNestedLoops(nParticles, secondHarmonicArray, angles, weights, fProfileCosineTwoNestedLoops[1]);
    fProfileCosineTwoNestedLoops[1]->SetTitle(Form("#LT2#GT_{l,-l}, l = %d", secondHarmonicArray[0]));
    fProfileCosineTwoNestedLoops[1]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", secondHarmonicArray[0], secondHarmonicArray[1]));

  // 3rd harmonic: <2>_{m,-m}.
    Double_t thirdHarmonicValueNested = 0.;
    if (fNumberHarmonicsInSC >= 3)
    {
      thirdHarmonicValueNested = ComputeTwoNestedLoops(nParticles, thirdHarmonicArray, angles, weights, fProfileCosineTwoNestedLoops[2]);
      fProfileCosineTwoNestedLoops[2]->SetTitle(Form("#LT2#GT_{m,-m}, m = %d", thirdHarmonicArray[0]));
      fProfileCosineTwoNestedLoops[2]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", thirdHarmonicArray[0], thirdHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC >= 3)
    tempoThirdHarmoNested = thirdHarmonicValueNested;

  // 4th harmonic: <2>_{n,-n}.
    Double_t fourthHarmonicValueNested = 0.;
    fourthHarmonicValueNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 4)
    {
      fourthHarmonicValueNested = ComputeTwoNestedLoops(nParticles, fourthHarmonicArray, angles, weights, fProfileCosineTwoNestedLoops[3]);
      fProfileCosineTwoNestedLoops[3]->SetTitle(Form("#LT2#GT_{n,-n}, n = %d", fourthHarmonicArray[0]));
      fProfileCosineTwoNestedLoops[3]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", fourthHarmonicArray[0], fourthHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC == 4)

  // Sum first two harmonics: <2>_{k+l,-k-l}.
    Double_t sumTwoHarmonicsValueNested = 0.;
    sumTwoHarmonicsValueNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      sumTwoHarmonicsValueNested = ComputeTwoNestedLoops(nParticles, sumTwoHarmonicsArray, angles, weights, fProfileCosineTwoNestedLoops[4]);
      fProfileCosineTwoNestedLoops[4]->SetTitle(Form("#LT2#GT_{n,-n}, n = %d", sumTwoHarmonicsArray[0]));
      fProfileCosineTwoNestedLoops[4]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", sumTwoHarmonicsArray[0], sumTwoHarmonicsArray[1]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // Difference first two harmonics: <2>_{k-l,-k+l}.
    Double_t diffTwoHarmonicsValueNested = 0.;
    diffTwoHarmonicsValueNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      diffTwoHarmonicsValueNested = ComputeTwoNestedLoops(nParticles, diffTwoHarmonicsArray, angles, weights, fProfileCosineTwoNestedLoops[5]);
      fProfileCosineTwoNestedLoops[5]->SetTitle(Form("#LT2#GT_{n,-n}, n = %d", diffTwoHarmonicsArray[0]));
      fProfileCosineTwoNestedLoops[5]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d, %d}", diffTwoHarmonicsArray[0], diffTwoHarmonicsArray[1]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // <<2>_{k,-k}<2>_{l,-l}>.
    Double_t doubleCosineFirstCaseNested = firstHarmonicValueNested*secondHarmonicValueNested;
    fProfileTwoCosineNestedLoops[0]->Fill(0.5, doubleCosineFirstCaseNested);
    fProfileTwoCosineNestedLoops[0]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{l,-l}, k = %d, l = %d", firstHarmonicArray[0], secondHarmonicArray[0]));
    fProfileTwoCosineNestedLoops[0]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], secondHarmonicArray[0], secondHarmonicArray[1]));

  // <2>_{k,-k}<2>_{m,-m}.
    Double_t doubleCosineSecondCaseNested = 0.;
    if (fNumberHarmonicsInSC >= 3)
    {
      doubleCosineSecondCaseNested = firstHarmonicValueNested*thirdHarmonicValueNested;
      fProfileTwoCosineNestedLoops[1]->Fill(0.5, doubleCosineSecondCaseNested);
      fProfileTwoCosineNestedLoops[1]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{m,-m}, k = %d, m = %d", firstHarmonicArray[0], thirdHarmonicArray[0]));
      fProfileTwoCosineNestedLoops[1]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC >= 3)

  // <2>_{l,-l}<2>_{m,-m}.
    Double_t doubleCosineThirdCaseNested = 0.;
    if (fNumberHarmonicsInSC >= 3)
    {
      doubleCosineThirdCaseNested = secondHarmonicValueNested*thirdHarmonicValueNested;
      fProfileTwoCosineNestedLoops[2]->Fill(0.5, doubleCosineThirdCaseNested);
      fProfileTwoCosineNestedLoops[2]->SetTitle(Form("#LT2#GT_{l,-l}#LT2#GT_{m,-m}, l = %d, m = %d", secondHarmonicArray[0], thirdHarmonicArray[0]));
      fProfileTwoCosineNestedLoops[2]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", secondHarmonicArray[0], secondHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC >= 3)

  // <2>_{k,-k}<2>_{l,-l}<2>_{m,-m}.
    Double_t doubleCosineFourthCaseNested = 0.;
    if (fNumberHarmonicsInSC >= 3)
    {
      doubleCosineFourthCaseNested = firstHarmonicValueNested*secondHarmonicValueNested*thirdHarmonicValueNested;
      fProfileTwoCosineNestedLoops[3]->Fill(0.5, doubleCosineFourthCaseNested);
      fProfileTwoCosineNestedLoops[3]->SetTitle(Form("#LT2#GT_{k,-k}#LT2#GT_{l,-l}#LT2#GT_{m,-m}, k = %d, l = %d, m = %d", firstHarmonicArray[0], secondHarmonicArray[0], thirdHarmonicArray[0]));
      fProfileTwoCosineNestedLoops[3]->GetYaxis()->SetTitle(Form("#LT#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#LT2#GT_{%d,%d}#GT", firstHarmonicArray[0], firstHarmonicArray[1], secondHarmonicArray[0], secondHarmonicArray[1], thirdHarmonicArray[0], thirdHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC >= 3)

  // Reset of the observables for the two nested loops.
    firstHarmonicValueNested = 0.;
    secondHarmonicValueNested = 0.;
    thirdHarmonicValueNested = 0.;
    fourthHarmonicValueNested = 0.;
    sumTwoHarmonicsValueNested = 0.;
    diffTwoHarmonicsValueNested = 0.;
    doubleCosineFirstCaseNested = 0.;
    doubleCosineSecondCaseNested = 0.;
    doubleCosineThirdCaseNested = 0.;
    doubleCosineFourthCaseNested = 0.;
  } // End: if (fComputeNestedLoops)

//**************************************************************************************//
// Compute the 3-particle correlations.
/// Case one: <3>_{k+l,-k,-l}.
  Int_t threeParticleCaseOneArray[3] = {harmonicsArray[0]+harmonicsArray[2], harmonicsArray[1], harmonicsArray[3]};
  TComplex threeParticleCaseOneComplex = TComplex(0.,0.); // Complex value of <3>_{k+l,-k,-l}.
  Double_t threeParticleCaseOneValue = 0.;  // Value of <3>_{k+l,-k,-l}.
  if (fNumberHarmonicsInSC == 2)
  {
    threeParticleCaseOneComplex = (CalculateRecursion(3, threeParticleCaseOneArray))/threeParticleDenominator;
    threeParticleCaseOneValue = threeParticleCaseOneComplex.Re();
    fProfileCosineThreeParticles[0]->Fill(0.5, threeParticleCaseOneValue, threeParticleDenominator);
    fProfileCosineThreeParticles[0]->SetTitle(Form("#LT3#GT_{k+l,-k,-l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
    fProfileCosineThreeParticles[0]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseOneArray[0], threeParticleCaseOneArray[1], threeParticleCaseOneArray[2]));
  } // End: if (fNumberHarmonicsInSC == 2)

/// Case two: <3>_{k-l,-k,+l}.
  Int_t threeParticleCaseTwoArray[3] = {harmonicsArray[0]+harmonicsArray[3], harmonicsArray[1], harmonicsArray[2]};
  TComplex threeParticleCaseTwoComplex = TComplex(0.,0.); // Complex value of <3>_{k+l,-k,-l}.
  Double_t threeParticleCaseTwoValue = 0.;  // Value of <3>_{k-l,-k,+l}.
  if (fNumberHarmonicsInSC == 2)
  {
    threeParticleCaseTwoComplex = (CalculateRecursion(3, threeParticleCaseTwoArray))/threeParticleDenominator;
    threeParticleCaseTwoValue = threeParticleCaseTwoComplex.Re();
    fProfileCosineThreeParticles[1]->Fill(0.5, threeParticleCaseTwoValue, threeParticleDenominator);
    fProfileCosineThreeParticles[1]->SetTitle(Form("#LT3#GT_{k-l,-k,l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
    fProfileCosineThreeParticles[1]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseTwoArray[0], threeParticleCaseTwoArray[1], threeParticleCaseTwoArray[2]));
  } // End: if (fNumberHarmonicsInSC == 2)

/// Case three: <3>_{k,l-k,-l}.
  Int_t threeParticleCaseThreeArray[3] = {harmonicsArray[0], harmonicsArray[2]+harmonicsArray[1], harmonicsArray[3]};
  TComplex threeParticleCaseThreeComplex = TComplex(0.,0.); // Complex value of <3>_{k,l-k,-l}.
  Double_t threeParticleCaseThreeValue = 0.;  // Value of <3>_{k,l-k,-l}.
  if (fNumberHarmonicsInSC == 2)
  {
    threeParticleCaseThreeComplex = (CalculateRecursion(3, threeParticleCaseThreeArray))/threeParticleDenominator;
    threeParticleCaseThreeValue = threeParticleCaseThreeComplex.Re();
    fProfileCosineThreeParticles[2]->Fill(0.5, threeParticleCaseThreeValue, threeParticleDenominator);
    fProfileCosineThreeParticles[2]->SetTitle(Form("#LT3#GT_{k,l-k,-l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
    fProfileCosineThreeParticles[2]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseThreeArray[0], threeParticleCaseThreeArray[1], threeParticleCaseThreeArray[2]));
  } // End: if (fNumberHarmonicsInSC == 2)

/// Case four: <3>_{k,-k-l,l}.
  Int_t threeParticleCaseFourArray[3] = {harmonicsArray[0], harmonicsArray[1]+harmonicsArray[3], harmonicsArray[2]};
  TComplex threeParticleCaseFourComplex = TComplex(0.,0.); // Complex value of <3>_{k,-k-l,l}.
  Double_t threeParticleCaseFourValue = 0.;  // Value of <3>_{k,-k-l,l}.
  if (fNumberHarmonicsInSC == 2)
  {
    threeParticleCaseFourComplex = (CalculateRecursion(3, threeParticleCaseFourArray))/threeParticleDenominator;
    threeParticleCaseFourValue = threeParticleCaseFourComplex.Re();
    fProfileCosineThreeParticles[3]->Fill(0.5, threeParticleCaseFourValue, threeParticleDenominator);
    fProfileCosineThreeParticles[3]->SetTitle(Form("#LT3#GT_{k,-k-l,l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
    fProfileCosineThreeParticles[3]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseFourArray[0], threeParticleCaseFourArray[1], threeParticleCaseFourArray[2]));
  } // End: if (fNumberHarmonicsInSC == 2)
    
/// Cross-check with nested loops for three-particle correlations.
  if (fComputeNestedLoops)
  {
  // Case one: <3>_{k+l,-k,-l}.
    Double_t threeParticleCaseOneNested = 0.;
    threeParticleCaseOneNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      threeParticleCaseOneNested = ComputeThreeNestedLoops(nParticles, threeParticleCaseOneArray, angles, weights, fProfileCosineThreeNestedLoops[0]);
      fProfileCosineThreeNestedLoops[0]->SetTitle(Form("#LT3#GT_{k+l,-k,-l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
      fProfileCosineThreeNestedLoops[0]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseOneArray[0], threeParticleCaseOneArray[1], threeParticleCaseOneArray[2]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // Case two: <3>_{k-l,-k,+l}.
    Double_t threeParticleCaseTwoNested = 0.;
    threeParticleCaseTwoNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      threeParticleCaseTwoNested = ComputeThreeNestedLoops(nParticles, threeParticleCaseTwoArray, angles, weights, fProfileCosineThreeNestedLoops[1]);
      fProfileCosineThreeNestedLoops[1]->SetTitle(Form("#LT3#GT_{k-l,-k,l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
      fProfileCosineThreeNestedLoops[1]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseTwoArray[0], threeParticleCaseTwoArray[1], threeParticleCaseTwoArray[2]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // Case three: <3>_{k,l-k,-l}.
    Double_t threeParticleCaseThreeNested = 0.;
    threeParticleCaseThreeNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      threeParticleCaseThreeNested = ComputeThreeNestedLoops(nParticles, threeParticleCaseThreeArray, angles, weights, fProfileCosineThreeNestedLoops[2]);
      fProfileCosineThreeNestedLoops[2]->SetTitle(Form("#LT3#GT_{k,l-k,-l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
      fProfileCosineThreeNestedLoops[2]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseThreeArray[0], threeParticleCaseThreeArray[1], threeParticleCaseThreeArray[2]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // Case four: <3>_{k,-k-l,l}.
    Double_t threeParticleCaseFourNested = 0.;
    threeParticleCaseFourNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC == 2)
    {
      threeParticleCaseFourNested = ComputeThreeNestedLoops(nParticles, threeParticleCaseFourArray, angles, weights, fProfileCosineThreeNestedLoops[3]);
      fProfileCosineThreeNestedLoops[3]->SetTitle(Form("#LT3#GT_{k,-k-l,l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
      fProfileCosineThreeNestedLoops[3]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", threeParticleCaseFourArray[0], threeParticleCaseFourArray[1], threeParticleCaseFourArray[2]));
    } // End: if (fNumberHarmonicsInSC == 2)

  // Reset of the observables for the three nested loops.
    threeParticleCaseOneNested = 0.;
    threeParticleCaseTwoNested = 0.;
    threeParticleCaseThreeNested = 0.;
    threeParticleCaseFourNested = 0.;
  } // End: if (fComputeNestedLoops)

//**************************************************************************************//
// Compute the 4-particle correlations.
/// Case one: <4>_{k,-k,l,-l}.
  Int_t fourParticleCaseOneArray[4] = {harmonicsArray[0], harmonicsArray[1], harmonicsArray[2], harmonicsArray[3]};
  TComplex fourParticleCaseOneComplex = (CalculateRecursion(4, fourParticleCaseOneArray))/fourParticleDenominator; // Complex value of <4>_{k,-k,l,-l}.
  Double_t fourParticleCaseOneValue = fourParticleCaseOneComplex.Re();  // Value of <4>_{k,-k,l,-l}.
  fProfileCosineFourParticles[0]->Fill(0.5, fourParticleCaseOneValue, fourParticleDenominator);
  fProfileCosineFourParticles[0]->SetTitle(Form("#LT4#GT_{k,-k,l,-l}, k = %d, l = %d", fourParticleCaseOneArray[0], fourParticleCaseOneArray[2]));
  fProfileCosineFourParticles[0]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseOneArray[0], fourParticleCaseOneArray[1],fourParticleCaseOneArray[2],fourParticleCaseOneArray[3]));

/// Case two: <4>_{k,-k,m,-m}.
  Int_t fourParticleCaseTwoArray[4] = {harmonicsArray[0], harmonicsArray[1], harmonicsArray[4], harmonicsArray[5]};
  TComplex fourParticleCaseTwoComplex = TComplex(0.,0.); // Complex value of <4>_{k,-k,m,-m}.
  Double_t fourParticleCaseTwoValue = 0.;  // Value of <4>_{k,-k,m,-m}.
  if (fNumberHarmonicsInSC >= 3)
  {
    fourParticleCaseTwoComplex = (CalculateRecursion(4, fourParticleCaseTwoArray))/fourParticleDenominator;
    fourParticleCaseTwoValue = fourParticleCaseTwoComplex.Re();
    fProfileCosineFourParticles[1]->Fill(0.5, fourParticleCaseTwoValue, fourParticleDenominator);
    fProfileCosineFourParticles[1]->SetTitle(Form("#LT4#GT_{k,-k,m,-m}, k = %d, m = %d", fourParticleCaseTwoArray[0], fourParticleCaseTwoArray[2]));
    fProfileCosineFourParticles[1]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseTwoArray[0], fourParticleCaseTwoArray[1],fourParticleCaseTwoArray[2],fourParticleCaseTwoArray[3]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// Case three: <4>_{l,-l,m,-m}.
  Int_t fourParticleCaseThreeArray[4] = {harmonicsArray[2], harmonicsArray[3], harmonicsArray[4], harmonicsArray[5]};
  TComplex fourParticleCaseThreeComplex = TComplex(0.,0.);  // Complex value of <4>_{l,-l,m,-m}.
  Double_t fourParticleCaseThreeValue = 0.; // Value of <4>_{l,-l,m,-m}.
  if (fNumberHarmonicsInSC >= 3)
  {
    fourParticleCaseThreeComplex = (CalculateRecursion(4, fourParticleCaseThreeArray))/fourParticleDenominator;
    fourParticleCaseThreeValue = fourParticleCaseThreeComplex.Re();
    fProfileCosineFourParticles[2]->Fill(0.5, fourParticleCaseThreeValue, fourParticleDenominator);
    fProfileCosineFourParticles[2]->SetTitle(Form("#LT4#GT_{l,-l,m,-m}, l = %d, m = %d", fourParticleCaseThreeArray[0], fourParticleCaseThreeArray[2]));
    fProfileCosineFourParticles[2]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseThreeArray[0], fourParticleCaseThreeArray[1],fourParticleCaseThreeArray[2],fourParticleCaseThreeArray[3]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// TBI: cases with four different harmonics...

/// <4>_{k,-k,l,-l}<2>_{m,-m}.
  Double_t doubleCosineValueFifthCase = 0.;
  if (fNumberHarmonicsInSC >= 3)
  {
    doubleCosineValueFifthCase = fourParticleCaseOneValue*thirdHarmonicValue;
    fProfileTwoCosine[4]->Fill(0.5, doubleCosineValueFifthCase);
    fProfileTwoCosine[4]->SetTitle(Form("#LT4#GT_{k,-k,l,-l}#LT2#GT_{m,-m}, k = %d, l = %d, m = %d", fourParticleCaseOneArray[0], fourParticleCaseOneArray[2], thirdHarmonicArray[0]));
    fProfileTwoCosine[4]->GetYaxis()->SetTitle(Form("#LT#LT4#GT_{%d,%d,%d,%d}#LT2#GT_{%d,%d}#GT", fourParticleCaseOneArray[0], fourParticleCaseOneArray[1], fourParticleCaseOneArray[2], fourParticleCaseOneArray[3], thirdHarmonicArray[0], thirdHarmonicArray[1]));
  } // End: if (fNumberHarmonicsInSC >= 3)

/// Cross-check with nested loops for four-particle correlations.
  if (fComputeNestedLoops)
  {
  // Case one: <4>_{k,-k,l,-l}.
    Double_t fourParticleCaseOneNested = ComputeFourNestedLoops(nParticles, fourParticleCaseOneArray, angles, weights, fProfileCosineFourNestedLoops[0]);
    fProfileCosineFourNestedLoops[0]->SetTitle(Form("#LT4#GT_{k,-k,l,-l}, k = %d, l = %d", harmonicsArray[0], harmonicsArray[2]));
    fProfileCosineFourNestedLoops[0]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseOneArray[0], fourParticleCaseOneArray[1], fourParticleCaseOneArray[2],fourParticleCaseOneArray[3]));

  // Case two: <4>_{k,m,-k,-m}.
    Double_t fourParticleCaseTwoNested = 0.;
    fourParticleCaseTwoNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC >= 3)
    {
      fourParticleCaseTwoNested = ComputeFourNestedLoops(nParticles, fourParticleCaseTwoArray, angles, weights, fProfileCosineFourNestedLoops[1]);
      fProfileCosineFourNestedLoops[1]->SetTitle(Form("#LT4#GT_{k,-k,m,-m}, k = %d, m = %d", harmonicsArray[0], harmonicsArray[4]));
      fProfileCosineFourNestedLoops[1]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseTwoArray[0], fourParticleCaseTwoArray[1], fourParticleCaseTwoArray[2],fourParticleCaseTwoArray[3]));
    } // End: if (fNumberHarmonicsInSC >= 3) 

  // Case three: <4>_{l,m,-l,-m}.
    Double_t fourParticleCaseThreeNested = 0.;
    fourParticleCaseThreeNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC >= 3)
    {
      fourParticleCaseThreeNested = ComputeFourNestedLoops(nParticles, fourParticleCaseThreeArray, angles, weights, fProfileCosineFourNestedLoops[2]);
      fProfileCosineFourNestedLoops[2]->SetTitle(Form("#LT4#GT_{l,-l,m,-m}, l = %d, m = %d", harmonicsArray[2], harmonicsArray[4]));
      fProfileCosineFourNestedLoops[2]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", fourParticleCaseThreeArray[0], fourParticleCaseThreeArray[1], fourParticleCaseThreeArray[2],fourParticleCaseThreeArray[3]));
    } // End: if (fNumberHarmonicsInSC >= 3)

  // TBI: cases with four harmonics...

  // <4>_{k,-k,l,-l}<2>_{m,-m}.
    Double_t doubleCosineFifthCaseNested = 0.;
    doubleCosineFifthCaseNested++;  // Dummy line since this variable is not useful yet.
    if (fNumberHarmonicsInSC >= 3)
    {
      doubleCosineFifthCaseNested = fourParticleCaseOneNested*tempoThirdHarmoNested;
      fProfileTwoCosineNestedLoops[4]->Fill(0.5, doubleCosineFifthCaseNested);
      fProfileTwoCosineNestedLoops[4]->SetTitle(Form("#LT4#GT_{k,-k,l,-l}#LT2#GT_{m,-m}, k = %d, l = %d, m = %d", fourParticleCaseOneArray[0], fourParticleCaseOneArray[2], thirdHarmonicArray[0]));
      fProfileTwoCosineNestedLoops[4]->GetYaxis()->SetTitle(Form("#LT#LT4#GT_{%d,%d,%d,%d}#LT2#GT_{%d,%d}#GT", fourParticleCaseOneArray[0], fourParticleCaseOneArray[1], fourParticleCaseOneArray[2], fourParticleCaseOneArray[3], thirdHarmonicArray[0], thirdHarmonicArray[1]));
    } // End: if (fNumberHarmonicsInSC >= 3)
      
  // Reset of the observables for the four nested loops.
    fourParticleCaseOneNested = 0.;
    fourParticleCaseTwoNested = 0.;
    fourParticleCaseThreeNested = 0.;
    //
    doubleCosineFifthCaseNested = 0.;
  } // End: if (fComputeNestedLoops)

//**************************************************************************************//
// Compute the 6-particle correlations.
/// Case one: <6>_{k,-k,l,-l,m,-m}.
  Int_t sixParticleCaseOneArray[6] = {harmonicsArray[0], harmonicsArray[1], harmonicsArray[2], harmonicsArray[3], harmonicsArray[4], harmonicsArray[5]};
  TComplex sixParticleCaseOneComplex = TComplex(0.,0.); // Complex value of <6>_{k,-k,l,-l,m,-m}.
  Double_t sixParticleCaseOneValue = 0.;  // Value of <6>_{k,-k,l,-l,m,-m}.
  if (fNumberHarmonicsInSC >= 3)
  {
    sixParticleCaseOneComplex = (CalculateRecursion(6, sixParticleCaseOneArray))/sixParticleDenominator;
    sixParticleCaseOneValue = sixParticleCaseOneComplex.Re();
    fProfileCosineSixParticles[0]->Fill(0.5, sixParticleCaseOneValue, sixParticleDenominator);
    fProfileCosineSixParticles[0]->SetTitle(Form("#LT6#GT_{k,-k,l,-l,m,-m}, k = %d, l = %d, m = %d", sixParticleCaseOneArray[0], sixParticleCaseOneArray[2], sixParticleCaseOneArray[4]));
    fProfileCosineSixParticles[0]->GetYaxis()->SetTitle(Form("#LT#LT6#GT#GT_{%d,%d,%d,%d,%d,%d}", sixParticleCaseOneArray[0], sixParticleCaseOneArray[1], sixParticleCaseOneArray[2], sixParticleCaseOneArray[3], sixParticleCaseOneArray[4], sixParticleCaseOneArray[5]));
  } // End: if (fNumberHarmonicsInSC >= 3)

//**************************************************************************************//
// Compute the 8-particle correlations.
/// TBI

//**************************************************************************************//
// Reset of the variables to zero.
  twoParticleDenominator = 0.;
  threeParticleDenominator = 0.;
  fourParticleDenominator = 0.;
  sixParticleDenominator = 0.;
  eightParticleDenominator = 0.;
  tempoThirdHarmoNested = 0.;

  firstHarmonicComplex = TComplex(0.,0.);
  firstHarmonicValue = 0.;
  secondHarmonicComplex = TComplex(0.,0.);
  secondHarmonicValue = 0.;
  thirdHarmonicComplex = TComplex(0.,0.);
  thirdHarmonicValue = 0.;
  fourthHarmonicComplex = TComplex(0.,0.);
  fourthHarmonicValue = 0.;

  doubleCosineValueFirstCase = 0.,
  doubleCosineValueSecondCase = 0.;
  doubleCosineValueThirdCase = 0.;

  threeParticleCaseOneComplex = TComplex(0.,0.);
  threeParticleCaseOneValue = 0.;
  threeParticleCaseTwoComplex = TComplex(0.,0.);
  threeParticleCaseTwoValue = 0.;
  threeParticleCaseThreeComplex = TComplex(0.,0.);
  threeParticleCaseThreeValue = 0.;
  threeParticleCaseFourComplex = TComplex(0.,0.);
  threeParticleCaseFourValue = 0.;

  fourParticleCaseOneComplex = TComplex(0.,0.);
  fourParticleCaseOneValue = 0.;
  fourParticleCaseTwoComplex = TComplex(0.,0.);
  fourParticleCaseTwoValue = 0.;
  fourParticleCaseThreeComplex = TComplex(0.,0.);
  fourParticleCaseThreeValue = 0.;
  // TBI: cases with four possible harmonics: kn, ln, mn...
  doubleCosineValueFifthCase = 0.;

  sixParticleCaseOneComplex = TComplex(0.,0.);
  sixParticleCaseOneValue = 0.;
} // End: GSCfullAnalysis().

//======================================================================================//
TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
/* Simplify the use of the Q-vectors. */
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]);
} // End: Q().

//======================================================================================//
TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
/* Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). */
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
} // End: CalculateRecursion().

//======================================================================================//
Double_t AliAnalysisTaskTwoMultiCorrelations::ComputeTwoNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile)
{
/* Calculate the two-particle correlations using two nested loops. */
  Double_t twoParticleCosine = 0.;  // cos(k(phi_1 - phi_2)).
  Double_t twoParticleWeight = 1.;  // Total particle weights.
  TProfile twoProfile;  // Returns <cos(k(phi1-phi2))>.
  twoProfile.Sumw2();

  for (long long m = 0; m < nParticles; m++)
  {
    for (long long n = 0; n < nParticles; n++)
    {
    // Remove the autocorrelations m==n.
      if (m == n) {continue;}
    // Compute cos(k(phi_1-phi_2)).
      twoParticleCosine = TMath::Cos(harmonic[0]*angles[m] + harmonic[1]*angles[n]);
      twoParticleWeight = weights[m] + weights[n];
      profile->Fill(0.5, twoParticleCosine, twoParticleWeight);
      twoProfile.Fill(0.5, twoParticleCosine, twoParticleWeight);

    // Reset the variables to default.
      twoParticleCosine = 0.;
      twoParticleWeight = 1.;
    } // End: for (long long n = 0; n < nParticles; n++)
  } // End: for (long long m = 0; m < nParticles; m++)

// Return <cos(k(phi_1-phi_2))>.
  return twoProfile.GetBinContent(1);
} // End: ComputeTwoNestedLoops(long long, Int_t*, Double_t[], Double_t[], TProfile*)

//======================================================================================//
Double_t AliAnalysisTaskTwoMultiCorrelations::ComputeThreeNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile)
{
/* Calculate the three-particle correlations using three nested loops. */
  Double_t threeParticleCosine = 0.;  // cos(kphi_1 +lphi_2 +mphi_3).
  Double_t threeParticleWeight = 1.;  // Total particle weights.
  TProfile threeProfile;  // Returns <cos(kphi_1 +lphi_2 +mphi_3)>.
  threeProfile.Sumw2();

  for (long long k = 0; k < nParticles; k++)
  {
    for (long long l = 0; l < nParticles; l++)
    {
      if (k == l) {continue;}

      for (long long m = 0; m < nParticles; m++)
      {
        if ((k == m) || (l == m)) {continue;}

      // Compute cos(kphi_1 +lphi_2 +mphi_3).
        threeParticleCosine = TMath::Cos(harmonic[0]*angles[k] + harmonic[1]*angles[l] + harmonic[2]*angles[m]);
        threeParticleWeight = weights[k] + weights[l] + weights[m];
        profile->Fill(0.5, threeParticleCosine, threeParticleWeight);
        threeProfile.Fill(0.5, threeParticleCosine, threeParticleWeight);

      // Reset the variables to default.
        threeParticleCosine = 0.;
        threeParticleWeight = 1.;
      } // End: for (long long m = 0; m < nParticles; m++)
    } // End: for (long long l = 0; l < nParticles; l++)
} // End: for (long long k = 0; k < nParticles; k++)

  // Return <cos(kphi_1 +lphi_2 +mphi_3)>.
  return threeProfile.GetBinContent(1);
} // End: ComputeThreeNestedLoops().

//======================================================================================//
Double_t AliAnalysisTaskTwoMultiCorrelations::ComputeFourNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile)
{
// Calculate the four-particle correlations using four nested loops.
  Double_t fourParticleCosine = 0.; // cos(kphi_1 +lphi_2 +mphi_3 +nphi_4).
  Double_t fourParticleWeight = 1.; // Total particle weights.
  TProfile fourProfile; // Returns <cos(kphi_1 +lphi_2 +mphi_3 +nphi_4)>.
  fourProfile.Sumw2();

  for (long long k = 0; k < nParticles; k++)
  {
    for (long long l = 0; l < nParticles; l++)
    {
      if (k == l) {continue;}
      for (long long m = 0; m < nParticles; m++)
      {
        if ((k == m) || (l == m)) {continue;}
        for (long long n = 0; n < nParticles; n++)
        {
          if ((k == n) || (l == n) || (m == n)) {continue;}
        // Computate cos(k(phi_1 - phi_2) + l(phi_3-phi_4)).
          fourParticleCosine = TMath::Cos(harmonic[0]*angles[k] + harmonic[1]*angles[l] + harmonic[2]*angles[m] + harmonic[3]*angles[n]);
          fourParticleWeight = weights[k] + weights[l] + weights[m] + weights[n];
          profile->Fill(0.5, fourParticleCosine, fourParticleWeight);
          fourProfile.Fill(0.5, fourParticleCosine, fourParticleWeight);

        // Reset the variables to default.
          fourParticleCosine = 0.;
          fourParticleWeight = 1.;
        } // End: for (long long n = 0; n < nParticles; n++)
      } // End: for (long long m = 0; m < nParticles; m++)
    } // End: for (long long l = 0; l < nParticles; l++)
  } // End: for (long long k = 0; k < nParticles; k++)

  // Return <cos(kphi_1 +lphi_2 +mphi_3 +nphi_4)>.
  return fourProfile.GetBinContent(1);
} // End: Double_t ComputeFourNestedLoops().

//======================================================================================//
// Methods called in Terminate(Option_t *).
//--------------------------------------------------------------------------------------//

//======================================================================================//

