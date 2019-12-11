
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

#include "AliAnalysisTaskTwoMultiCorrelations.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "Riostream.h"
#include "TList.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
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
#include <vector>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskTwoMultiCorrelations)

/* ========================================================================== /
/ Mandatory methods needed for AliAnalysisTaskSE.                             /
/ ========================================================================== */
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations() :
  AliAnalysisTaskSE(),
// General parameters for the analysis.
  fMainList(NULL),
  fHighestHarmonic(8),
  fLargestCorrelators(8),
  fComputeEtaGaps(kFALSE),
  fDoCorrelationsHisto(kFALSE),
  fUseNonUnitParticleWeights(kFALSE),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
// Type of files used in the analysis.
  fAnalyseOnlyAOD(kFALSE),
  fAnalyseOnlyMC(kFALSE),
  fAnalyseBothMCandAOD(kFALSE),
  fAODevent(NULL),
  fMCevent(NULL),
// Parameters and histograms related to the centrality.
  fMultiplicityList(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE), 
  fCentralityMin(0.),
  fCentralityMax(100.),
  fMultSelection(NULL),
// Parameters and histograms related to the number of tracks.
  fInitialMultiplicity(0),
// Parameters and histograms related to the event selection.
  fEventSelectionList(NULL),
  fCutOnPVx(kFALSE),
  fPVxMin(-44.),
  fPVxMax(44.),
  fCutOnPVy(kFALSE),
  fPVyMin(-44.),
  fPVyMax(44.),
  fCutOnPVz(kFALSE),
  fPVzMin(-10.),
  fPVzMax(-10.),
  fMultiplicityMin(6),
  fMainFilter(128),
  fGlobalFilter(256),
  fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.),
  fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.),
  fMultiplicityMaxB(0.),
// Parameters and histograms related to the track selection.
  fTrackSelectionList(NULL),
  fFilterbitIndex(0),
  fCutOnPt(kFALSE),
  fPtMin(0.2),
  fPtMax(5.),
  fHistoEfficiency(NULL),
  fHistoEffInverse(NULL),
  fCutOnEta(kFALSE),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fCutOnNTPC(kFALSE),
  fNTPCMin(70),
  fCutOnChi(kFALSE),
  fChiMin(0.1),
  fChiMax(4.),
  fCutOnNITS(kFALSE),
  fNITSMin(2),
  fCutOnDCAxy(kFALSE),
  fDCAxyMax(3.2),
  fCutOnDCAz(kFALSE),
  fDCAzMax(2.4),
  fCutOnCharge(kFALSE),
  fKeepPositiveCharges(kFALSE),
// Parameters related to the multi-particle correlations.
  fMPCList(NULL),
  fReducedQPower(0),
  fProfileTwoPartCorrel(NULL),
  fProfileFourPartCorrel(NULL),
  fProfileFourPartCorrelCheck(NULL),
  fProfileSixPartCorrel(NULL),
// Parameters related to the 2-particle correlations with eta gaps.
  fTPCEtaList(NULL)
{
/* Dummy constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations()");

// Initialise all elements in the arrays.
  InitialiseArraysOfQvectors();
  InitialiseArraysOfHistos();
  InitialiseArraysOfTProfiles();

}   // End of the dummy constructor.

/* ------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
  AliAnalysisTaskSE(name),
// General parameters for the analysis.
  fMainList(NULL),
  fHighestHarmonic(8),
  fLargestCorrelators(8),
  fComputeEtaGaps(kFALSE),
  fDoCorrelationsHisto(kFALSE),
  fUseNonUnitParticleWeights(kFALSE),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
// Type of files used in the analysis.
  fAnalyseOnlyAOD(kFALSE),
  fAnalyseOnlyMC(kFALSE),
  fAnalyseBothMCandAOD(kFALSE),
  fAODevent(NULL),
  fMCevent(NULL),
// Parameters and histograms related to the centrality.
  fMultiplicityList(NULL),
  fCentralityFromVZero(kFALSE),
  fCentralityFromSPD(kFALSE), 
  fCentralityMin(0.),
  fCentralityMax(100.),
  fMultSelection(NULL),
// Parameters and histograms related to the number of tracks.
  fInitialMultiplicity(0),
// Parameters and histograms related to the event selection.
  fEventSelectionList(NULL),
  fCutOnPVx(kFALSE),
  fPVxMin(-44.),
  fPVxMax(44.),
  fCutOnPVy(kFALSE),
  fPVyMin(-44.),
  fPVyMax(44.),
  fCutOnPVz(kFALSE),
  fPVzMin(-10.),
  fPVzMax(-10.),
  fMultiplicityMin(6),
  fMainFilter(128),
  fGlobalFilter(256),
  fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.),
  fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.),
  fMultiplicityMaxB(0.),
// Parameters and histograms related to the track selection.
  fTrackSelectionList(NULL),
  fFilterbitIndex(0),
  fCutOnPt(kFALSE),
  fPtMin(0.2),
  fPtMax(5.),
  fHistoEfficiency(NULL),
  fHistoEffInverse(NULL),
  fCutOnEta(kFALSE),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fCutOnNTPC(kFALSE),
  fNTPCMin(70),
  fCutOnChi(kFALSE),
  fChiMin(0.1),
  fChiMax(4.),
  fCutOnNITS(kFALSE),
  fNITSMin(2),
  fCutOnDCAxy(kFALSE),
  fDCAxyMax(3.2),
  fCutOnDCAz(kFALSE),
  fDCAzMax(2.4),
  fCutOnCharge(kFALSE),
  fKeepPositiveCharges(kFALSE),
// Parameters related to the multi-particle correlations.
  fMPCList(NULL),
  fReducedQPower(0),
  fProfileTwoPartCorrel(NULL),
  fProfileFourPartCorrel(NULL),
  fProfileFourPartCorrelCheck(NULL),
  fProfileSixPartCorrel(NULL),
// Parameters related to the 2-particle correlations with eta gaps.
  fTPCEtaList(NULL)
{
/* Constructor of the class. */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Create the mother list with the rights on everything contained inside it.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE);

// Define the input and output slots.
  DefineOutput(1, TList::Class());

// Initialise all elements in the arrays.
  InitialiseArraysOfQvectors();
  InitialiseArraysOfHistos();
  InitialiseArraysOfTProfiles();

}   // End of the constructor.

/* ------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
/* Destructor of the class. */
  if (fMainList) {delete fMainList;}
}   // End of the destructor.

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Define the outputs of the task at the beginning of the analysis. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()";

// Avoid name clashes.
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

// JEfficiency for NUA correction : DongJo
  fEfficiency = new AliJEfficiency();
  fEfficiency->SetMode(1); // 1: priod should work for you
  fEfficiency->SetDataPath( "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data" );

// Book the lists, histograms and profiles.
  this->BookAllLists();
  this->BookMultiplicityList();
  this->BookEventSelectionList();
  this->BookTrackSelectionList();
  this->BookMPCList();
  this->BookTPCEtaList();

// Continue to avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1, fMainList);
  fFirstEvent = kTRUE;
}   // End of "void UserCreateOutputObjects()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the analysis for each event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";
  fAODevent = dynamic_cast<AliAODEvent*>(InputEvent());   // Pointer to an AOD event.
  fMCevent = MCEvent();                                   // Pointer to a MC event.

// Select the correct method according to the given type of files.
  if ((Int_t)fAnalyseOnlyAOD + (Int_t)fAnalyseOnlyMC + (Int_t)fAnalyseBothMCandAOD != 1)
  {Fatal(sMethodName.Data(), "ERROR: only one 'fProcess' must be set to kTRUE.");}
  else if (fAnalyseOnlyAOD) {AnalyseAODevent();}
  else if (fAnalyseOnlyMC) {AnalyseMCevent();}
  else if (fAnalyseBothMCandAOD) {Fatal(sMethodName.Data(),"TBA: analysis with both AOD and MC files.");}

// PostData.
  PostData(1, fMainList);
}   // End of "void UserExec()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
/* Save the outputs at the end of the execution of the script. */
// Access the mother list.
  fMainList = (TList*)GetOutputData(1);
  if (!fMainList) {exit(1);}

// Create the output file and save the mother list inside.
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fMainList->Write(fMainList->GetName(),TObject::kSingleKey);
  delete outputFile;
}   // End of "void Terminate()".

/* ========================================================================== /
/ Methods called in the constructors.                                         /
/ ========================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfQvectors()
{
/* Initialise to zero all the elements of "fQvectors". */
  for (Int_t iHarmo = 0; iHarmo < 65; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }
}   // End of "void InitialiseArraysOfQvectors()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfHistos()
{
/* Initialise to NULL all the histograms inside arrays. */
  for (Int_t i = 0; i < 2; i++)
  {
    fHistoCentrality[i] = NULL;
    fHistoMultiplicity[i] = NULL;
    fHistoMultiplicityMain[i] = NULL;
    fHistoMultiplicityGlobal[i] = NULL;
    fHistoCorrelatedFilters[i] = NULL;
    fHistoPVx[i] = NULL;
    fHistoPVy[i] = NULL;
    fHistoPVz[i] = NULL;
    fHistoPt[i] = NULL;
    fHistoEta[i] = NULL;
    fHistoPhi[i] = NULL;
    fHistoNTPC[i] = NULL;
    fHistoChiSquare[i] = NULL;
    fHistoNITS[i] = NULL;
    fHistoDCAxy[i] = NULL;
    fHistoDCAz[i] = NULL;
    fHistoCharge[i] = NULL;
  }

  for (Int_t j = 0; j < 8; j++) {fHistoReducedQvectors[j] = NULL;}

}   // End of "void InitialiseArraysOfHistos()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfTProfiles()
{
/* Initialise to NULL all the profiles inside arrays. */
  for (Int_t i = 0; i < 11; i++) {fProfileTPCEta[i] = NULL;}
}   // End of "void InitialiseArraysOfTProfiles()".

/* ========================================================================== /
/ Methods called in "UserExec".                                               /
/ ========================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent()
{
/* Execute the analysis for a given AOD event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseAODevent()";

// Check the presence of an AOD event.
  if (!fAODevent) {Fatal(sMethodName.Data(), "ERROR: no AOD event found.");}

// Set the detector for the centrality estimation as set in the task.
  TString centralityEstimator = "centralityEstimator";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Select only the events belonging to the current centrality range.
  fMultSelection = (AliMultSelection*)fAODevent->FindListObject("MultSelection");
  if (!fMultSelection) {return;}

  Float_t centrality = fMultSelection->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));
  if ( (centrality < fCentralityMin) || (centrality >= fCentralityMax) ) {return;}

  fHistoCentrality[0]->Fill(centrality);

// Look at the initial observables.
  fInitialMultiplicity = fAODevent->GetNumberOfTracks();
  AliAODVertex *primaryVertex = (AliAODVertex*)fAODevent->GetPrimaryVertex();   // 3D position of the PV.

  fHistoMultiplicity[0]->Fill(fInitialMultiplicity);
  fHistoPVx[0]->Fill(primaryVertex->GetX());
  fHistoPVy[0]->Fill(primaryVertex->GetY());
  fHistoPVz[0]->Fill(primaryVertex->GetZ());

// Apply the event selection.
  if(!ApplyPhysicsEventSelection()) {return;}
  if (!RemoveHMOs()) {return;}

  fHistoCentrality[1]->Fill(centrality);
  fHistoPVx[1]->Fill(primaryVertex->GetX());
  fHistoPVy[1]->Fill(primaryVertex->GetY());
  fHistoPVz[1]->Fill(primaryVertex->GetZ());

// Get the run number and load the correct NUE correction (from DongJo).
  if (fFirstEvent)
  {
    fEfficiency->SetRunNumber(fAODevent->GetRunNumber());
    fEfficiency->Load();
    fFirstEvent = kFALSE;
  }

// Apply the track selection to classify the tracks.
  long long multiplicity = fAODevent->GetNumberOfTracks();    // Current number of tracks.
  long long finalMultiplicity = 0;                            // Number of tracks selected for the analysis.
  Bool_t *isTrackSelected = new Bool_t[multiplicity]();       // kTRUE: the track passed the selection, kFALSE: the track is rejected.

  for (Int_t iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // Get the current track belonging to the main filter.
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(iTrack));
    if (!aTrack) {continue;}    // No track found.
    if (!aTrack->TestFilterBit(fMainFilter)) {continue;}    // The track does not belong to the filter.

  // Apply the track selection (histograms are filled inside).
    if (!ApplyTrackSelectionAOD(aTrack)) {isTrackSelected[iTrack] = kFALSE;}
    else {isTrackSelected[iTrack] = kTRUE; finalMultiplicity++;}
  }   // End of the first loop over the tracks.

// Select only the events with enough tracks for the event weight.
  if (finalMultiplicity < fMultiplicityMin) {return;}
  fHistoMultiplicity[1]->Fill(finalMultiplicity);

// Prepare the variables needed for the rest of the analysis.
  Float_t   iPt = 0.;                                               // Transverse momentum.
  Float_t   *iEta = new Float_t[finalMultiplicity]();               // Pseudorapidity for the eta gaps.
  Float_t   *iPhi = new Float_t[finalMultiplicity]();               // Azimuthal angles.
  Float_t   *iParticleWeights = new Float_t[finalMultiplicity]();   // Particle weights.
  Int_t     iIndex = 0;                                             // Index of the selected track in the final arrays.
  Float_t   iEffCorr = 0.;                                          // Efficiency correction.
  Float_t   iEffInverse = 0.;                                       // Inverse of the efficiency correction.

  for (Int_t iiTrack = 0; iiTrack < multiplicity; iiTrack++)
  {
  // Check if iiTrack corresponds to a selected track or not.
    if (!isTrackSelected[iiTrack]) {continue;}

    AliAODTrack *aaTrack = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(iiTrack));
    if (!aaTrack) {continue;}
    if (!aaTrack->TestFilterBit(fMainFilter)) {continue;}

    iPt           = aaTrack->Pt();
    iEta[iIndex]  = aaTrack->Eta();
    iPhi[iIndex]  = aaTrack->Phi();
    iEffCorr      = fEfficiency->GetCorrection(iPt, fFilterbitIndex, centrality);
    iEffInverse   = 1.0/iEffCorr;
    iParticleWeights[iIndex] = 1.;    // Unit particle weights are used by default.
    if (fUseNonUnitParticleWeights) {iParticleWeights[iIndex] = iEffInverse;}

    fHistoEfficiency->Fill(iEffCorr);
    fHistoEffInverse->Fill(iEffInverse);

    iIndex++;
  }   // End of the second loop over the tracks.

// Compute the multiparticle correlations for the current event.
  CalculateQvectors(finalMultiplicity, iPhi, iParticleWeights);
  ComputeReducedQvectors(finalMultiplicity);
  ComputeAllCorrelators(finalMultiplicity, iPhi, iParticleWeights);
  if (fComputeEtaGaps) {ComputeTPCWithEtaGaps(finalMultiplicity, iPhi, iParticleWeights, iEta);}

// Reset the EbE variables.
  centrality              = 0.;
  fInitialMultiplicity    = 0;
  multiplicity            = 0;
  finalMultiplicity       = 0;
  delete [] isTrackSelected;
  iPt                     = 0.;
  delete [] iEta;
  delete [] iPhi;
  delete [] iParticleWeights;
  iIndex                  = 0;
  iEffCorr                = 0.;
  iEffInverse             = 0.;
}   // End of "void AnalyseAODevent()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent()
{
/* Execute the analysis for a given MC event. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseMCevent()";

// Check the presence of an AOD event.
  if (!fMCevent) {Fatal(sMethodName.Data(), "ERROR: no MC event found.");}

// Set the detector for the centrality estimation as set in the task.
  TString centralityEstimator = "centralityEstimator";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  {
    Fatal(sMethodName.Data(), "ERROR: only one detector must be set to kTRUE.");
  }
  else if (fCentralityFromVZero) {centralityEstimator = "V0M";}
  else if (fCentralityFromSPD) {centralityEstimator = "CL1";}

// Select only the events belonging to the current centrality range.
  Float_t centrality = 0.;
  if (!fAnalyseOnlyMC)    // Centrality has no sense for Kine events.
  {
    fMultSelection = (AliMultSelection*)fMCevent->FindListObject("MultSelection");
    if (!fMultSelection) {return;}

    Double_t centrality = fMultSelection->GetMultiplicityPercentile(Form("%s", centralityEstimator.Data()));
    if ( (centrality < fCentralityMin) || (centrality >= fCentralityMax) ) {return;}
  }
  fHistoCentrality[0]->Fill(centrality);

// Look at the initial observables.
  fInitialMultiplicity = fMCevent->GetNumberOfTracks();
  AliMCVertex *primaryVertex = (AliMCVertex*)fMCevent->GetPrimaryVertex();

  fHistoMultiplicity[0]->Fill(fInitialMultiplicity);
  fHistoPVx[0]->Fill(primaryVertex->GetX());
  fHistoPVy[0]->Fill(primaryVertex->GetY());
  fHistoPVz[0]->Fill(primaryVertex->GetZ());

// Apply the physics event selection. (No need to remove HMOs in MC events.)
  if(!ApplyPhysicsEventSelection()) {return;}

  fHistoCentrality[1]->Fill(centrality);
  fHistoPVx[1]->Fill(primaryVertex->GetX());
  fHistoPVy[1]->Fill(primaryVertex->GetY());
  fHistoPVz[1]->Fill(primaryVertex->GetZ());

// Apply the track selection to classify the tracks.
  long long multiplicity = fMCevent->GetNumberOfTracks();   // Current number of tracks.
  long long finalMultiplicity = 0;                          // Number of tracks selected for the analysis.
  Bool_t *isTrackSelected = new Bool_t[multiplicity]();     // kTRUE: the track passed the selection, kFALSE: the track is rejected.

  for (Int_t iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // Get the current track belonging to the main filter.
    AliAODMCParticle *aTrack = dynamic_cast<AliAODMCParticle*>(fMCevent->GetTrack(iTrack));
    if (!aTrack) {continue;}    // No track found.

  // Apply the track selection (histograms are filled inside).
    if (!ApplyTrackSelectionMC(aTrack)) {isTrackSelected[iTrack] = kFALSE;}
    else {isTrackSelected[iTrack] = kTRUE; finalMultiplicity++;}
  }   // End of the first loop over the tracks.

// Select only the events with enough tracks for the event weight.
  if (finalMultiplicity < fMultiplicityMin) {return;}
  fHistoMultiplicity[1]->Fill(finalMultiplicity);

// Prepare the variables needed for the rest of the analysis.
  Float_t   iPt = 0.;                                               // Transverse momentum.
  Float_t   *iEta = new Float_t[finalMultiplicity]();               // Pseudorapidity for the eta gaps.
  Float_t   *iPhi = new Float_t[finalMultiplicity]();               // Azimuthal angles.
  Float_t   *iParticleWeights = new Float_t[finalMultiplicity]();   // Particle weights.
  Int_t     iIndex = 0;                                             // Index of the selected track in the final arrays.

  for (Int_t iiTrack = 0; iiTrack < multiplicity; iiTrack++)
  {
  // Check if iiTrack corresponds to a selected track or not.
    if (!isTrackSelected[iiTrack]) {continue;}

    AliAODMCParticle *aaTrack = dynamic_cast<AliAODMCParticle*>(fMCevent->GetTrack(iiTrack));
    if (!aaTrack) {continue;}

    iPt           = aaTrack->Pt();
    iEta[iIndex]  = aaTrack->Eta();
    iPhi[iIndex]  = aaTrack->Phi();
    iParticleWeights[iIndex] = 1.;    // Unit particle weights are used by default.
    if (fUseNonUnitParticleWeights) {Fatal(sMethodName.Data(), "ERROR: TBA use of non-unit particle weights.");}

    iIndex++;
  }   // End of the second loop over the tracks.

// Compute the multiparticle correlations for the current event.
  CalculateQvectors(finalMultiplicity, iPhi, iParticleWeights);
  ComputeReducedQvectors(finalMultiplicity);
  ComputeAllCorrelators(finalMultiplicity, iPhi, iParticleWeights);
  if (fComputeEtaGaps) {ComputeTPCWithEtaGaps(finalMultiplicity, iPhi, iParticleWeights, iEta);}

// Do a dummy expression for iPt.
  iPt++;

// Reset the EbE variables.
  centrality              = 0.;
  fInitialMultiplicity    = 0;
  multiplicity            = 0;
  finalMultiplicity       = 0;
  delete [] isTrackSelected;
  iPt                     = 0.;
  delete [] iEta;
  delete [] iPhi;
  delete [] iParticleWeights;
  iIndex                  = 0;
}   // End of "void AnalyseMCevent()".

/* ------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyPhysicsEventSelection()
{
/* Apply the part of the event selection related to the physics. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::ApplyPhysicsEventSelection()";

// Select the right method according to the types of files.
  Float_t PVx = 0.;
  Float_t PVy = 0.;
  Float_t PVz = 0.;

  if (fAnalyseOnlyAOD)
  {
    AliAODVertex *primaryVertex = (AliAODVertex*)fAODevent->GetPrimaryVertex();
    PVx = primaryVertex->GetX();
    PVy = primaryVertex->GetY();
    PVz = primaryVertex->GetZ();
  }
  else if (fAnalyseOnlyMC)
  {
    AliMCVertex *primaryVertex = (AliMCVertex*)fMCevent->GetPrimaryVertex();
    PVx = primaryVertex->GetX();
    PVy = primaryVertex->GetY();
    PVz = primaryVertex->GetZ();
  }

// Apply the selection.
  if (fCutOnPVx)
  {if ( (PVx < fPVxMin) || (PVx > fPVxMax) ) {return kFALSE;}}
  if (fCutOnPVy)
  {if ( (PVy < fPVyMin) || (PVy > fPVyMax) ) {return kFALSE;}}
  if (fCutOnPVz)
  {if ( (PVz < fPVzMin) || (PVz > fPVzMax) ) {return kFALSE;}}

// Reset the variables.
  PVx = 0.;
  PVy = 0.;
  PVz = 0.;

  return kTRUE;
}   // End of "Bool_t ApplyPhysicsEventSelection()".

/* ------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::RemoveHMOs()
{
/* Get and apply the criteria to remove the HMOs in AOD events. */
  long long multiplicity = fAODevent->GetNumberOfTracks();    // Number of tracks after the physics event selection.
  long long globalMultiplicity = 0;                           // Number of tracks in the global filter.
  long long mainMultiplicity = 0;                             // Number of tracks in the main filter.

  for (Int_t iTrack = 0; iTrack < multiplicity; iTrack++)
  {
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(iTrack));
    if (!aTrack) {continue;}
    if (aTrack->TestFilterBit(fGlobalFilter)) {globalMultiplicity++;}
    if (aTrack->TestFilterBit(fMainFilter)) {mainMultiplicity++;}
  }

  fHistoMultiplicityGlobal[0]->Fill(globalMultiplicity);
  fHistoMultiplicityMain[0]->Fill(mainMultiplicity);
  if (fDoCorrelationsHisto) {fHistoCorrelatedFilters[0]->Fill(globalMultiplicity, mainMultiplicity);}

// Apply the HMOs cuts.
  if (fCutOnHMOs)
  {
    if ( (Float_t)mainMultiplicity < (fMultiplicityMinA*(Float_t)globalMultiplicity + fMultiplicityMinB) ) {return kFALSE;} // The number of tracks in the main filter is under the accepted band.
    if ( (Float_t)mainMultiplicity > (fMultiplicityMaxA*(Float_t)globalMultiplicity + fMultiplicityMaxB) ) {return kFALSE;} // The number of tracks in the main filter is above the accepted band.
  }

  fHistoMultiplicityGlobal[1]->Fill(globalMultiplicity);
  fHistoMultiplicityMain[1]->Fill(mainMultiplicity);
  if (fDoCorrelationsHisto) {fHistoCorrelatedFilters[1]->Fill(globalMultiplicity, mainMultiplicity);}

// Reset the variables.
  multiplicity        = 0;
  globalMultiplicity  = 0;
  mainMultiplicity    = 0;

  return kTRUE;
}   // End of "Bool_t RemoveHMOs()".

/* ------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionAOD(AliAODTrack *aAODtrack)
{
/* Apply the track selection criteria and return if the track passed it or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)";

// Get all the observables for the given AOD track.
  Float_t   pT          = aAODtrack->Pt();                            // Transverse momentum.
  Float_t   eta         = aAODtrack->Eta();                           // Pseudorapidity.
  Float_t   phi         = aAODtrack->Phi();                           // Azimuthal angle.
  Int_t     NTPC        = aAODtrack->GetTPCNcls();                    // Number of TPC clusters.
  Float_t   chiSquare   = (aAODtrack->GetTPCchi2())/(aAODtrack->GetNcls(1));    // Chi square of p in the TPC.
  Int_t     NITS        = aAODtrack->GetITSNcls();                    // Number of ITS clusters.
  Float_t   DCAx        = aAODtrack->XAtDCA();                        // DCA along x.
  Float_t   DCAy        = aAODtrack->YAtDCA();                        // DCA along y.
  Float_t   DCAz        = aAODtrack->ZAtDCA();                        // DCA along z.
  Int_t     charge      = aAODtrack->Charge();                        // Electric charge.
  Float_t   DCAxy       = TMath::Sqrt( (DCAx*DCAx) + (DCAy*DCAy) );   // DCA in the xy-plane.

// Draw the initial distributions.
  fHistoPt[0]->Fill(pT);
  fHistoEta[0]->Fill(eta);
  fHistoPhi[0]->Fill(phi);
  fHistoNTPC[0]->Fill(NTPC);
  fHistoChiSquare[0]->Fill(chiSquare);
  fHistoNITS[0]->Fill(NITS);
  fHistoDCAxy[0]->Fill(DCAxy);
  fHistoDCAz[0]->Fill(DCAz);
  fHistoCharge[0]->Fill(charge);

// Apply the selection.
  if (fCutOnPt)
  {if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;}}
  if (fCutOnEta)
  {if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;}}
  if (fCutOnNTPC)
  {if (NTPC < fNTPCMin) {return kFALSE;}}
  if (fCutOnChi)
  {if ( (chiSquare < fChiMin) || (chiSquare > fChiMax) ) {return kFALSE;}}
  if (fCutOnNITS)
  {
    if (fMainFilter == 128) {Fatal(sMethodName.Data(), "ERROR: TPConly tracks do not have ITS information.");}
    else if (NITS < fNITSMin) {return kFALSE;}
  }
  if (fCutOnDCAxy)
  {if (DCAxy > fDCAxyMax) {return kFALSE;}}
  if (fCutOnDCAz)
  {if (DCAz > fDCAzMax) {return kFALSE;}}
  if (fCutOnCharge)
  {
    if (fKeepPositiveCharges) {if (charge < 0) {return kFALSE;}}          // Keep only the positive tracks.
    else if (!fKeepPositiveCharges) {if (charge > 0) {return kFALSE;}}    // Keep only the negative tracks.
    else {Fatal(sMethodName.Data(), "ERROR: Select the sign of the charges to keep.");}
  }

// Draw the final distributions.
  fHistoPt[1]->Fill(pT);
  fHistoEta[1]->Fill(eta);
  fHistoPhi[1]->Fill(phi);
  fHistoNTPC[1]->Fill(NTPC);
  fHistoChiSquare[1]->Fill(chiSquare);
  fHistoNITS[1]->Fill(NITS);
  fHistoDCAxy[1]->Fill(DCAxy);
  fHistoDCAz[1]->Fill(DCAz);
  fHistoCharge[1]->Fill(charge);

// Reset the variables.
  pT        = 0.;
  eta       = 0.;
  phi       = 0.;
  NTPC      = 0;
  chiSquare = 0.;
  NITS      = 0;
  DCAx      = 0.;
  DCAy      = 0.;
  DCAz      = 0.;
  charge    = 0;
  DCAxy     = 0.;

  return kTRUE;
}   // End of "Bool_t ApplyTrackSelectionAOD()".

/* ------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionMC(AliAODMCParticle *aMCtrack)
{
/* Apply the track selection criteria to the given MC track and return if it passed or not. */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelectionMC(AliAODMCParticle *aMCtrack)";

// Get all the observables for the given MC track.
  Float_t   pT          = aMCtrack->Pt();       // Transverse momentum.
  Float_t   eta         = aMCtrack->Eta();      // Pseudorapidity.
  Float_t   phi         = aMCtrack->Phi();      // Azimuthal angle.
  Int_t     charge      = aMCtrack->Charge();   // Electric charge.

// Draw the initial distributions.
  fHistoPt[0]->Fill(pT);
  fHistoEta[0]->Fill(eta);
  fHistoPhi[0]->Fill(phi);
  fHistoCharge[0]->Fill(charge);

// Apply the track selection.
  if (fCutOnPt)
  {if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;}}
  if (fCutOnEta)
  {if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;}}
  if (fCutOnCharge)
  {
    if (fKeepPositiveCharges) {if (charge < 0) {return kFALSE;}}          // Keep only the positive tracks.
    else if (!fKeepPositiveCharges) {if (charge > 0) {return kFALSE;}}    // Keep only the negative tracks.
    else {Fatal(sMethodName.Data(), "ERROR: Select the sign of the charges to keep.");}
  }

// Draw the final distributions.
  fHistoPt[1]->Fill(pT);
  fHistoEta[1]->Fill(eta);
  fHistoPhi[1]->Fill(phi);
  fHistoCharge[1]->Fill(charge);

// Reset the variables.
  pT        = 0.;
  eta       = 0.;
  phi       = 0.;
  charge    = 0;

  return kTRUE;
}   // End of "Bool_t ApplyTrackSelectionMC()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Calculate all the Q-vectors for the given arrays of azimuthal angles and particle weights. */
  Float_t   iAngle            = 0.;   // Azimuthal angle of the current track.
  Float_t   iWeight           = 0.;   // Particle weight of the current track.
  Float_t   iWeightToPowerP   = 0.;   // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 65; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (long long iTrack = 0; iTrack < numberOfParticles; iTrack++)
  {
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 65; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < 9; iPower++)
      {
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle), iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }   // End of the loop over the powers.
    }   // End of the loop over the harmonics.
  }   // End of the loop over the tracks.

// Reset the variables.
  iAngle            = 0.;
  iWeight           = 0.;
  iWeightToPowerP   = 0.;
}   // End of "void CalculateQvectors()".

/* ------------------------------------------------------------------------- */
TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
/* Simplify the use of the Q-vectors in the next methods. */
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]);   // Use that Q*(n,p) = Q(-n,p).
}   // End of "TComplex Q()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::ComputeReducedQvectors(long long numberOfParticles)
{
/* Compute the modulus of reduced Q-vector for the current event and save it in the correct histogram. */
  Double_t  reducedQ    = 0.;   // Modulus of reduced Q-vector.
  Int_t     iHarmo      = 0;    // Harmonic corresponding to the current histogram.

// Loop over the harmonics to get each q_n.
  for (Int_t i = 0; i < 8; i++)
  {
    iHarmo = i + 1;
    reducedQ = Q(iHarmo, fReducedQPower).Rho()/(TMath::Sqrt(1.*numberOfParticles));
    fHistoReducedQvectors[i]->Fill(reducedQ);
  }

// Reset the variables.
  reducedQ  = 0.;
  iHarmo    = 0;
}   // End of "void ComputeReducedQvectors()".

/* ------------------------------------------------------------------------- */
TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
/* Calculate the multi-particle correlators by using the recursion method.                  /
/ Improved, faster version originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). */
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
}   // End of "TComplex CalculateRecursion()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::ComputeAllCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Compute and save all the multiparticle correlators for the current event. */
  Int_t       twoZerosArray[2]    = {0};                // 2d array with zeros for the denominator.
  Int_t       twoPartArray[2]     = {0};                // 2d array for the current harmonic.
  Int_t       fourZerosArray[4]   = {0};                // 4d array with zeros for the denominator.
  Int_t       fourPartArray[4]    = {0};                // 4d array for the current couple of harmonics.
  Int_t       sixZerosArray[6]    = {0};                // 6d array with zeros for the denominator.
  Int_t       sixPartArray[6]     = {0};                // 6d array for the current triplet of harmonics.
  Double_t    eventWeight         = 0.;                 // Event weight and denominator of the correlator.
  TComplex    complexCorrelator   = TComplex(0., 0.);   // Complex form of the correlator.
  Double_t    realCorrelator      = 0.;                 // Real part of the correlator.
  Int_t       iBin                = 1;                  // Bin for the current correlator.
  Double_t    iMiddleBin          = 0.;                 // Middle of the bin for the current correlator.

// Compute the 2-particle correlators for v_1 to v_8.
  eventWeight = (CalculateRecursion(2, twoZerosArray)).Re();
  for (Int_t n = 1; n < 9; n++)   // Loop over the harmonics.
  {
  // Compute the correlator.
    twoPartArray[0]     = n;
    twoPartArray[1]     = -1*n;
    complexCorrelator   = (CalculateRecursion(2, twoPartArray))/eventWeight;
    realCorrelator      = complexCorrelator.Re();

  // Fill the corresponding bin in the profile.
    iMiddleBin          = (1.*n) - 0.5;
    fProfileTwoPartCorrel->Fill(iMiddleBin, realCorrelator, eventWeight);
    fProfileTwoPartCorrel->GetXaxis()->SetBinLabel(n, Form("%d", n));

  // Reset the variables for the next harmonic.
    complexCorrelator   = TComplex(0., 0.);
    realCorrelator      = 0.;
    iMiddleBin          = 0.;
  }   // End of the loop over the harmonics for the 2-p correlators.

  eventWeight           = 0.;   // Reset of the denominator.

// Compute the 4-particle correlators.
  eventWeight           = (CalculateRecursion(4, fourZerosArray)).Re();
  for(Int_t m = 1; m < 9; m++)
  {
    for (Int_t n = m; n < 9; n++)
    {
    // Compute the correlator.
      fourPartArray[0]  = m;
      fourPartArray[1]  = n;
      fourPartArray[2]  = -1*m;
      fourPartArray[3]  = -1*n;
      complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
      realCorrelator    = complexCorrelator.Re();

    // Fill the corresponding bin in the profile.
      iMiddleBin        = (1.*iBin) - 0.5;
      fProfileFourPartCorrel->Fill(iMiddleBin, realCorrelator, eventWeight);
      fProfileFourPartCorrel->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

    // Reset of the variables for the next pair of harmonics.
      iBin++;
      complexCorrelator = TComplex(0., 0.);
      realCorrelator    = 0.;
      iMiddleBin        = 0.;
    }   // End of the loop over the second harmonic for the 4-p correlators.
  }   // End of the loop over the first harmonic for the 4-p correlators.

  iBin                  = 1;

  for(Int_t m = 2; m < 9; m++)
  {
    for (Int_t n = 1; n < m; n++)
    {
    // Compute the correlator.
      fourPartArray[0]  = m;
      fourPartArray[1]  = n;
      fourPartArray[2]  = -1*m;
      fourPartArray[3]  = -1*n;
      complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
      realCorrelator    = complexCorrelator.Re();

    // Fill the corresponding bin in the profile.
      iMiddleBin        = (1.*iBin) - 0.5;
      fProfileFourPartCorrelCheck->Fill(iMiddleBin, realCorrelator, eventWeight);
      fProfileFourPartCorrelCheck->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

    // Reset of the variables for the next pair of harmonics.
     iBin++;
     complexCorrelator  = TComplex(0., 0.);
     realCorrelator     = 0.;
     iMiddleBin         = 0.;
    }   // End of the loop over the second harmonic for cross-check.
  }   // End of the loop over the first harmonic for cross-check.

  eventWeight           = 0.;
  iBin                  = 1;

// Compute the needed 6-particle correlators.
  eventWeight           = (CalculateRecursion(6, sixZerosArray)).Re();
  for (Int_t l = 1; l < 7; l++)
  {
    for (Int_t m = 2; m < 8; m++)
    {
      if (l >= m) {continue;}   // Remove the cases in (l,l,n).
      for (Int_t n = 3; n < 9; n++)
      {
        if ((l >= n) || (m >= n)) {continue;}   // Remove the cases in (l,m,l) or (l,m,m).

      // Compute the correlator.
        sixPartArray[0] = l;
        sixPartArray[1] = m;
        sixPartArray[2] = n;
        sixPartArray[3] = -1*l;
        sixPartArray[4] = -1*m;
        sixPartArray[5] = -1*n;
        complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
        realCorrelator  = complexCorrelator.Re();

      // Fill the corresponding bin in the TProfile.
        iMiddleBin = (1.*iBin) - 0.5;
        fProfileSixPartCorrel->Fill(iMiddleBin, realCorrelator, eventWeight);
        fProfileSixPartCorrel->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d,%d)", l, m, n));

      // Reset of the variables for the next pair of harmonics.
        iBin++;
        complexCorrelator = TComplex(0., 0.);
        realCorrelator  = 0.;
        iMiddleBin      = 0.; 
      }   // End of the loop over the third harmonic of the triplet.
    }   // End of the loop over the second harmonic of the triplet.
  }   // End of the loop over for the first harmonic of the triplet.

// Reset the variables for the next event.
  eventWeight           = 0.;
  complexCorrelator     = TComplex(0., 0.);
  realCorrelator        = 0.;
  iBin                  = 1;
  iMiddleBin            = 0.;

}   // End of "void ComputeAllCorrelators()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::ComputeTPCWithEtaGaps(long long numberOfParticles, Float_t angles[], Float_t pWeights[], Float_t pseudorapidity[])
{
/* Compute the 2-particle correlators using eta gaps for the current event. */
  TComplex  Qminus[11][8]   = {{TComplex(0., 0.)}};   // Q-vectors for the negative subset of the eta range, for v_1 to v_8.
  TComplex  Qplus[11][8]    = {{TComplex(0., 0.)}};   // Q-vectors for the positive subset of the eta range, for v_1 to v_8.
  Float_t   Mminus[11][8]   = {{0.}};                 // Multiplicity in the negative subset of the eta range.
  Float_t   Mplus[11][8]    = {{0.}};                 // Multiplicity in the positive subset of the eta range.
  Float_t   etaGaps[11]     = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};  // Possible values for the gaps.
  Float_t   iAngle          = 0.;                     // Azimuthal angle of the current particle.
  Float_t   iWeight         = 1.;                     // Particle weight of the current particle (default: unit weight).
  Float_t   iEta            = 0.;                     // Pseudorapidity of the current particle.
  Float_t   iWeightToP      = 1.;                     // Particle weight rised to the power p.
  TComplex  complexCorrel   = TComplex(0., 0.);       // Complex value of the 2-p correlator.
  Double_t  realCorrel      = 0.;                     // Real value of the 2-p correlator.

// Compute the Q-vectors for the negative and positive subsets of the eta range.
  for (Int_t iPart = 0; iPart < numberOfParticles; iPart++)
  {
  // Read the right elements in the provided arrays.
    iAngle  = angles[iPart];
    iWeight = pWeights[iPart];
    iEta    = pseudorapidity[iPart];
    if (fUseNonUnitParticleWeights) {iWeightToP = iWeight;}   // All weights are multiplied to get the final one.

  // Compute the Q-vectors.
    if (iEta < 0.)    // Negative subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta < ((-0.5)*etaGaps[iGap]))    // Compute only if the particle is in the range.
          {
            Qminus[iGap][iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mminus[iGap][iHarmo] += iWeightToP;
          }
          else {continue;}
        }   // End of the loop over the gaps.
      }   // End of the loop over the harmonics.
    }   // End of the condition "negative subset".
    else if (iEta > 0.)   // Positive subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta > (0.5*etaGaps[iGap]))   // Compute only if the particle is in the range.
          {
            Qplus[iGap][iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mplus[iGap][iHarmo] += iWeightToP;
          }
        }   // End of the loop over the gaps.
      }   // End of the loop over the harmonics.
    }   // End of the condition "positive subset".
    else {continue;}    // Particle with iEta = 0.
  }   // End of the loop over the particles for the Q-vectors.

// Compute the 2-p correlators using Qminus and Qplus.
  for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
  {
    for (Int_t iGap = 0; iGap < 11; iGap++)
    {
      if (!( (Qminus[iGap][iHarmo].TComplex::Rho() > 0.) && (Qplus[iGap][iHarmo].TComplex::Rho() > 0.) )) {continue;}
      if (!( (Mminus[iGap][iHarmo] > 0.) && (Mplus[iGap][iHarmo] > 0.) )) {continue;}

      complexCorrel = Qminus[iGap][iHarmo]*TComplex::Conjugate(Qplus[iGap][iHarmo]);
      realCorrel    = (complexCorrel.Re())/(Mminus[iGap][iHarmo]*Mplus[iGap][iHarmo]);
      fProfileTPCEta[iGap]->Fill(iHarmo + 0.5, realCorrel, Mminus[iGap][iHarmo]*Mplus[iGap][iHarmo]);

    // Reset the 2-particle correlator.
      complexCorrel = TComplex(0.,0.);
      realCorrel    = 0.;
    }   // End of the loop over the gaps.
  }   // End of the loop over the harmonics.

}   // End of "void ComputeTPCWithEtaGaps()".

/* ========================================================================== /
/ Methods called in "Terminate".                                              /
/ ========================================================================== */

/* ========================================================================== /
/ Methods called in "UserCreateOutputObjects".                                /
/ ========================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists in fMainList. */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()";

// Check if the mother list exists.
  if (!fMainList) {Fatal(sMethodName.Data(), "ERROR: 'fMainList' does not exist.");}

// Daughter list with the multiplicity histograms.
  fMultiplicityList = new TList();
  fMultiplicityList->SetName("fMultiplicityList");
  fMultiplicityList->SetOwner(kTRUE);
  fMainList->Add(fMultiplicityList);

// Daughter list for the event histograms.
  fEventSelectionList = new TList();
  fEventSelectionList->SetName("fEventSelectionList");
  fEventSelectionList->SetOwner(kTRUE);
  fMainList->Add(fEventSelectionList);

// Daughter list for the track histograms.
  fTrackSelectionList = new TList();
  fTrackSelectionList->SetName("fTrackSelectionList");
  fTrackSelectionList->SetOwner(kTRUE);
  fMainList->Add(fTrackSelectionList);

// Daughter list for the MPC histograms.
  fMPCList = new TList();
  fMPCList->SetName("fMPCList");
  fMPCList->SetOwner(kTRUE);
  fMainList->Add(fMPCList);

// Daughter list for the eta gaps profiles.
  fTPCEtaList = new TList();
  fTPCEtaList->SetName("fTPCEtaList");
  fTPCEtaList->SetOwner(kTRUE);
  fMainList->Add(fTPCEtaList);
}   // End of "void BookAllLists()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMultiplicityList()
{
/* Book the histograms related to the multiplicity distributions. */
  for (Int_t i = 0; i < 2; i++)
  {
    TString iState = "Final";
    if (i == 0) {iState = "Initial";}

  // Centrality distributions.
    fHistoCentrality[i] = new TH1F("", "", 100, 0., 100.);
    fHistoCentrality[i]->SetName(Form("fHistoCentrality%s", iState.Data()));
    fHistoCentrality[i]->SetTitle(Form("%s distribution of the centrality", iState.Data()));
    fHistoCentrality[i]->SetStats(kTRUE);
    fHistoCentrality[i]->GetXaxis()->SetTitle("Centrality percentile");
    fHistoCentrality[i]->GetYaxis()->SetTitle("Number of events");
    fMultiplicityList->Add(fHistoCentrality[i]);

  // Multiplicity distributions.
    fHistoMultiplicity[i] = new TH1I("", "", 30000, 0., 30000.);
    fHistoMultiplicity[i]->SetName(Form("fHistoMultiplicity%s", iState.Data()));
    fHistoMultiplicity[i]->SetTitle(Form("%s distribution of the number of tracks", iState.Data()));
    fHistoMultiplicity[i]->SetStats(kTRUE);
    fHistoMultiplicity[i]->GetXaxis()->SetTitle("Number of tracks");
    fHistoMultiplicity[i]->GetYaxis()->SetTitle("Number of events");
    fMultiplicityList->Add(fHistoMultiplicity[i]);

  // Multiplicity distributions meaningful only for AOD files.
    if (fAnalyseOnlyAOD)
    {
    // Main filter bit.
      fHistoMultiplicityMain[i] = new TH1I("", "", 30000, 0., 30000.);
      fHistoMultiplicityMain[i]->SetName(Form("fHistoMultiplicityMain%s", iState.Data()));
      fHistoMultiplicityMain[i]->SetTitle(Form("%s distribution of the number of tracks (filter %d)", iState.Data(), fMainFilter));
      fHistoMultiplicityMain[i]->SetStats(kTRUE);
      fHistoMultiplicityMain[i]->GetXaxis()->SetTitle("Number of tracks");
      fHistoMultiplicityMain[i]->GetYaxis()->SetTitle("Number of events");
      fMultiplicityList->Add(fHistoMultiplicityMain[i]);

    // Multiplicity distributions for the global filter bit.
      fHistoMultiplicityGlobal[i] = new TH1I("", "", 30000, 0., 30000.);
      fHistoMultiplicityGlobal[i]->SetName(Form("fHistoMultiplicityGlobal%s", iState.Data()));
      fHistoMultiplicityGlobal[i]->SetTitle(Form("%s distribution of the number of tracks (filter %d)", iState.Data(), fGlobalFilter));
      fHistoMultiplicityGlobal[i]->SetStats(kTRUE);
      fHistoMultiplicityGlobal[i]->GetXaxis()->SetTitle("Number of tracks");
      fHistoMultiplicityGlobal[i]->GetYaxis()->SetTitle("Number of events");
      fMultiplicityList->Add(fHistoMultiplicityGlobal[i]);
    }

  // Correlations between filters.
    if (fDoCorrelationsHisto)   // Manually enabled only to decide the cuts, as it is memory-greedy.
    {
      fHistoCorrelatedFilters[i] = new TH2I("", "", 5000, 0., 5000., 5000, 0., 5000.);
      fHistoCorrelatedFilters[i]->SetName(Form("fHistoCorrelatedFilters%s", iState.Data()));
      fHistoCorrelatedFilters[i]->SetTitle(Form("%s correlated tracks", iState.Data()));
      fHistoCorrelatedFilters[i]->SetStats(kTRUE);
      fHistoCorrelatedFilters[i]->GetXaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fGlobalFilter));
      fHistoCorrelatedFilters[i]->GetYaxis()->SetTitle(Form("Multiplicity_{Filter %d}", fMainFilter));
      fHistoCorrelatedFilters[i]->SetMarkerSize(0.5);
      fHistoCorrelatedFilters[i]->SetMarkerColor(kBlue);
      fHistoCorrelatedFilters[i]->SetMarkerStyle(kFullCircle);
      fMultiplicityList->Add(fHistoCorrelatedFilters[i]);
    }
  }   // End of the loop over the states (initial, final).
}   // End of "void BookMultiplicityList()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookEventSelectionList()
{
/* Book the control histograms for the event selection criteria. */
  for (Int_t i = 0; i < 2; i++)
  {
    TString iState = "Final";
    if (i == 0) {iState = "Initial";}

  // Primary Vertex x.
    fHistoPVx[i] = new TH1F("", "", 1000, -20., 20.);
    fHistoPVx[i]->SetName(Form("fHistoPVx%s", iState.Data()));
    fHistoPVx[i]->SetTitle(Form("%s distribution of PV_{x}", iState.Data()));
    fHistoPVx[i]->SetStats(kTRUE);
    fHistoPVx[i]->GetXaxis()->SetTitle("PV_{x} [cm]");
    fHistoPVx[i]->GetYaxis()->SetTitle("dN/dPV_{x}");
    fEventSelectionList->Add(fHistoPVx[i]);

  // Primary Vertex y.
    fHistoPVy[i] = new TH1F("", "", 1000, -20., 20.);
    fHistoPVy[i]->SetName(Form("fHistoPVy%s", iState.Data()));
    fHistoPVy[i]->SetTitle(Form("%s distribution of PV_{y}", iState.Data()));
    fHistoPVy[i]->SetStats(kTRUE);
    fHistoPVy[i]->GetXaxis()->SetTitle("PV_{y} [cm]");
    fHistoPVy[i]->GetYaxis()->SetTitle("dN/dPV_{y}");
    fEventSelectionList->Add(fHistoPVy[i]);

  // Primary Vertex z.
    fHistoPVz[i] = new TH1F("", "", 1000, -20., 20.);
    fHistoPVz[i]->SetName(Form("fHistoPVz%s", iState.Data()));
    fHistoPVz[i]->SetTitle(Form("%s distribution of PV_{z}", iState.Data()));
    fHistoPVz[i]->SetStats(kTRUE);
    fHistoPVz[i]->GetXaxis()->SetTitle("PV_{z} [cm]");
    fHistoPVz[i]->GetYaxis()->SetTitle("dN/dPV_{z}");
    fEventSelectionList->Add(fHistoPVz[i]);
  }   // End of the loop over the states (initial, final).
}   // End of "void BookEventSelectionList()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookTrackSelectionList()
{
/* Book the control histograms for the track selection criteria. */
  for (Int_t i = 0; i < 2; i++)
  {
    TString iState = "Final";
    if (i == 0) {iState = "Initial";}

  // Transverse momentum.
    fHistoPt[i] = new TH1F("", "", 1000, 0., 20.);
    fHistoPt[i]->SetName(Form("fHistoPt%s", iState.Data()));
    fHistoPt[i]->SetTitle(Form("%s distribution of p_{T}", iState.Data()));
    fHistoPt[i]->SetStats(kTRUE);
    fHistoPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistoPt[i]->GetYaxis()->SetTitle("dN/dp_{T}");
    fTrackSelectionList->Add(fHistoPt[i]);

  // Pseudorapidity.
    fHistoEta[i] = new TH1F("", "", 1000, -1., 1.);
    fHistoEta[i]->SetName(Form("fHistoEta%s", iState.Data()));
    fHistoEta[i]->SetTitle(Form("%s distribution of #eta", iState.Data()));
    fHistoEta[i]->SetStats(kTRUE);
    fHistoEta[i]->GetXaxis()->SetTitle("#eta");
    fHistoEta[i]->GetYaxis()->SetTitle("dN/d#eta");
    fTrackSelectionList->Add(fHistoEta[i]);

  // Azimuthal angles.
    fHistoPhi[i] = new TH1F("", "", 1000, 0., 6.3);
    fHistoPhi[i]->SetName(Form("fHistoPhi%s", iState.Data()));
    fHistoPhi[i]->SetTitle(Form("%s distribution of #phi", iState.Data()));
    fHistoPhi[i]->SetStats(kTRUE);
    fHistoPhi[i]->GetXaxis()->SetTitle("#phi");
    fHistoPhi[i]->GetYaxis()->SetTitle("dN/d#phi");
    fTrackSelectionList->Add(fHistoPhi[i]);

  // Number of TPC clusters.
    fHistoNTPC[i] = new TH1I("", "", 170, 0., 170.);
    fHistoNTPC[i]->SetName(Form("fHistoNTPC%s", iState.Data()));
    fHistoNTPC[i]->SetTitle(Form("%s distribution of N_{TPC}", iState.Data()));
    fHistoNTPC[i]->SetStats(kTRUE);
    fHistoNTPC[i]->GetXaxis()->SetTitle("N_{TPC}");
    fHistoNTPC[i]->GetYaxis()->SetTitle("dN/dN_{TPC}");
    fTrackSelectionList->Add(fHistoNTPC[i]);

  // chi^2/NDF of p per TPC clusters.
    fHistoChiSquare[i] = new TH1F("", "", 1000, 0., 20.);
    fHistoChiSquare[i]->SetName(Form("fHistoChiSquare%s", iState.Data()));
    fHistoChiSquare[i]->SetTitle(Form("%s distribution of #chi^{2}/NDF", iState.Data()));
    fHistoChiSquare[i]->SetStats(kTRUE);
    fHistoChiSquare[i]->GetXaxis()->SetTitle("#chi^{2}/NDF of p in TPC");
    fHistoChiSquare[i]->GetYaxis()->SetTitle("Number of tracks");
    fTrackSelectionList->Add(fHistoChiSquare[i]);

  // Number of ITS clusters.
    fHistoNITS[i] = new TH1I("", "", 14, 0., 7.);
    fHistoNITS[i]->SetName(Form("fHistoNITS%s", iState.Data()));
    fHistoNITS[i]->SetTitle(Form("%s distribution of N_{ITS}", iState.Data()));
    fHistoNITS[i]->SetStats(kTRUE);
    fHistoNITS[i]->GetXaxis()->SetTitle("Number of ITS clusters");
    fHistoNITS[i]->GetYaxis()->SetTitle("dN/dN_{ITS}");
    fTrackSelectionList->Add(fHistoNITS[i]);

  // DCA in the xy-plane.
    fHistoDCAxy[i] = new TH1F("", "", 1000, 0., 10.);
    fHistoDCAxy[i]->SetName(Form("fHistoDCAxy%s", iState.Data()));
    fHistoDCAxy[i]->SetTitle(Form("%s distribution of DCA_{xy}", iState.Data()));
    fHistoDCAxy[i]->SetStats(kTRUE);
    fHistoDCAxy[i]->GetXaxis()->SetTitle("DCA_{xy} [cm]");
    fHistoDCAxy[i]->GetYaxis()->SetTitle("dN/dDCA_{xy}");
    fTrackSelectionList->Add(fHistoDCAxy[i]);

  // DCA along z.
    fHistoDCAz[i] = new TH1F("", "", 1000, 0., 10.);
    fHistoDCAz[i]->SetName(Form("fHistoDCAz%s", iState.Data()));
    fHistoDCAz[i]->SetTitle(Form("%s distribution of DCA_{z}", iState.Data()));
    fHistoDCAz[i]->SetStats(kTRUE);
    fHistoDCAz[i]->GetXaxis()->SetTitle("DCA_{z} [cm]");
    fHistoDCAz[i]->GetYaxis()->SetTitle("dN/dDCA_{z}");
    fTrackSelectionList->Add(fHistoDCAz[i]);

  // Electric charge of the tracks.
    fHistoCharge[i] = new TH1I("", "", 2, -2, 2);
    fHistoCharge[i]->SetName(Form("fHistoCharge%s", iState.Data()));
    fHistoCharge[i]->SetTitle(Form("%s distribution of the charges", iState.Data()));
    fHistoCharge[i]->SetStats(kTRUE);
    fHistoCharge[i]->GetXaxis()->SetTitle("Charge");
    fHistoCharge[i]->GetYaxis()->SetTitle("Number of tracks");
    fTrackSelectionList->Add(fHistoCharge[i]);
  }   // End of the loop over the states (initial, final).

  fHistoEfficiency = new TH1F("fHistoEfficiency",
    "Distribution of the efficiency correction", 1000, 0., 20.);
  fHistoEfficiency->SetStats(kTRUE);
  fHistoEfficiency->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoEfficiency->GetYaxis()->SetTitle("dN/dp_{T}");
  fTrackSelectionList->Add(fHistoEfficiency);

  fHistoEffInverse = new TH1F("fHistoEffInverse",
    "Distribution of the inverse of the efficiency correction", 1000, 0., 20.);
  fHistoEffInverse->SetStats(kTRUE);
  fHistoEffInverse->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoEffInverse->GetYaxis()->SetTitle("dN/dp_{T}");
  fTrackSelectionList->Add(fHistoEffInverse);

}   // End of "void BookTrackSelectionList()".
/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMPCList()
{
/* Book the histograms and profiles for the MPC. */
// Modulus of reduced Q-vectors.
  for (Int_t i = 0; i < 8; i++)
  {
    fHistoReducedQvectors[i] = new TH1F("", "", 1000, 0., 10.);
    fHistoReducedQvectors[i]->SetName(Form("fHistoReducedQvectors%d", i));
    fHistoReducedQvectors[i]->SetTitle(Form("Distribution of q_{%d}", i+1));
    fHistoReducedQvectors[i]->SetStats(kTRUE);
    fHistoReducedQvectors[i]->GetXaxis()->SetTitle(Form("q_{%d}", i+1));
    fHistoReducedQvectors[i]->GetYaxis()->SetTitle(Form("dN/dq_{%d}", i+1));
    fMPCList->Add(fHistoReducedQvectors[i]);
  }   // Enf of the loop for the q_n distributions.

// Profile with the 2-particle correlators.
  fProfileTwoPartCorrel = new TProfile("fProfileTwoPartCorrel",
    "#LT#LT2#GT#GT_{n,-n} (no #eta gap)", 8, 0., 8.);
  fProfileTwoPartCorrel->Sumw2();
  fProfileTwoPartCorrel->SetStats(kTRUE);
  fProfileTwoPartCorrel->GetXaxis()->SetTitle("n");
  fMPCList->Add(fProfileTwoPartCorrel);

// Profiles with the 4-particle correlators.
  fProfileFourPartCorrel = new TProfile("fProfileFourPartCorrel",
    "#LT#LT4#GT#GT_{m,n,-m,-n}", 36, 0., 36.);
  fProfileFourPartCorrel->Sumw2();
  fProfileFourPartCorrel->SetStats(kTRUE);
  fProfileFourPartCorrel->GetXaxis()->SetTitle("(m,n)");
  fMPCList->Add(fProfileFourPartCorrel);

  fProfileFourPartCorrelCheck = new TProfile("fProfileFourPartCorrelCheck",
    "#LT#LT4#GT#GT_{m,n,-m,-n}", 28, 0., 28.);
  fProfileFourPartCorrelCheck->Sumw2();
  fProfileFourPartCorrelCheck->SetStats(kTRUE);
  fProfileFourPartCorrelCheck->GetXaxis()->SetTitle("(m,n)");
  fMPCList->Add(fProfileFourPartCorrelCheck);

// Profile with the needed 6-particle correlators.
  fProfileSixPartCorrel = new TProfile("fProfileSixPartCorrel",
    "#LT#LT6#GT#GT_{m,n,o,-m,-n,-o}", 56, 0., 56.);
  fProfileSixPartCorrel->Sumw2();
  fProfileSixPartCorrel->SetStats(kTRUE);
  //fProfileSixPartCorrel->GetXaxis()->SetTitle("(m,n,o)");
  fMPCList->Add(fProfileSixPartCorrel);
}   // End of "void BookMPCList()".

/* ------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookTPCEtaList()
{
/* Book the profiles for the 2-particle correlators with eta gaps. */
  if (fComputeEtaGaps)
  {
    Float_t etaGaps[11] = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};  // Possible values for the eta gap.
    for (Int_t i = 0; i < 11; i++)
    {
      fProfileTPCEta[i] = new TProfile("", "", 8, 0., 8.);
      fProfileTPCEta[i]->SetName(Form("fProfileTPCEta_EtaGap%d", i));
      fProfileTPCEta[i]->SetTitle(Form("#LT#LT2#GT#GT_{n,-n} (#eta gap: %.1f)", etaGaps[i]));
      fProfileTPCEta[i]->SetStats(kTRUE);
      fProfileTPCEta[i]->Sumw2();
      fProfileTPCEta[i]->GetXaxis()->SetTitle("n");
      fTPCEtaList->Add(fProfileTPCEta[i]);

    // Set bin labels.
      for (Int_t gap = 1; gap < 9; gap++)
      {
       fProfileTPCEta[i]->GetXaxis()->SetBinLabel(gap, Form("%d", gap));
      }   // End of the loop to set the labels.
    }   // End of the loop to setup the profiles.
  }   // End of the if.
}   // End of "void BookTPCEtaList()".



