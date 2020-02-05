
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
#include <TExMap.h>
#include "TDirectoryFile.h"
#include "TSystem.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskTwoMultiCorrelations)

/* ========================================================================================== /
/ Mandatory methods needed for AliAnalysisTaskSE.                                             /
/ ========================================================================================== */
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations() :
  AliAnalysisTaskSE(),
// 1. General parameters for the configuration of the analysis.
  fMainList(NULL),
  fHighestHarmonic(8), fLargestCorrelators(8),
  fDoKineAnalysis(kFALSE),
  fDoCorrelationsHisto(kFALSE),
  fComputeEtaGaps(kFALSE),
  fRecoEvent(NULL), fKineEvent(NULL),
// 2. Parameters related to the centrality.
  fMultiplicityList(NULL),
  fCentralityMin(0.), fCentralityMax(100.),
  fCentralityFromSPD(kFALSE), fCentralityFromVZero(kFALSE),
  fCentrality(-1.),
  fHistoEvents(NULL),
// 3. Parameters related to the multiplicity.
  fInitialMultiplicity(0),
// 4. Parameters related to the physics event selection.
  fEventSelectionList(NULL),
  fCutOnPVx(kFALSE), fPVxMin(-44.), fPVxMax(44.),
  fCutOnPVy(kFALSE), fPVyMin(-44.), fPVyMax(44.),
  fCutOnPVz(kFALSE), fPVzMin(-10.), fPVzMax(-10.),
// 5. Parameters related to the HMOs selection.
  fMultiplicityMin(6),
  fMainFilter(128), fGlobalFilter(256),
  fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.), fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.), fMultiplicityMaxB(0.),
// 6. Parameters related to the track selection.
  fTrackSelectionList(NULL),
  fCutOnPt(kFALSE), fPtMin(0.2), fPtMax(5.),
  fCutOnEta(kFALSE), fEtaMin(-0.8), fEtaMax(0.8),
  fCutOnNTPC(kFALSE), fNTPCMin(70),
  fCutOnChi(kFALSE), fChiMin(0.1), fChiMax(4.),
  fCutOnNITS(kFALSE), fNITSMin(2),
  fCutOnDCAxy(kFALSE), fDCAxyMax(3.2),
  fCutOnDCAz(kFALSE), fDCAzMax(2.4),
  fCutOnCharge(kFALSE), fKeepPositiveCharges(kFALSE),
  fNumberBinsPt(1000), fNumberBinsEta(1000), fNumberBinsPhi(1000),
// 7. Parameters related to the efficiency and acceptance weights.
  fUseKineRecoTable(kFALSE),
  fUseNonUnitParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fUseLocalFiles(kFALSE),
  fPeriodUsedForWeight(""),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kFALSE),
  fFilterbitIndex(0),
  fHistoEfficiency(NULL), fHistoEffInverse(NULL),
// 8. Parameters related to the multi-particle correlations.
  fMPCList(NULL),
  fReducedQPower(0),
  fProfileTwoPartCorrel(NULL), fProfileFourPartCorrel(NULL),
  fProfileFourPartCorrelCheck(NULL), fProfileSixPartCorrel(NULL),
// 9. Parameters related to the 2-particle correlations with eta gaps.
  fTPCEtaList(NULL)
{
/* Dummy constructor of the class.                                                            /
/ 1. Initialise the arrays of data members.                                                  */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations()");

// 1. Initialise the arrays of data members.
  InitialiseArraysOfQvectors();
  InitialiseArraysOfHistos();
  InitialiseArraysOfTProfiles();

} // End: AliAnalysisTaskTwoMultiCorrelations().

/* ----------------------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
  AliAnalysisTaskSE(name),
// 1. General parameters for the configuration of the analysis.
  fMainList(NULL),
  fHighestHarmonic(8), fLargestCorrelators(8),
  fDoKineAnalysis(kFALSE),
  fDoCorrelationsHisto(kFALSE),
  fComputeEtaGaps(kFALSE),
  fRecoEvent(NULL), fKineEvent(NULL),
// 2. Parameters related to the centrality.
  fMultiplicityList(NULL),
  fCentralityMin(0.), fCentralityMax(100.),
  fCentralityFromSPD(kFALSE), fCentralityFromVZero(kFALSE),
  fCentrality(-1.),
  fHistoEvents(NULL),
// 3. Parameters related to the multiplicity.
  fInitialMultiplicity(0),
// 4. Parameters related to the physics event selection.
  fEventSelectionList(NULL),
  fCutOnPVx(kFALSE), fPVxMin(-44.), fPVxMax(44.),
  fCutOnPVy(kFALSE), fPVyMin(-44.), fPVyMax(44.),
  fCutOnPVz(kFALSE), fPVzMin(-10.), fPVzMax(-10.),
// 5. Parameters related to the HMOs selection.
  fMultiplicityMin(6),
  fMainFilter(128), fGlobalFilter(256),
  fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.), fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.), fMultiplicityMaxB(0.),
// 6. Parameters related to the track selection.
  fTrackSelectionList(NULL),
  fCutOnPt(kFALSE), fPtMin(0.2), fPtMax(5.),
  fCutOnEta(kFALSE), fEtaMin(-0.8), fEtaMax(0.8),
  fCutOnNTPC(kFALSE), fNTPCMin(70),
  fCutOnChi(kFALSE), fChiMin(0.1), fChiMax(4.),
  fCutOnNITS(kFALSE), fNITSMin(2),
  fCutOnDCAxy(kFALSE), fDCAxyMax(3.2),
  fCutOnDCAz(kFALSE), fDCAzMax(2.4),
  fCutOnCharge(kFALSE), fKeepPositiveCharges(kFALSE),
  fNumberBinsPt(1000), fNumberBinsEta(1000), fNumberBinsPhi(720),
// 7. Parameters related to the efficiency and acceptance weights.
  fUseKineRecoTable(kFALSE),
  fUseNonUnitParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fUseLocalFiles(kFALSE),
  fPeriodUsedForWeight(""),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kFALSE),
  fFilterbitIndex(0),
  fHistoEfficiency(NULL), fHistoEffInverse(NULL),
// 8. Parameters related to the multi-particle correlations.
  fMPCList(NULL),
  fReducedQPower(0),
  fProfileTwoPartCorrel(NULL), fProfileFourPartCorrel(NULL),
  fProfileFourPartCorrelCheck(NULL), fProfileSixPartCorrel(NULL),
// 9. Parameters related to the 2-particle correlations with eta gaps.
  fTPCEtaList(NULL)
{
/* Constructor of the class.                                                                  /
/ 1. Create the mother list. (The rights on everything it holds are given to the list.)       /
/ 2. Define the input and output slots.                                                       /
/ 3. Initialise the arrays of data members.                                                  */
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// 1. Create the mother list.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE);

// 2. Define the input and output slots.
  DefineOutput(1, TList::Class());

// 3. Initialise the arrays of data members.
  InitialiseArraysOfQvectors();
  InitialiseArraysOfHistos();
  InitialiseArraysOfTProfiles();

} // End: AliAnalysisTaskTwoMultiCorrelations(const char *, Bool_t)

/* ----------------------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
/* Destructor of the class.                                                                  */
  if (fMainList) { delete fMainList; }
} // End: ~AliAnalysisTaskTwoMultiCorrelations()

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Define the outputs of the task at the beginning of the analysis.                           /
/ 1. Avoid name clashes. (Part 1.)                                                            /
/ 2. JEfficiency for NUA correction : DongJo (If chosen).                                     /
/ 3. Book the lists and their content.                                                        /
/ 4. Avoid name clashes. (Part 2.)                                                           */
  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()";

// 1. Avoid name clashes.
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

// 2. JEfficiency for NUA correction : DongJo
  if (fUseJEfficiency)
  {
    fEfficiency = new AliJEfficiency();
    fEfficiency->SetMode(1);  // 1: priod should work for you
    fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data");
  }

// 3. Book the lists and their content.
  this->BookAllLists();
  this->BookMultiplicityList();
  this->BookEventSelectionList();
  this->BookTrackSelectionList();
  this->BookMPCList();
  this->BookTPCEtaList();

// 4. Avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1, fMainList);
  fFirstEvent = kTRUE;

} // End: void UserCreateOutputObjects().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the chosen analysis for each event of the input dataset.                           /
/ 1. Get the pointers for the current event at reco and kine levels. (The kine level exists   /
/   only for MC-Generator AOD files.)                                                         /
/ 2. Select the correct analysis according to the AddTask configuration.                      /
/ 3. PostData.                                                                               */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)";

// 1. Get the pointers for the current event at reco and kine levels.
  fRecoEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  fKineEvent = MCEvent();

// 2. Select the correct analysis according to the AddTask configuration.
  if (fDoKineAnalysis) { GetRatioDistributions(); }
  else { AnalyseRecoEvent(); }

// 3. PostData.
  PostData(1, fMainList);

} // End: void UserExec(Option_t *).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
/* Save the outputs at the end of the execution of the script.                                /
/ 1. Access the mother list.                                                                  /
/ 2. Create the output file and save inside the mother list.                                 */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)";
// 1. Access the mother list.
  fMainList = (TList*)GetOutputData(1);
  if (!fMainList) { Fatal(sMethod.Data(), "ERROR: fMainList not found."); }

// 2. Create the output file and save inside the mother list.
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fMainList->Write(fMainList->GetName(), TObject::kSingleKey);
  delete outputFile;

} // End: void Terminate(Option_t *).

/* ========================================================================================== /
/ Methods called in the constructors.                                                         /
/ ========================================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfQvectors()
{
/* Initialise to zero all the elements of "fQvectors".                                       */
  for (Int_t iHarmo = 0; iHarmo < 65; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 9; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    } // End: iPower.
  } // End: iHarmo.

} // End: void InitialiseArraysOfQvectors().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfHistos()
{
/* Initialise to NULL all the histograms inside arrays.                                      */
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

  for (Int_t i = 0; i < 3; i++) { fHistoWeights[i] = NULL; }

  for (Int_t i = 0; i < 8; i++) { fHistoReducedQvectors[i] = NULL; }

} // End: void InitialiseArraysOfHistos().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfTProfiles()
{
/* Initialise to NULL all the profiles inside arrays. */
  for (Int_t i = 0; i < 11; i++) {fProfileTPCEta[i] = NULL;}
} // End: void InitialiseArraysOfTProfiles().

/* ========================================================================================== /
/ Methods called in "UserExec".                                                               /
/ ========================================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::AnalyseRecoEvent()
{
/* Do the analysis at reco level to get the multiparticle correlations.                       /
/ 1. Check the presence of the pointer to an AOD track.                                       /
/ 2. Check if the event belongs to the current centrality bin.                                /
/ 3. Apply the event selection with the HMOs criteria.                                        /
/ 4. If chosen: get the run number and load the correct NUE correction.                       /
/ 5. Identify which tracks pass the selection.                                                /
/ 6. Keep the events with enough tracks to have a meaningful event weight.                    /
/ 7. Get the arrays needed to compute the multiparticle correlations.                         /
/ 8. Get the particle weights if chosen.                                                      /
/ 9. Compute all the possible multi-particle correlations.                                    /
/ 10. Reset the variables for the next event*/
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::AnalyseRecoEvent()";

// 1. Check the presence of the pointer to an AOD track.
  if (!fRecoEvent) { Fatal(sMethod.Data(), "ERROR: no AOD event found."); }

// 2. Check if the event belongs to the current centrality bin.
//    Save the number of events in this centrality bin.
  if (!ApplyCentralitySelection()) {return;}
  fHistoEvents->Fill(0.5);

// 3. Apply the event selection with the HMOs criteria.
//    Save the number of events after the event selection.
  if(!ApplyEventSelection(kTRUE)) {return;}
  fHistoEvents->Fill(1.5);

// 4. If chosen: get the run number and load the correct NUE correction.
  if (fUseJEfficiency && fFirstEvent)
  {
    fEfficiency->SetRunNumber(fRecoEvent->GetRunNumber());
    fEfficiency->Load();
    fFirstEvent = kFALSE;
  }

// 5. Identify which tracks pass the selection.
  long long multiplicity = 0;         // Current number of tracks in the event.
  long long finalMultiplicity = 0;    // Multiplicity after the track selection.

// 5.1 Get the current number of events.
  multiplicity = fRecoEvent->GetNumberOfTracks();
  Bool_t *isTrackSelected = new Bool_t[multiplicity](); // Flag if the track passed or not.
    // kTRUE: the track passed the selection, kFALSE: the track is rejected.

// 5.2
  for (long long iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // 5.2.1 Get the pointer to the current track.
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(fRecoEvent->GetTrack(iTrack));
    if (!aTrack) { continue; }
    if (!aTrack->TestFilterBit(fMainFilter)) { continue; }

  // 5.2.2 Apply the track selection at reco level.
    if (ApplyTrackSelection(aTrack))  // The track passed the selection.
    { finalMultiplicity++; isTrackSelected[iTrack] = kTRUE; }
    else { isTrackSelected[iTrack] = kFALSE; }  // The track failed the selection.
  } // End of the loop over the tracks.

// 5.3 Save the number of events after the track selection.
  fHistoEvents->Fill(2.5);

// 6. Keep the events with enough tracks to have a meaningful event weight.
//    Save the final number of tracks and final number of events.
  if (finalMultiplicity < fMultiplicityMin) { return; }
  fHistoMultiplicity[1]->Fill(finalMultiplicity);
  fHistoEvents->Fill(3.5);

// 7. Get the arrays needed to compute the multiparticle correlations.
  Int_t iIndex = 0;           // Index of the selected track in the final arrays.
  Float_t iEffCorr = 1.;      // Efficiency correction.
  Float_t iEffInverse = 1.;   // Inverse of the efficiency correction.
  Float_t *iPt  = new Float_t[finalMultiplicity](); // Transverse momentum.
  Float_t *iEta = new Float_t[finalMultiplicity](); // Pseudorapidity for the eta gaps.
  Float_t *iPhi = new Float_t[finalMultiplicity](); // Azimuthal angles.
  Float_t *iParticleWeights = new Float_t[finalMultiplicity](); // Particle weights.

  for (Int_t iiTrack = 0; iiTrack < multiplicity; iiTrack++)
  {
  // 7.1 Check if iiTrack corresponds to a selected track or not.
    if (!isTrackSelected[iiTrack]) { continue; }

  // 7.2 Get a pointer to the selected track.
    AliAODTrack *aaTrack = dynamic_cast<AliAODTrack*>(fRecoEvent->GetTrack(iiTrack));
    if (!aaTrack) { continue; }
    if (!aaTrack->TestFilterBit(fMainFilter)) { continue; }

  // 7.3 Get all the needed variables.
    iPt[iIndex] = aaTrack->Pt();
    iEta[iIndex] = aaTrack->Eta();
    iPhi[iIndex] = aaTrack->Phi();
    iParticleWeights[iIndex] = 1.;
    if (fUseJEfficiency)
    {
      iEffCorr = fEfficiency->GetCorrection(iPt[iIndex], fFilterbitIndex, fCentrality);
      iEffInverse = 1.0/iEffCorr;
      iParticleWeights[iIndex] = iEffInverse;
      fHistoEfficiency->Fill(iEffCorr);
      fHistoEffInverse->Fill(iEffInverse);
    }

  // 7.4 Increase the index in the observables' arrays.
    iIndex++;
  }   // End of the second loop over the tracks.

// 8. Get the particle weights if chosen.
  if (fUseNonUnitParticleWeights)
  {
    Int_t eventRunNumber = fRecoEvent->GetRunNumber();
    //printf("Run number: %d \n", eventRunNumber);
    CalculateWeight(eventRunNumber, finalMultiplicity, iParticleWeights, iPhi, iPt, iEta);
  }

// 9. Compute all the possible multi-particle correlations.
  CalculateQvectors(finalMultiplicity, iPhi, iParticleWeights);
  ComputeReducedQvectors(finalMultiplicity);
  ComputeAllCorrelators(finalMultiplicity, iPhi, iParticleWeights);
  if (fComputeEtaGaps)
  { ComputeTPCWithEtaGaps(finalMultiplicity, iPhi, iParticleWeights, iEta); }

// 10. Reset the variables for the next event.
  fInitialMultiplicity = 0;
  multiplicity = 0;
  finalMultiplicity = 0;
  delete [] isTrackSelected;
  iIndex = 0;
  iEffCorr = 0.;
  iEffInverse = 0.;
  delete [] iPt;
  delete [] iEta;
  delete [] iPhi;
  delete [] iParticleWeights;

} // End: void AnalyseRecoEvent().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::GetRatioDistributions()
{
/* Do the event/track selection at kine level to get the distributions.                       /
/ 1. Check the presence of the pointer to the reco and kine events.                           /
/ 2. Keep only the reco events belonging to the current centrality bin (The                   /
/   centrality estimators are based on detectors' informations, so it does not                /
/   work for kine events).                                                                    /
/ 3. Apply the event selection at kine level.                                                 /
/ 4. Get the mapping table between reco and kine levels.                                      /
/ 5. Keep the events with enough tracks to have a meaningful event weight.                    /
/ 6. Get the distributions at kine level.                                                     /
/ 7. Reset the variables for the next event.                                                 */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::GetRatioDistributions()";

// 1. Check the presence of the pointer to the reco and kine events.
  if (!fRecoEvent) { Fatal(sMethod.Data(), "ERROR: no AOD event found."); }
  if (!fKineEvent) { Fatal(sMethod.Data(), "ERROR: no MC event found."); }

// 2. Keep only the reco events belonging to the current centrality bin.
//    Save the number of events in this centrality bin.
  if (!ApplyCentralitySelection()) {return;}
  fHistoEvents->Fill(0.5);

// 3.1 Apply the event selection at reco level.
  if(!ApplyEventSelection(kTRUE)) {return;}

// 3.2 Apply the event selection at kine level.
//    Save the number of events after the event selection.
  if(!ApplyEventSelection(kFALSE)) {return;}
  fHistoEvents->Fill(1.5);

// 4. Get the mapping table between reco and kine levels.
  long long multiplicity = 0;           // Current number of reco tracks.
  long long finalMultiplicity = 0;      // Multiplicity after the track selection.
  TExMap *kineRecoMap = new TExMap();   // Mapping table between kine and reco tracks.
    // Key: label of the kine track ; value: selection status of the reco track.

// 4.1 Get the number of reco tracks.
  multiplicity = fRecoEvent->GetNumberOfTracks();

// 4.2
  for (long long iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // 4.2.1 Get the pointer to the current track.
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(fRecoEvent->GetTrack(iTrack));
    if (!aTrack) { continue; }
    if (!aTrack->TestFilterBit(fMainFilter)) { continue; }

  // 4.2.2 Apply the track selection at reco level.
    if (ApplyTrackSelection(aTrack))  // The track passed the selection: MC-key tagged with "1".
    {
      if (aTrack->GetLabel() >= 0)
      {
        if ( kineRecoMap->GetValue(aTrack->GetLabel()) == 0 ) {kineRecoMap->Add(aTrack->GetLabel(), 1);}
      }
      finalMultiplicity++;
    }
    else
    {
      if (aTrack->GetLabel() >= 0)
      {
        if ( kineRecoMap->GetValue(aTrack->GetLabel()) == 0 ) {kineRecoMap->Add(aTrack->GetLabel(), -1);}
      }
    }
      // The track failed the selection: MC-key tagged with "-1".
  } // End of the loop over the tracks.

// 5. Keep the events with enough tracks to have a meaningful event weight.
//    Save the final number of tracks and final number of events.
  if (finalMultiplicity < fMultiplicityMin) { return; }
  fHistoMultiplicity[1]->Fill(finalMultiplicity);
  fHistoEvents->Fill(3.5);

// 6. Get the distributions at kine level.
  multiplicity = 0;
  multiplicity = fKineEvent->GetNumberOfTracks();

  for (long long iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // 6.1 Get the pointer to the current track.
    AliAODMCParticle *aMCtrack = dynamic_cast<AliAODMCParticle*>(fKineEvent->GetTrack(iTrack));
    if (!aMCtrack) { continue; }

	  Float_t phi = aMCtrack->Phi();  // Azimuthal angle.
    Float_t eta = aMCtrack->Eta();  // Pseudorapidity.
    Float_t pt = aMCtrack->Pt();    // Tranverse momentum.
    Int_t   charge = aMCtrack->Charge();    // Charge.

  // 4.2.2 Apply the track selection at kine level.
    if (!ApplyTrackSelection(aMCtrack)) { continue; }
    if (fUseKineRecoTable)
	  {
		  if (kineRecoMap->GetValue(iTrack)==-1) {continue;}
        // GetValue = -1: the corresponding reco track has not been selected by its track selection.
		  else
		  {
      // GetValue = 0: the kine track has been lost in ALICE --> Included in the pT distribution.
      // GetValue = 1: the reco track has been selected.
        fHistoPt[1]->Fill(pt);
        fHistoEta[1]->Fill(eta);
        fHistoPhi[1]->Fill(phi);
        fHistoCharge[1]->Fill(charge);
		  }
	  }
	  else
	  {
        fHistoPt[1]->Fill(pt);
        fHistoEta[1]->Fill(eta);
        fHistoPhi[1]->Fill(phi);
        fHistoCharge[1]->Fill(charge);
	  }
  } // End of the loop over the tracks.

// 7. Reset the variables for the next event.
  fInitialMultiplicity = 0;
  multiplicity = 0;
  finalMultiplicity = 0;
  kineRecoMap->Delete();  // Empty the mapping table.
  delete kineRecoMap; // Delete the pointer to the mapping table.

} // End: void GetRatioDistributions().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyCentralitySelection()
{
/* Check if the reconstructed event belongs to the current centrality bin.                    /
/ 1. Set the centrality estimator according to user's settings.                               /
/ 2. Fill the initial centrality distribution.                                                /
/ 3. Apply the centrality selection.                                                          /
/ 4. Fill the final centrality distribution.                                                 */
  TString sMethod = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyCentralitySelection()";

// 1. Set the centrality estimator according to user's settings.
  TString centrEstimator = "centralityEstimator";
  if ( (Int_t)fCentralityFromVZero + (Int_t)fCentralityFromSPD != 1 )
  { Fatal(sMethod.Data(), "ERROR: only one estimator must be set to kTRUE."); }
  else if (fCentralityFromSPD) {centrEstimator = "CL1";}
  else if (fCentralityFromVZero) {centrEstimator = "V0M";}

// 2. Fill the initial centrality distribution.
  AliMultSelection *multSelection = (AliMultSelection*)fRecoEvent->FindListObject("MultSelection");
  if (!multSelection) {Fatal(sMethod.Data(), "ERROR: no multSelection found.");}
  fCentrality = multSelection->GetMultiplicityPercentile(Form("%s", centrEstimator.Data()));
  fHistoCentrality[0]->Fill(fCentrality);

// 3. Apply the centrality selection.
  if ( (fCentrality < fCentralityMin) || (fCentrality >= fCentralityMax) )
  { return kFALSE; }

// 4. Fill the final centrality distribution.
  fHistoCentrality[1]->Fill(fCentrality);

  return kTRUE;
} // End: Bool_t ApplyCentralitySelection().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelection(Bool_t isRecoEvent)
{
/* Apply the event selection to the current event.                                            /
/ 1. Get the right variables depending on the type of event.                                  /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the physics event selection.                                                       /
/ 4. Apply the HMOs selection (for reco only).                                                /
/ 5. Fill the final QA histograms.                                                            /
/ 6. Reset the variables for the next event.                                                 */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelection(Bool_t)";
  Float_t PVx = 0.;   // x-position of the Primary Vertex.
  Float_t PVy = 0.;   // y-position of the Primary Vertex.
  Float_t PVz = 0.;   // z-position of the Primary Vertex.

// 1. Get the right variables depending on the type of event.
  if (isRecoEvent)  // The event selection is done at reco level.
  {
    fInitialMultiplicity = fRecoEvent->GetNumberOfTracks();
    AliAODVertex *primaryVertex = (AliAODVertex*)fRecoEvent->GetPrimaryVertex();
    PVx = primaryVertex->GetX();
    PVy = primaryVertex->GetY();
    PVz = primaryVertex->GetZ();
  }
  else  // The event selection is done at kine level.
  {
    fInitialMultiplicity = fKineEvent->GetNumberOfTracks();
    AliMCVertex *primaryVertex = (AliMCVertex*)fKineEvent->GetPrimaryVertex();
    PVx = primaryVertex->GetX();
    PVy = primaryVertex->GetY();
    PVz = primaryVertex->GetZ();
  }

// 2. Fill the initial QA histograms.
//    Reset the histograms if the filling needs to be done at kine level.
  if (!isRecoEvent)
  {
    fHistoMultiplicity[0]->Reset();
    fHistoPVx[0]->Reset();
    fHistoPVy[0]->Reset();
    fHistoPVz[0]->Reset();
  }

  fHistoMultiplicity[0]->Fill(fInitialMultiplicity);
  fHistoPVx[0]->Fill(PVx);
  fHistoPVy[0]->Fill(PVy);
  fHistoPVz[0]->Fill(PVz);

// 3. Apply the physics criteria of the event selection.
  if (fCutOnPVx)
  { if ( (PVx < fPVxMin) || (PVx > fPVxMax) ) {return kFALSE;} }
  if (fCutOnPVy)
  { if ( (PVy < fPVyMin) || (PVy > fPVyMax) ) {return kFALSE;} }
  if (fCutOnPVz)
  { if ( (PVz < fPVzMin) || (PVz > fPVzMax) ) {return kFALSE;} }

// 4. Apply the HMOs selection (for reco only).
  if (isRecoEvent)
  { if (fCutOnHMOs) { if (!RemoveHMOs()) {return kFALSE;} } }

// 5. Fill the final QA histograms.
//    Reset the histogram for the filling at kine level.
  if (!isRecoEvent)
  {
    fHistoPVx[1]->Reset();
    fHistoPVy[1]->Reset();
    fHistoPVz[1]->Reset();
  }
  fHistoPVx[1]->Fill(PVx);
  fHistoPVy[1]->Fill(PVy);
  fHistoPVz[1]->Fill(PVz);

// 6. Reset the variables for the next event.
  PVx = 0.;
  PVy = 0.;
  PVz = 0.;

  return kTRUE;
} // End: Bool_t ApplyEventSelection(Bool_t).

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::RemoveHMOs()
{
/* Get and apply the criteria to remove the HMOs in AOD events.                               /
/ 1. Get the number of tracks after the physics event selection.                              /
/ 2.1 Get a pointer to a reconstructed track.                                                 /
/ 2.2 Get the number of tracks in the main and global filters.                                /
/ 3. Fill the initial QA histograms.                                                          /
/ 4. Apply the HMOs cuts.                                                                     /
/ 5. Fill the final QA histograms.                                                            /
/ 6. Reset the variables for the next event.                                                 */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::RemoveHMOs()";
  long long multiplicity = 0;         // Number of tracks after the selection.
  long long globalMultiplicity = 0;   // Number of tracks in the global filter.
  long long mainMultiplicity = 0;     // Number of tracks in the main filter.

// 1. Get the number of tracks after the physics event selection.
  multiplicity = fRecoEvent->GetNumberOfTracks();

// 2.
  for (long long iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // 2.1 Get a pointer to a reconstructed track.
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(fRecoEvent->GetTrack(iTrack));
    if (!aTrack) { continue; }

  // 2.2 Get the number of tracks in the main and global filters.
    if (aTrack->TestFilterBit(fGlobalFilter)) { globalMultiplicity++; }
    if (aTrack->TestFilterBit(fMainFilter)) { mainMultiplicity++; }
  } // End: iTrack.

// 3. Fill the initial QA histograms.
  fHistoMultiplicityGlobal[0]->Fill(globalMultiplicity);
  fHistoMultiplicityMain[0]->Fill(mainMultiplicity);
  if (fDoCorrelationsHisto)
  { fHistoCorrelatedFilters[0]->Fill(globalMultiplicity, mainMultiplicity); }

// 4. Apply the HMOs cuts.
  if (fCutOnHMOs)
  {
  // Minimum and maximum lines which delimite the accepted band.
    Float_t lineMin = fMultiplicityMinA*(Float_t)globalMultiplicity + fMultiplicityMinB;
    Float_t lineMax = fMultiplicityMaxA*(Float_t)globalMultiplicity + fMultiplicityMaxB;

  // The number of tracks in the main filter is under the accepted band.
    if ( (Float_t)mainMultiplicity < lineMin ) { return kFALSE; }
  // The number of tracks in the main filter is above the accepted band.
    if ( (Float_t)mainMultiplicity > lineMax ) { return kFALSE; }
  } // End: if (fCutOnHMOs).

// 5. Fill the final QA histograms.
  fHistoMultiplicityGlobal[1]->Fill(globalMultiplicity);
  fHistoMultiplicityMain[1]->Fill(mainMultiplicity);
  if (fDoCorrelationsHisto)
  { fHistoCorrelatedFilters[1]->Fill(globalMultiplicity, mainMultiplicity); }

// 6. Reset the variables for the next event.
  multiplicity = 0;
  globalMultiplicity = 0;
  mainMultiplicity = 0;

  return kTRUE;
} // End: Bool_t RemoveHMOs().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)
{
/* Apply the track selection criteria.                                                        /
/ 1. Declare all the observables for the reco track.                                          /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the track selection.                                                               /
/ 4. Fill the final QA histograms.                                                            /
/ 5. Reset the variables for the next track.                                                 */
  TString sMethod = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)";

// 1. Declare all the observables for the reco track.
  Float_t pT = aAODtrack->Pt();           // Transverse momentum.
  Float_t eta = aAODtrack->Eta();         // Pseudorapidity.
  Float_t phi = aAODtrack->Phi();         // Azimuthal angle.
  Int_t NTPC = aAODtrack->GetTPCNcls();   // Number of TPC clusters.
  Float_t chiSquare = (aAODtrack->GetTPCchi2())/(aAODtrack->GetNcls(1));  // Chi square of p in the TPC.
  Int_t NITS = aAODtrack->GetITSNcls();   // Number of ITS clusters.
  Float_t DCAx = aAODtrack->XAtDCA();     // DCA along x.
  Float_t DCAy = aAODtrack->YAtDCA();     // DCA along y.
  Float_t DCAz = aAODtrack->ZAtDCA();     // DCA along z.
  Int_t charge = aAODtrack->Charge();     // Electric charge.
  Float_t DCAxy = TMath::Sqrt( (DCAx*DCAx) + (DCAy*DCAy) ); // DCA in the xy-plane.

// 2. Fill the initial QA histograms.
//    Don't fill them if it is for computing the weights.
  if (!fDoKineAnalysis)
  {
    fHistoPt[0]->Fill(pT);
    fHistoEta[0]->Fill(eta);
    fHistoPhi[0]->Fill(phi);
    fHistoNTPC[0]->Fill(NTPC);
    fHistoChiSquare[0]->Fill(chiSquare);
    fHistoNITS[0]->Fill(NITS);
    fHistoDCAxy[0]->Fill(DCAxy);
    fHistoDCAz[0]->Fill(DCAz);
    fHistoCharge[0]->Fill(charge);
  }

// 3. Apply the track selection.
  if (fCutOnPt)
  { if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;} }
  if (fCutOnEta)
  { if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;} }
  if (fCutOnNTPC)
  { if (NTPC < fNTPCMin) {return kFALSE;} }
  if (fCutOnChi)
  { if ( (chiSquare < fChiMin) || (chiSquare > fChiMax) ) {return kFALSE;} }
  if (fCutOnNITS)
  {
    if (fMainFilter == 128) {Fatal(sMethod.Data(), "ERROR: TPConly tracks do not have ITS information.");}
    else if (NITS < fNITSMin) {return kFALSE;}
  }
  if (fCutOnDCAxy)
  { if (DCAxy > fDCAxyMax) {return kFALSE;} }
  if (fCutOnDCAz)
  { if (DCAz > fDCAzMax) {return kFALSE;} }
  if (fCutOnCharge)
  {
  // Keep only the positive tracks.
    if (fKeepPositiveCharges) { if (charge < 0) {return kFALSE;} }
  // Keep only the negative tracks.
    else if (!fKeepPositiveCharges) { if (charge > 0) {return kFALSE;} }
    else { Fatal(sMethod.Data(), "ERROR: Select the sign of the charges to keep."); }
  }

// 4. Fill the final QA histograms.
//    Don't fill them if it is for computing the weights.
  if (!fDoKineAnalysis)
  {
    fHistoPt[1]->Fill(pT);
    fHistoEta[1]->Fill(eta);
    fHistoPhi[1]->Fill(phi);
    fHistoNTPC[1]->Fill(NTPC);
    fHistoChiSquare[1]->Fill(chiSquare);
    fHistoNITS[1]->Fill(NITS);
    fHistoDCAxy[1]->Fill(DCAxy);
    fHistoDCAz[1]->Fill(DCAz);
    fHistoCharge[1]->Fill(charge);
  } // End: if (!fDoKineAnalysis).

// 5. Reset the variables for the next track.
  pT = 0.;
  eta = 0.;
  phi = 0.;
  NTPC = 0;
  chiSquare = 0.;
  NITS = 0;
  DCAx = 0.;
  DCAy = 0.;
  DCAz = 0.;
  charge = 0;
  DCAxy = 0.;

  return kTRUE;
} // End: Bool_t ApplyTrackSelection().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODMCParticle *aMCtrack)
{
/* Apply the track selection criteria to the MC particle.                                     /
/ 1. Declare all the observables for the kine track.                                          /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the track selection.                                                               /
/ 4. Reset the variables for the next track.                                                 */
  TString sMethodName = "Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODMCParticle *aMCtrack)";

// 1. Declare all the observables for the kine track.
  Float_t pT = aMCtrack->Pt();        // Transverse momentum.
  Float_t eta = aMCtrack->Eta();      // Pseudorapidity.
  Float_t phi = aMCtrack->Phi();      // Azimuthal angle.
  Int_t charge = aMCtrack->Charge();  // Electric charge.
  Bool_t isPrimary = aMCtrack->IsPhysicalPrimary();
    // kTRUE: the particle is a primary particle.
	Bool_t isWeakSecondary = aMCtrack->IsSecondaryFromWeakDecay();
  //Bool_t isMaterialSecondary = aMCtrack->IsSecondaryFromMaterial();

// 2. Fill the initial QA histograms.
  fHistoPt[0]->Fill(pT);
  fHistoEta[0]->Fill(eta);
  fHistoPhi[0]->Fill(phi);
  fHistoCharge[0]->Fill(charge);

// 3. Apply the track selection.
  if (charge == 0) {return kFALSE;}
  if (!isPrimary && !isWeakSecondary) {return kFALSE;}
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

// 4. Reset the variables for the next track.
  pT        = 0.;
  eta       = 0.;
  phi       = 0.;
  charge    = 0;

  return kTRUE;
} // End: Bool_t ApplyTrackSelection(AliAODMCParticle *).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::CalculateWeight(Int_t RunNumber, long long numberOfParticles, Float_t* pWeights, Float_t* angles, Float_t* pt, Float_t* eta)
{
/* Get the particle weights for all the tracks in the current event.                          /
/ 1. Get the correct histograms.                                                              /
/ 2. Get the phi, pT and eta weights.                                                         /
/ 3. Fill "pWeights".                                                                        */

// 1. Get the correct histograms.
  if (fUsePhiWeights) { fHistoWeights[0] = GetHistogramWithWeights(RunNumber, "phi"); }
  if (fUsePtWeights) { fHistoWeights[1] = GetHistogramWithWeights(RunNumber, "pt"); }
  if (fUseEtaWeights) { fHistoWeights[2] = GetHistogramWithWeights(RunNumber, "eta"); }

//
  for (long long i = 0; i < numberOfParticles; i++)
  {
    Int_t   iBin = 0;         // Current bin.
    Float_t weight_phi = 1.;  // phi weight, unit by default.
    Float_t weight_pt = 1.;   // pT weight, unit by default.
    Float_t weight_eta = 1.;  // eta weight, unit by default.

  // 2. Get the phi, pT and eta weight.
    if (fUsePhiWeights)
    {
      iBin = fHistoWeights[0]->FindBin(angles[i]);
      weight_phi = fHistoWeights[0]->GetBinContent(iBin);
    }
    if (fUsePtWeights)
    {
      iBin = fHistoWeights[1]->FindBin(pt[i]);
      weight_pt = fHistoWeights[1]->GetBinContent(iBin);
    }
    if (fUseEtaWeights)
    {
      iBin = fHistoWeights[2]->FindBin(eta[i]);
      weight_eta = fHistoWeights[2]->GetBinContent(eta[i]);
    }

  // 3. Fill "pWeights".
    pWeights[i] = weight_phi*weight_pt*weight_eta;
  } //

} // End of "void CalculateWeight(Int_t, long long, Float_t*, Float_t*, Float_t*, Float_t*)".

/* ----------------------------------------------------------------------------------------- */
TH1F* AliAnalysisTaskTwoMultiCorrelations::GetHistogramWithWeights(Int_t RunNumber, const char *variable)
{
/* Access the external ROOT file with the needed weights.                                     /
/ 1. Insanity check for the arguments.                                                        /
/ 2. Access the external ROOT file to fetch the desired histogram with weights.               /
/ 3. Close the external ROOT file.                                                           */
  TString sMethod = "void AliAnalysisTaskTwoMultiCorrelations::GetHistogramWithWeights()";
  TH1F* hist = NULL;

// 1. Insanity check for the arguments.
  if ( !(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")) )
  { Fatal(sMethod.Data(), "ERROR: not the correct variable."); }

// 2. Access the external ROOT file to fetch the desired histogram with weights.
  TFile *weightsFile =  NULL;
  if (fUseLocalFiles) { weightsFile = TFile::Open(Form("/home/cindy/TestAliAnalysisTask/Test_UseWeights/%s/%d/Weights.root", fPeriodUsedForWeight.Data(), RunNumber), "READ"); }
  else { weightsFile = TFile::Open(Form("/alice/cern.ch/user/c/cimordas/Weights/%s/%d/Weights.root", fPeriodUsedForWeight.Data(), RunNumber), "READ"); }

  if (!weightsFile) { Fatal(sMethod.Data(), "ERROR 404: File not found"); } 
 
  TDirectoryFile *weightsDirectoryFile = dynamic_cast<TDirectoryFile*>(weightsFile->Get(Form("%s_Weights", variable)));
  if (!weightsDirectoryFile) { Fatal(sMethod.Data(), "ERROR: Directory not found"); }

  TList *weightsList = dynamic_cast<TList*>(weightsDirectoryFile->Get(Form("%s_Weight=>%.1f-%.1f", variable, fCentralityMin, fCentralityMax))); 
  if (!weightsList) { Fatal(sMethod.Data(), "ERROR: List not found"); } 

  hist = dynamic_cast<TH1F*>(weightsList->FindObject(Form("%s_Weight_Hist", variable)));
  if (!hist) { Fatal(sMethod.Data(), "ERROR: Hist not found"); }
  else { hist->SetDirectory(0); } // Kill the default ownership.

// 3. Close the external ROOT file.
  weightsFile->Close();
  delete weightsFile;

  return hist;
} // End: void GetHistogramWithWeights(Int_t, const char*).

/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
/* Simplify the use of the Q-vectors in the next methods. */
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]);   // Use that Q*(n,p) = Q(-n,p).
}   // End of "TComplex Q()".

/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
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

/* ========================================================================================== /
/ Methods called in "Terminate".                                              /
/ ========================================================================================== */

/* ========================================================================================== /
/ Methods called in "UserCreateOutputObjects".                                /
/ ========================================================================================== */
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

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMultiplicityList()
{
/* Book the histograms related to the multiplicity distributions. */
// Number of events at each selection step.
  fHistoEvents = new TH1I("", "", 4, 0., 4.);
  fHistoEvents->SetName("fHistoEvents");
  fHistoEvents->SetTitle("Number of events at each step");
  fHistoEvents->SetStats(kTRUE);
  fHistoEvents->GetXaxis()->SetBinLabel(1, "After centrality selection");
  fHistoEvents->GetXaxis()->SetBinLabel(2, "After event selection");
  fHistoEvents->GetXaxis()->SetBinLabel(3, "After track selection");
  fHistoEvents->GetXaxis()->SetBinLabel(4, "Final number of events");
  fHistoEvents->GetYaxis()->SetTitle("Number of events");
  fMultiplicityList->Add(fHistoEvents);

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

/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookTrackSelectionList()
{
/* Book the control histograms for the track selection criteria. */
  for (Int_t i = 0; i < 2; i++)
  {
    TString iState = "Final";
    if (i == 0) {iState = "Initial";}

  // Transverse momentum.
    fHistoPt[i] = new TH1F("", "", fNumberBinsPt, 0., 20.);
    fHistoPt[i]->SetName(Form("fHistoPt%s", iState.Data()));
    fHistoPt[i]->SetTitle(Form("%s distribution of p_{T}", iState.Data()));
    fHistoPt[i]->SetStats(kTRUE);
    fHistoPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistoPt[i]->GetYaxis()->SetTitle("dN/dp_{T}");
    fTrackSelectionList->Add(fHistoPt[i]);

  // Pseudorapidity.
    fHistoEta[i] = new TH1F("", "", fNumberBinsEta, -1., 1.);
    fHistoEta[i]->SetName(Form("fHistoEta%s", iState.Data()));
    fHistoEta[i]->SetTitle(Form("%s distribution of #eta", iState.Data()));
    fHistoEta[i]->SetStats(kTRUE);
    fHistoEta[i]->GetXaxis()->SetTitle("#eta");
    fHistoEta[i]->GetYaxis()->SetTitle("dN/d#eta");
    fTrackSelectionList->Add(fHistoEta[i]);

  // Azimuthal angles.
    fHistoPhi[i] = new TH1F("", "", fNumberBinsPhi, 0., 6.3);
    fHistoPhi[i]->SetName(Form("fHistoPhi%s", iState.Data()));
    fHistoPhi[i]->SetTitle(Form("%s distribution of #varphi", iState.Data()));
    fHistoPhi[i]->SetStats(kTRUE);
    fHistoPhi[i]->GetXaxis()->SetTitle("#varphi");
    fHistoPhi[i]->GetYaxis()->SetTitle("dN/d#varphi");
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
/* ----------------------------------------------------------------------------------------- */
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

/* ----------------------------------------------------------------------------------------- */
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



