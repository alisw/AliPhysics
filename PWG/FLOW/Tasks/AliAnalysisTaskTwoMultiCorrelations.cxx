
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
  fACHarmoOne(-1), fACHarmoTwo(-1), fACHarmoThree(-1),
  fDoKineAnalysis(kFALSE), fComputeACs(kFALSE), fWriteMinimum(kFALSE),
  fMainList(NULL),
  fRecoEvent(NULL), fKineEvent(NULL),
  fProfileEventCuts(NULL), fProfileTrackCuts(NULL),
// 2. Parameters related to the centrality.
  fNumberBinsMulti(30000), fTotalCentralityBin(9),
  fCentralityBin(-1), fInitialMultiplicity(0),
  fCentralityMin(0.), fCentralityMax(100.),
  fFstCentrality(-1.), fSndCentrality(-1.),
  fGetEstimCorrel(kFALSE),
  fFstCentralityEstim(""), fSndCentralityEstim(""),
  fMinimumQAList(NULL), fMultiplicityList(NULL),
  fHistoInitCorrelEstim(NULL),
// 3. Parameters related to the physics event selection.
  fCutOnPVx(kFALSE), fPVxMin(-44.), fPVxMax(44.),
  fCutOnPVy(kFALSE), fPVyMin(-44.), fPVyMax(44.),
  fCutOnPVz(kFALSE), fPVzMin(-10.), fPVzMax(-10.),
  fNumberBinsPV(1000), fMaxHistoPV(20.),
  fEventQAList(NULL),
// 4. Parameters related to the HMOs selection.
  fMultiplicityMin(6),
  fFstFilter(128), fSndFilter(256),
  fGetFiltersCorrel(kFALSE), fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.), fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.), fMultiplicityMaxB(0.),
// 5. Parameters related to the track selection.
  fCutOnPt(kFALSE), fPtMin(0.2), fPtMax(5.),
  fCutOnEta(kFALSE), fEtaMin(-0.8), fEtaMax(0.8),
  fCutOnNTPC(kFALSE), fNTPCMin(70),
  fPersoChiSquare(1), fCutOnChi(kFALSE), fChiMin(0.1), fChiMax(4.),
  fCutOnNITS(kFALSE), fNITSMin(2),
  fCutOnDCAxy(kFALSE), fDCAxyMax(2.4),
  fCutOnDCAz(kFALSE), fDCAzMax(3.2),
  fCutOnCharge(kFALSE), fKeepPosCharges(kFALSE),
  fKeepWeakSecondaries(kTRUE),
  fNumberBinsPt(1000), fNumberBinsEta(1000), fNumberBinsPhi(720),
  fTrackQAList(NULL),
// 6. Parameters related to the efficiency and acceptance weights.
  fUseKineRecoTable(kFALSE), fUseParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fNumberRuns(90),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kFALSE),
  fFilterbitIndex(0),
  fHistoEfficiency(NULL), fHistoEffInverse(NULL),
// 7. Parameters related to the multi-particle correlations.
  fHighestHarmonic(8), fLargestCorrelators(10),
  fReducedQPower(0),
  fMPCList(NULL),
// 8. Parameters related to the 2-particle correlations with eta gaps.
  fComputeEtaGaps(kFALSE),
  fTPCEtaList(NULL)
{
/* Dummy constructor of the class............................................................ /
/ 1. Initialise the arrays of data members.                                                  */
  AliDebug(2, "Ali-TwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations()");

// 1. Initialise the arrays of data members.
  InitialiseArraysOfDataMembers();

} // End: AliAnalysisTaskTwoMultiCorrelations().

/* ----------------------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
  AliAnalysisTaskSE(name),
// 1. General parameters for the configuration of the analysis.
  fACHarmoOne(-1), fACHarmoTwo(-1), fACHarmoThree(-1),
  fDoKineAnalysis(kFALSE), fComputeACs(kFALSE), fWriteMinimum(kFALSE),
  fMainList(NULL),
  fRecoEvent(NULL), fKineEvent(NULL),
  fProfileEventCuts(NULL), fProfileTrackCuts(NULL),
// 2. Parameters related to the centrality.
  fNumberBinsMulti(30000), fTotalCentralityBin(9),
  fCentralityBin(-1), fInitialMultiplicity(0),
  fCentralityMin(0.), fCentralityMax(100.),
  fFstCentrality(-1.), fSndCentrality(-1.),
  fGetEstimCorrel(kFALSE),
  fFstCentralityEstim(""), fSndCentralityEstim(""),
  fMinimumQAList(NULL), fMultiplicityList(NULL),
  fHistoInitCorrelEstim(NULL),
// 3. Parameters related to the physics event selection.
  fCutOnPVx(kFALSE), fPVxMin(-44.), fPVxMax(44.),
  fCutOnPVy(kFALSE), fPVyMin(-44.), fPVyMax(44.),
  fCutOnPVz(kFALSE), fPVzMin(-10.), fPVzMax(-10.),
  fNumberBinsPV(1000), fMaxHistoPV(20.),
  fEventQAList(NULL),
// 4. Parameters related to the HMOs selection.
  fMultiplicityMin(6),
  fFstFilter(128), fSndFilter(256),
  fGetFiltersCorrel(kFALSE), fCutOnHMOs(kFALSE),
  fMultiplicityMinA(0.), fMultiplicityMinB(0.),
  fMultiplicityMaxA(0.), fMultiplicityMaxB(0.),
// 5. Parameters related to the track selection.
  fCutOnPt(kFALSE), fPtMin(0.2), fPtMax(5.),
  fCutOnEta(kFALSE), fEtaMin(-0.8), fEtaMax(0.8),
  fCutOnNTPC(kFALSE), fNTPCMin(70),
  fPersoChiSquare(1), fCutOnChi(kFALSE), fChiMin(0.1), fChiMax(4.),
  fCutOnNITS(kFALSE), fNITSMin(2),
  fCutOnDCAxy(kFALSE), fDCAxyMax(2.4),
  fCutOnDCAz(kFALSE), fDCAzMax(3.2),
  fCutOnCharge(kFALSE), fKeepPosCharges(kFALSE),
  fKeepWeakSecondaries(kTRUE),
  fNumberBinsPt(1000), fNumberBinsEta(1000), fNumberBinsPhi(720),
  fTrackQAList(NULL),
// 6. Parameters related to the efficiency and acceptance weights.
  fUseKineRecoTable(kFALSE), fUseParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fNumberRuns(90),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kFALSE),
  fFilterbitIndex(0),
  fHistoEfficiency(NULL), fHistoEffInverse(NULL),
// 7. Parameters related to the multi-particle correlations.
  fHighestHarmonic(8), fLargestCorrelators(10),
  fReducedQPower(0),
  fMPCList(NULL),
// 8. Parameters related to the 2-particle correlations with eta gaps.
  fComputeEtaGaps(kFALSE),
  fTPCEtaList(NULL)
{
/* Constructor of the class.................................................................. /
/ 1. Create the mother list. (The rights on everything it holds are given to the list.)       /
/ 2. Define the input and output slots.                                                       /
/ 3. Initialise the arrays of data members.                                                  */
  AliDebug(2, "Ali-TwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// 1. Create the mother list.
  fMainList = new TList();
  fMainList->SetName("outputAnalysis");
  fMainList->SetOwner(kTRUE);

// 2. Define the input and output slots.
  DefineOutput(1, TList::Class());

// 3. Initialise the arrays of data members.
  InitialiseArraysOfDataMembers();

} // End: AliAnalysisTaskTwoMultiCorrelations(const char *, Bool_t)

/* ----------------------------------------------------------------------------------------- */
AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
/* Destructor of the class. */
  if (fMainList) {delete fMainList;}
} // End: ~AliAnalysisTaskTwoMultiCorrelations()

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Define the outputs of the task at the beginning of the analysis........................... /
/ 1. Avoid name clashes. (Part 1.)                                                            /
/ 2. JEfficiency for NUA correction : DongJo (If chosen).                                     /
/ 3. Book the lists and their content.                                                        /
/ 4. 
/ 5. Avoid name clashes. (Part 2.)                                                           */
  TString sMethodName = "void Ali-TwoMultiCorrelations::UserCreateOutputObjects()";

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
  this->BookMinimumQAList();
  if (!fWriteMinimum) // All the QA histograms must be saved.
  {
    this->BookMultiplicityList();
    this->BookEventQAList();
    if (!fGetFiltersCorrel) // We keep only the multiplicity and event histograms in analysis of HMOs.
    {this->BookTrackQAList();}
  } // End: if (!fWriteMinimum).

  if (!fGetFiltersCorrel)
  {
    this->BookMPCList();
    if (fComputeEtaGaps) {this->BookTPCEtaList();}
  }

// 4. Fill the configuration profiles for the event and track selections.
  if (fFstCentralityEstim == "CL1") {fProfileEventCuts->Fill(0.5, 1);}
  else {fProfileEventCuts->Fill(0.5, 2);} // 1: CL1, 2: V0M
  if (fSndCentralityEstim == "CL1") {fProfileEventCuts->Fill(1.5, 1);}
  else {fProfileEventCuts->Fill(1.5, 2);}
  if (fCutOnPVx) {fProfileEventCuts->Fill(2.5, fPVxMin); fProfileEventCuts->Fill(3.5, fPVxMax);}
  if (fCutOnPVy) {fProfileEventCuts->Fill(4.5, fPVyMin); fProfileEventCuts->Fill(5.5, fPVyMax);}
  if (fCutOnPVz) {fProfileEventCuts->Fill(6.5, fPVzMin); fProfileEventCuts->Fill(7.5, fPVzMax);}
  fProfileEventCuts->Fill(8.5, fMultiplicityMin);
  fProfileEventCuts->Fill(9.5, fFstFilter);
  fProfileEventCuts->Fill(10.5, fSndFilter);
  if (fCutOnHMOs)
  {
    fProfileEventCuts->Fill(11.5, fMultiplicityMinA);
    fProfileEventCuts->Fill(12.5, fMultiplicityMinB);
    fProfileEventCuts->Fill(13.5, fMultiplicityMaxA);
    fProfileEventCuts->Fill(14.5, fMultiplicityMaxB);
  }

  if (fCutOnPt) {fProfileTrackCuts->Fill(0.5, fPtMin); fProfileTrackCuts->Fill(1.5, fPtMax);}
  if (fCutOnEta) {fProfileTrackCuts->Fill(2.5, fEtaMin); fProfileTrackCuts->Fill(3.5, fEtaMax);}
  if (fCutOnNTPC) {fProfileTrackCuts->Fill(4.5, fNTPCMin);}
  fProfileTrackCuts->Fill(5.5, fPersoChiSquare);
  if (fCutOnChi) {fProfileTrackCuts->Fill(6.5, fChiMin); fProfileTrackCuts->Fill(7.5, fChiMax);}
  if (fCutOnNITS) {fProfileTrackCuts->Fill(8.5, fNITSMin);}
  if (fCutOnDCAxy) {fProfileTrackCuts->Fill(9.5, fDCAxyMax);}
  if (fCutOnDCAz) {fProfileTrackCuts->Fill(10.5, fDCAzMax);}
  if (fCutOnCharge)
  {
    if (fKeepPosCharges) {fProfileTrackCuts->Fill(11.5, 1);}
    else {fProfileTrackCuts->Fill(11.5, -1);}
  }
  if (fUsePtWeights) {fProfileTrackCuts->Fill(12.5, 1);}
  if (fUseEtaWeights) {fProfileTrackCuts->Fill(13.5, 1);}
  if (fUsePhiWeights) {fProfileTrackCuts->Fill(14.5, 1);}

// 5. Avoid name clashes.
  TH1::AddDirectory(oldHistAddStatus);
  PostData(1, fMainList);
  fFirstEvent = kTRUE;

} // End: void UserCreateOutputObjects().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the chosen analysis for each event of the input dataset........................... /
/ 1. Get the pointers for the current event at reco and kine levels. (The kine level exists   /
/   only for MC-Generator AOD files.)                                                         /
/ 2. Select the correct analysis according to the AddTask configuration.                      /
/ 3. PostData.                                                                               */
  TString sMethod = "void Ali-TwoMultiCorrelations::UserExec(Option_t *)";

// 1. Get the pointers for the current event at reco and kine levels.
  fRecoEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  fKineEvent = MCEvent();

// 2. Select the correct analysis according to the AddTask configuration.
  if (fDoKineAnalysis) {GetRatioDistributions();}
  else {AnalyseRecoEvent();}

// 3. PostData.
  PostData(1, fMainList);

} // End: void UserExec(Option_t *).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
/* Save the outputs at the end of the execution of the script................................ /
/ 1. Access the mother list.                                                                  /
/ 2. Create the output file and save inside the mother list.                                 */
  TString sMethod = "void Ali-TwoMultiCorrelations::Terminate(Option_t *)";

// 1. Access the mother list.
  fMainList = (TList*)GetOutputData(1);
  if (!fMainList) {Fatal(sMethod.Data(), "ERROR: fMainList not found.");}

} // End: void Terminate(Option_t *).

/* ========================================================================================== /
/ Methods called in the constructors and setters.                                             /
/ ========================================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::SetListOfRuns(TString dataPeriod)
{
/* Set the list of runs to use according to the chosen data-taking period................... */
  TString sMethod = "void Ali-TwoMultiCorrelations::SetListOfRuns()";
// TODO: extend to Run 2 run lists.

  if (dataPeriod == "LHC10h")
  {
    fNumberRuns = 90;
    Int_t listRuns[90] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
    for (Int_t i = 0; i < fNumberRuns; i++) {fListRuns[i] = listRuns[i];}
  } // End: if (dataPeriod == "LHC10h").
  else {Fatal(sMethod.Data(), "FATAL: not a valid data period!");}

} // End: void SetListOfRuns(TString).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::SetInputParticleWeights(TString fileWeight)
{
/* Setter to open the external file with the particle weights and import them in the task.... /
/ 1. Open the external file.                                                                  /
/ 2. Parse the runs.                                                                          /
/ 2.1 Open the TDirectoryFile for the current run.                                            /
/ 2.2 Open the list for the current centrality range.                                         /
/ 2.3 Fill the pT-weight histogram if needed.                                                 /
/ 2.4 Fill the eta-weight histogram if needed.                                                /
/ 2.5 Fill the phi-weight histogram if needed.                                                /
/ 3. Close the external file.                                                                */
  TString sMethod = "void Ali-TwoMultiCorrelations::SetInputParticleWeights()";
  Float_t centralityArray[10] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80. };

// 1. Open the external file.
  TFile *weightsFile = TFile::Open(Form("%s", fileWeight.Data()), "READ");
  if (!weightsFile) {Fatal(sMethod.Data(), "ERROR 404: File not found");}

// 2. Parse the runs.
  for (Int_t iRun = 0; iRun < fNumberRuns; iRun++)
  {
  // 2.1 Open the TDirectoryFile for the current run.
    Int_t runNumber = fListRuns[iRun];
    //printf("Run number: %d\n", runNumber);
    TDirectoryFile *runTDF = dynamic_cast<TDirectoryFile*>(weightsFile->Get(Form("%d", runNumber)));
    if (!runTDF) {Fatal(sMethod.Data(), "ERROR: Directory not found");}

    for (Int_t iCent = 0; iCent < fTotalCentralityBin; iCent++)
    {
    // 2.2 Open the list for the current centrality range.
      TList *centralityList = dynamic_cast<TList*>(runTDF->Get(Form("Centrality-%.1f-%.1f", centralityArray[iCent], centralityArray[iCent+1])));
      if (!centralityList) {Fatal(sMethod.Data(), "ERROR: List not found");}

    // 2.3 Fill the pT-weight histogram if needed.
      fHistoPtWeight[iRun][iCent] = dynamic_cast<TH1F*>(centralityList->FindObject("pt-weight"));
      if (!fHistoPtWeight[iRun][iCent]) {Fatal(sMethod.Data(), "ERROR: pt-weight histogram not found");}
      else {fHistoPtWeight[iRun][iCent]->SetDirectory(0);} // Kill the default ownership.

    // 2.4 Fill the eta-weight histogram if needed.
      fHistoEtaWeight[iRun][iCent] = dynamic_cast<TH1F*>(centralityList->FindObject("eta-weight"));
      if (!fHistoEtaWeight[iRun][iCent]) { Fatal(sMethod.Data(), "ERROR: eta-weight histogram not found"); }
      else {fHistoEtaWeight[iRun][iCent]->SetDirectory(0);}  // Kill the default ownership.

    // 2.5 Fill the phi-weight histogram if needed.
      fHistoPhiWeight[iRun][iCent] = dynamic_cast<TH1F*>(centralityList->FindObject("phi-weight"));
      if (!fHistoPhiWeight[iRun][iCent]) {Fatal(sMethod.Data(), "ERROR: phi-weight histogram not found");}
      else {fHistoPhiWeight[iRun][iCent]->SetDirectory(0);}  // Kill the default ownership.
    } // End: iCent.
  } // End: iRun.

// 3. Close the external file.
  weightsFile->Close();
  delete weightsFile;

} // End: void SetInputParticleWeights(TString).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfDataMembers()
{
/* Initialise to zero all the elements belonging to arrays of data members. ---------------- */
// Q-vectors.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 11; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    } // End: iPower.
  } // End: iHarmo.

// Run-by-run list and histograms.
  for (Int_t iRun = 0; iRun < 90; iRun++)
  {
    fListRuns[iRun] = 0;

    for (Int_t iCent = 0; iCent < 9; iCent++)
    {
      fHistoPtWeight[iRun][iCent] = NULL;
      fHistoEtaWeight[iRun][iCent] = NULL;
      fHistoPhiWeight[iRun][iCent] = NULL;
    } // End: iCent.
  } // End: iRun.

// Histograms.
  for (Int_t i = 0; i < 2; i++)
  {
    fHistoInitCentrality[i] = NULL;
    fHistoCorrelFilters[i] = NULL;
  }

  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    fHistoNumberEvents[iCent] = NULL;
    fHistoFstCentrality[iCent] = NULL;
    fHistoSndCentrality[iCent] = NULL;
    fHistoFinCorrelEstim[iCent] = NULL;

    for (Int_t iState = 0; iState < 2; iState++)
    {
      fHistoMultiplicity[iCent][iState] = NULL;
      fHistoPVx[iCent][iState] = NULL;
      fHistoPVy[iCent][iState] = NULL;
      fHistoPVz[iCent][iState] = NULL;
      fHistoMultiplicityMain[iCent][iState] = NULL;
      fHistoMultiplicityGlobal[iCent][iState] = NULL;
      fHistoPt[iCent][iState] = NULL;
      fHistoEta[iCent][iState] = NULL;
      fHistoPhi[iCent][iState] = NULL;
      fHistoNTPC[iCent][iState] = NULL;
      fHistoChiSquare[iCent][iState] = NULL;
      fHistoNITS[iCent][iState] = NULL;
      fHistoDCAxy[iCent][iState] = NULL;
      fHistoDCAz[iCent][iState] = NULL;
      fHistoCharge[iCent][iState] = NULL;
    } // End: iState.

    for (Int_t i = 0; i < 8; i++) {fHistoReducedQvectors[iCent][i] = NULL;}

  // Profiles.
    fProfileTwoPartCorrel[iCent] = NULL;
    fProfileFourPartCorrel[iCent] = NULL;
    fProfileFourPartCorrelCheck[iCent] = NULL;
    fProfileSixPartCorrel[iCent] = NULL;
    fProfileEightPartCorrel[iCent] = NULL;
    fProfileTenPartCorrel[iCent] = NULL;

    for (Int_t i = 0; i < 11; i++) {fProfileTPCEta[iCent][i] = NULL;}
  } // End: iCent.

} // End: void InitialiseArraysOfDataMembers().

/* ========================================================================================== /
/ Methods called in "UserExec".                                                               /
/ ========================================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::AnalyseRecoEvent()
{
/* Do the analysis at reco level to get the multiparticle correlations....................... /
/ 1. Check the presence of the pointer to an AOD event.                                       /
/ 2. Check if the event belongs to the current centrality bin.                                /
/ 3. Apply the event selection with the HMOs criteria.                                        /
/ 4. If chosen: get the run number and load the correct NUE correction.                       /
/ 5. Identify which tracks pass the selection.                                                /
/ 6. Keep the events with enough tracks to have a meaningful event weight.                    /
/ 7. Get the arrays needed to compute the multiparticle correlations.                         /
/ 8. Get the particle weights if chosen.                                                      /
/ 9. Compute all the possible multi-particle correlations.                                    /
/ 10. Reset the variables for the next event.                                                */
  TString sMethod = "void Ali-TwoMultiCorrelations::AnalyseRecoEvent()";

// 1. Check the presence of the pointer to an AOD event.
  if (!fRecoEvent) {Fatal(sMethod.Data(), "ERROR: no AOD event found.");}

// 2. Get the centrality of this event.
  if (!PassCentralitySelection()) {return;}
  fHistoNumberEvents[fCentralityBin]->Fill(0.5);

// 3. Apply the physics event selection and the HMO selection if chosen.
//    Save the number of events after the full event selection.
  if(!ApplyEventSelection(kTRUE)) {return;}
  if (fCutOnHMOs) { if (!RemoveHMOs()) {return;}}
  fHistoNumberEvents[fCentralityBin]->Fill(1.5);

// Stop the analysis here if the goal is only to get the correlations plots for HMOS.
  if (fGetFiltersCorrel) {return;}

// 4. If chosen: get the run number and load the correct NUE correction.
  if (fUseJEfficiency && fFirstEvent)
  {
    fEfficiency->SetRunNumber(fRecoEvent->GetRunNumber());
    fEfficiency->Load();
    fFirstEvent = kFALSE;
  } // End: if (fUseJEfficiency && fFirstEvent).

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
    if (!aTrack) {continue;}
    if (!aTrack->TestFilterBit(fFstFilter)) {continue;}

  // 5.2.2 Apply the track selection at reco level.
    if (ApplyTrackSelection(aTrack))  // The track passed the selection.
    {finalMultiplicity++; isTrackSelected[iTrack] = kTRUE;}
    else {isTrackSelected[iTrack] = kFALSE;}  // The track failed the selection.
  } // End: for (long long iTrack = 0; iTrack < multiplicity; iTrack++).

// 5.3 Save the number of events after the track selection.
  fHistoNumberEvents[fCentralityBin]->Fill(2.5);

// 6. Keep the events with enough tracks to have a meaningful event weight.
//    Save the final number of tracks and final number of events.
  if (finalMultiplicity < fMultiplicityMin) {return;}
  if (!fWriteMinimum) {fHistoMultiplicity[fCentralityBin][1]->Fill(finalMultiplicity);}
  fHistoNumberEvents[fCentralityBin]->Fill(3.5);

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
    if (!isTrackSelected[iiTrack]) {continue;}

  // 7.2 Get a pointer to the selected track.
    AliAODTrack *aaTrack = dynamic_cast<AliAODTrack*>(fRecoEvent->GetTrack(iiTrack));
    if (!aaTrack) {continue;}
    if (!aaTrack->TestFilterBit(fFstFilter)) {continue;}

  // 7.3 Get all the needed variables.
    iPt[iIndex] = aaTrack->Pt();
    iEta[iIndex] = aaTrack->Eta();
    iPhi[iIndex] = aaTrack->Phi();
    iParticleWeights[iIndex] = 1.;
    if (fUseJEfficiency)
    {
      iEffCorr = fEfficiency->GetCorrection(iPt[iIndex], fFilterbitIndex, fFstCentrality);
      iEffInverse = 1.0/iEffCorr;
      iParticleWeights[iIndex] = iEffInverse;
      fHistoEfficiency->Fill(iEffCorr);
      fHistoEffInverse->Fill(iEffInverse);
    } // End: if (fUseJEfficiency).

  // 7.4 Increase the index in the observables' arrays.
    iIndex++;
  } // End: for (Int_t iiTrack = 0; iiTrack < multiplicity; iiTrack++).

// 8. Get the particle weights if chosen.
  if (fUseParticleWeights)
  {
    Int_t eventRunNumber = fRecoEvent->GetRunNumber();
    CalculateWeight(eventRunNumber, finalMultiplicity, iParticleWeights, iPhi, iPt, iEta);
  } // End: if (fUseParticleWeights).

// 9. Compute all the possible multi-particle correlations.
  CalculateQvectors(finalMultiplicity, iPhi, iParticleWeights);
  ComputeReducedQvectors(finalMultiplicity);
  if (!fComputeACs) // Compute the correlators needed for SCs.
  ComputeSCsCorrelators(finalMultiplicity, iPhi, iParticleWeights);
  if (fComputeACs)  // Compute the correlators needed for ACs.
  {ComputeACsCorrelators(finalMultiplicity, iPhi, iParticleWeights);}
  if (fComputeEtaGaps)
  {ComputeTPCWithEtaGaps(finalMultiplicity, iPhi, iParticleWeights, iEta);}

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
/* Do the event/track selection at kine level to get the distributions....................... /
/ 1. Check the presence of the pointer to the reco and kine events.                           /
/ 2. Keep only the reco events belonging to the current centrality bin (The                   /
/   centrality estimators are based on detectors' informations, so it does not                /
/   work for kine events).                                                                    /
/ 3. Apply the event selection at kine level.                                                 /
/ 4. Get the mapping table between reco and kine levels.                                      /
/ 5. Keep the events with enough tracks to have a meaningful event weight.                    /
/ 6. Get the distributions at kine level.                                                     /
/ 7. Reset the variables for the next event.                                                 */
  TString sMethod = "void Ali-TwoMultiCorrelations::GetRatioDistributions()";

// 1. Check the presence of the pointer to the reco and kine events.
  if (!fRecoEvent) {Fatal(sMethod.Data(), "ERROR: no AOD event found.");}
  if (!fKineEvent) {Fatal(sMethod.Data(), "ERROR: no MC event found.");}

// 2. Keep only the reco events belonging to the current centrality bin.
//    Save the number of events in this centrality bin.
  if (!PassCentralitySelection()) {return;}
  fHistoNumberEvents[fCentralityBin]->Fill(0.5);

// 3.1 Apply the event selection at reco level if the table is used.
  if (fUseKineRecoTable) { if(!ApplyEventSelection(kTRUE)) {return;} }

// 3.2 Apply the event selection at kine level.
//    Save the number of events after the event selection.
  if(!ApplyEventSelection(kFALSE)) {return;}
  fHistoNumberEvents[fCentralityBin]->Fill(1.5);

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
    if (!aTrack) {continue;}
    if (!aTrack->TestFilterBit(fFstFilter)) {continue;}

  // 4.2.2 Apply the track selection at reco level.
    if (ApplyTrackSelection(aTrack))  // The track passed the selection: MC-key tagged with "1".
    {
      if (aTrack->GetLabel() >= 0)
      {
        if (kineRecoMap->GetValue(aTrack->GetLabel()) == 0) {kineRecoMap->Add(aTrack->GetLabel(), 1);}
      }
      finalMultiplicity++;
    } // End: the track passed the selection.
    else  // The track failed the selection: MC-key tagged with "-1".
    {
      if (aTrack->GetLabel() >= 0)
      {
        if (kineRecoMap->GetValue(aTrack->GetLabel()) == 0) {kineRecoMap->Add(aTrack->GetLabel(), -1);}
      }
    } // End: the track failed the selection.
  } // End: for (long long iTrack = 0; iTrack < multiplicity; iTrack++).

// 5. Keep the events with enough tracks to have a meaningful event weight.
//    Save the final number of tracks and final number of events.
  if (finalMultiplicity < fMultiplicityMin) {return;}
  fHistoMultiplicity[fCentralityBin][1]->Fill(finalMultiplicity);
  fHistoNumberEvents[fCentralityBin]->Fill(3.5);

// 6. Get the distributions at kine level.
  multiplicity = 0;
  multiplicity = fKineEvent->GetNumberOfTracks();

  for (long long iTrack = 0; iTrack < multiplicity; iTrack++)
  {
  // 6.1 Get the pointer to the current track.
    AliAODMCParticle *aMCtrack = dynamic_cast<AliAODMCParticle*>(fKineEvent->GetTrack(iTrack));
    if (!aMCtrack) {continue;}

	  Float_t phi = aMCtrack->Phi();  // Azimuthal angle.
    Float_t eta = aMCtrack->Eta();  // Pseudorapidity.
    Float_t pt = aMCtrack->Pt();    // Tranverse momentum.
    Int_t   charge = aMCtrack->Charge();    // Charge.

  // 6.2 Apply the track selection at kine level.
    if (!ApplyTrackSelection(aMCtrack)) {continue;}
    if (fUseKineRecoTable)
	  {
		  if (kineRecoMap->GetValue(iTrack)==-1) {continue;}
        // GetValue = -1: the corresponding reco track has not been selected by its track selection.
		  else
		  {
      // GetValue = 0: the kine track has been lost in ALICE --> Included in the pT distribution.
      // GetValue = 1: the reco track has been selected.
        fHistoPt[fCentralityBin][1]->Fill(pt);
        fHistoEta[fCentralityBin][1]->Fill(eta);
        fHistoPhi[fCentralityBin][1]->Fill(phi);
        fHistoCharge[fCentralityBin][1]->Fill(charge);
		  }
	  }
	  else
	  {
        fHistoPt[fCentralityBin][1]->Fill(pt);
        fHistoEta[fCentralityBin][1]->Fill(eta);
        fHistoPhi[fCentralityBin][1]->Fill(phi);
        fHistoCharge[fCentralityBin][1]->Fill(charge);
	  } // End: if (fUseKineRecoTable).
  } // End: for (long long iTrack = 0; iTrack < multiplicity; iTrack++).

// 7. Reset the variables for the next event.
  fInitialMultiplicity = 0;
  multiplicity = 0;
  finalMultiplicity = 0;
  kineRecoMap->Delete();  // Empty the mapping table.
  delete kineRecoMap; // Delete the pointer to the mapping table.

} // End: void GetRatioDistributions().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::PassCentralitySelection()
{
/* Check if the reconstructed event belongs to the current centrality bin.................... /
/ 1. Check the declaration of the two centrality estimators.                                  /
/ 2. Get the centrality for the two estimators.                                               /
/ 3. Fill the initial centrality distributions.                                               /
/ 4. Apply the centrality selection.                                                          /
/ 5. Determine the index of the centrality bin for the current event.                         /
/ 6. Fill the final centrality distributions.                                                */
  TString sMethod = "Bool_t Ali-TwoMultiCorrelations::PassCentralitySelection()";

// 1. Check the declaration of the two centrality estimators.
  if (fFstCentralityEstim == "")
  {Fatal(sMethod.Data(), "FATAL: no main centrality estimator has been set!");}
  if (fSndCentralityEstim == "")
  {Fatal(sMethod.Data(), "FATAL: no second centrality estimator has been set!");}

// 2. Get the centrality for the two estimators.
  AliMultSelection *multSelection = (AliMultSelection*)fRecoEvent->FindListObject("MultSelection");
  if (!multSelection) {Fatal(sMethod.Data(), "FATAL: no multSelection found!");}

  fFstCentrality = multSelection->GetMultiplicityPercentile(Form("%s", fFstCentralityEstim.Data()));
  fSndCentrality = multSelection->GetMultiplicityPercentile(Form("%s", fSndCentralityEstim.Data()));

// 3. Fill the initial centrality distributions.
  if (!fWriteMinimum)
  {
    fHistoInitCentrality[0]->Fill(fFstCentrality);
    fHistoInitCentrality[1]->Fill(fSndCentrality);
    if (fGetEstimCorrel)  // Fill the correlations between the estimators.
    {fHistoInitCorrelEstim->Fill(fFstCentrality, fSndCentrality);}
  } // End: if (!fWriteMinimum).

// 4. Reject the event if not in the analysis range.
  if ( (fFstCentrality < fCentralityMin) || (fFstCentrality >= fCentralityMax) )
  {return kFALSE;}
  // TODO: Insert here the cuts for the 2d centrality correlations?

// 5. Determine the index of the centrality bin for the current event.
//    The analysis always uses the same division of centralities.
  if ( (fFstCentrality >= 0.) && (fFstCentrality < 5.) ) {fCentralityBin = 0;}
  if ( (fFstCentrality >= 5.) && (fFstCentrality < 10.) ) {fCentralityBin = 1;}
  if ( (fFstCentrality >= 10.) && (fFstCentrality < 20.) ) {fCentralityBin = 2;}
  if ( (fFstCentrality >= 20.) && (fFstCentrality < 30.) ) {fCentralityBin = 3;}
  if ( (fFstCentrality >= 30.) && (fFstCentrality < 40.) ) {fCentralityBin = 4;}
  if ( (fFstCentrality >= 40.) && (fFstCentrality < 50.) ) {fCentralityBin = 5;}
  if ( (fFstCentrality >= 50.) && (fFstCentrality < 60.) ) {fCentralityBin = 6;}
  if ( (fFstCentrality >= 60.) && (fFstCentrality < 70.) ) {fCentralityBin = 7;}
  if ( (fFstCentrality >= 70.) && (fFstCentrality < 80.) ) {fCentralityBin = 8;}

// 6. Fill the final centrality distributions.
  if (!fWriteMinimum)
  {
    fHistoFstCentrality[fCentralityBin]->Fill(fFstCentrality);
    fHistoSndCentrality[fCentralityBin]->Fill(fSndCentrality);
    if (fGetEstimCorrel)
    {fHistoFinCorrelEstim[fCentralityBin]->Fill(fFstCentrality, fSndCentrality);}
  } // End: if (!fWriteMinimum).

  return kTRUE;
} // End: Bool_t PassCentralitySelection().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyEventSelection(Bool_t isRecoEvent)
{
/* Apply the event selection to the current event............................................ /
/ 1. Get the right variables depending on the type of event.                                  /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the physics event selection.                                                       /
/ 4. Fill the final QA histograms.                                                            /
/ 5. Reset the variables for the next event.                                                 */
  TString sMethod = "void Ali-TwoMultiCorrelations::ApplyEventSelection(Bool_t)";
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
  if (!fWriteMinimum)
  {
    if (!isRecoEvent)
    {
      fHistoMultiplicity[fCentralityBin][0]->Reset();
      fHistoPVx[fCentralityBin][0]->Reset();
      fHistoPVy[fCentralityBin][0]->Reset();
      fHistoPVz[fCentralityBin][0]->Reset();
    } // End: if (!isRecoEvent).

    fHistoMultiplicity[fCentralityBin][0]->Fill(fInitialMultiplicity);
    fHistoPVx[fCentralityBin][0]->Fill(PVx);
    fHistoPVy[fCentralityBin][0]->Fill(PVy);
    fHistoPVz[fCentralityBin][0]->Fill(PVz);
  } // End: if (!fWriteMinimum).

// 3. Apply the physics criteria of the event selection.
  if (fCutOnPVx)
  { if ( (PVx < fPVxMin) || (PVx > fPVxMax) ) {return kFALSE;} }
  if (fCutOnPVy)
  { if ( (PVy < fPVyMin) || (PVy > fPVyMax) ) {return kFALSE;} }
  if (fCutOnPVz)
  { if ( (PVz < fPVzMin) || (PVz > fPVzMax) ) {return kFALSE;} }

// 4. Fill the final QA histograms.
//    Reset the histogram for the filling at kine level.
  if (!fWriteMinimum)
  {
    if (!isRecoEvent)
    {
      fHistoPVx[fCentralityBin][1]->Reset();
      fHistoPVy[fCentralityBin][1]->Reset();
      fHistoPVz[fCentralityBin][1]->Reset();
    } // End: if (!isRecoEvent).

    fHistoPVx[fCentralityBin][1]->Fill(PVx);
    fHistoPVy[fCentralityBin][1]->Fill(PVy);
    fHistoPVz[fCentralityBin][1]->Fill(PVz);
  } // End: if (!fWriteMinimum).

// 5. Reset the variables for the next event.
  PVx = 0.;
  PVy = 0.;
  PVz = 0.;

  return kTRUE;
} // End: Bool_t ApplyEventSelection(Bool_t).

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::RemoveHMOs()
{
/* Get and apply the criteria to remove the HMOs in AOD events............................... /
/ 1. Get the number of tracks after the physics event selection.                              /
/ 2.1 Get a pointer to a reconstructed track.                                                 /
/ 2.2 Get the number of tracks in the main and global filters.                                /
/ 3. Fill the initial QA histograms.                                                          /
/ 4. Apply the HMOs cuts.                                                                     /
/ 5. Fill the final QA histograms.                                                            /
/ 6. Reset the variables for the next event.                                                 */
  TString sMethod = "void Ali-TwoMultiCorrelations::RemoveHMOs()";
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
    if (!aTrack) {continue;}

  // 2.2 Get the number of tracks in the main and global filters.
    if (aTrack->TestFilterBit(fSndFilter)) { globalMultiplicity++; }
    if (aTrack->TestFilterBit(fFstFilter)) { mainMultiplicity++; }
  } // End: iTrack.

// 3. Fill the initial QA histograms.
  if (!fWriteMinimum && !fDoKineAnalysis)
  {
    fHistoMultiplicityGlobal[fCentralityBin][0]->Fill(globalMultiplicity);
    fHistoMultiplicityMain[fCentralityBin][0]->Fill(mainMultiplicity);
    if (fGetFiltersCorrel)
    {fHistoCorrelFilters[0]->Fill(globalMultiplicity, mainMultiplicity);}
  } // End: if (!fWriteMinimum && !fDoKineAnalysis).

// 4. Apply the HMOs cuts.
  if (fCutOnHMOs)
  {
  // Minimum and maximum lines which delimite the accepted band.
    Float_t lineMin = fMultiplicityMinA*(Float_t)globalMultiplicity + fMultiplicityMinB;
    Float_t lineMax = fMultiplicityMaxA*(Float_t)globalMultiplicity + fMultiplicityMaxB;

  // The number of tracks in the main filter is under the accepted band.
    if ( (Float_t)mainMultiplicity < lineMin ) {return kFALSE;}
  // The number of tracks in the main filter is above the accepted band.
    if ( (Float_t)mainMultiplicity > lineMax ) {return kFALSE;}
  } // End: if (fCutOnHMOs).

// 5. Fill the final QA histograms.
  if (!fWriteMinimum && !fDoKineAnalysis)
  {
    fHistoMultiplicityGlobal[fCentralityBin][1]->Fill(globalMultiplicity);
    fHistoMultiplicityMain[fCentralityBin][1]->Fill(mainMultiplicity);
    if (fGetFiltersCorrel)
    {fHistoCorrelFilters[1]->Fill(globalMultiplicity, mainMultiplicity);}
  } // End: if (!fWriteMinimum && !fDoKineAnalysis).

// 6. Reset the variables for the next event.
  multiplicity = 0;
  globalMultiplicity = 0;
  mainMultiplicity = 0;

  return kTRUE;
} // End: Bool_t RemoveHMOs().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)
{
/* Apply the track selection criteria........................................................ /
/ 1. Declare all the observables for the reco track.                                          /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the track selection.                                                               /
/ 4. Fill the final QA histograms.                                                            /
/ 5. Reset the variables for the next track.                                                 */
  TString sMethod = "Bool_t Ali-TwoMultiCorrelations::ApplyTrackSelection(AliAODTrack *aAODtrack)";

// 1. Declare all the observables for the reco track.
  Float_t pT = aAODtrack->Pt();           // Transverse momentum.
  Float_t eta = aAODtrack->Eta();         // Pseudorapidity.
  Float_t phi = aAODtrack->Phi();         // Azimuthal angle.
  Int_t NTPC = aAODtrack->GetTPCNcls();   // Number of TPC clusters (default method).
  Int_t NITS = aAODtrack->GetITSNcls();   // Number of ITS clusters.
  Int_t charge = aAODtrack->Charge();     // Electric charge.

// Get the DCA information (cf PWGCF/EBYE/BalanceFunctions/AliAnalysisTaskBFPsi.cxx)
  Float_t DCAxy = 999.;   // DCA in the xy-plane.
  Float_t DCAz = 999.;    // DCA along z.

  if (fFstFilter == 128)  // These methods work only for constrained TPConly tracks.
  { // These two quantities are the DCA from global tracks but not what we will cut on.
    DCAxy = aAODtrack->DCA();
    DCAz = aAODtrack->ZAtDCA();
  }
  else  // For the unconstrained tracks.
  {
    AliAODVertex *primaryVertex = (AliAODVertex*)fRecoEvent->GetPrimaryVertex();
    Double_t v[3];    // Coordinates of the PV?
    Double_t pos[3];  // Coordinates of the track closest to PV?

    primaryVertex->GetXYZ(v);
    aAODtrack->GetXYZ(pos);
    DCAxy = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
    DCAz = pos[2] - v[2];
  }

// Get the chi^2 per TPC cluster.
  Float_t chiSquare = 999.;
  /// Personal method, should be equal to GetTPCchi2perCluster().
  if (fPersoChiSquare == 1) {chiSquare = (aAODtrack->GetTPCchi2())/(aAODtrack->GetNcls(1));}
  else if (fPersoChiSquare == 2) {chiSquare = aAODtrack->GetTPCchi2perCluster();}
  else if (fPersoChiSquare == 3) {chiSquare = aAODtrack->GetTPCchi2perNDF();}
  else {chiSquare = aAODtrack->Chi2perNDF();}

// 2. Fill the initial QA histograms.
//    Fill them only if it is full analysis with full writing mode.
  if (!fDoKineAnalysis && !fWriteMinimum)
  {
    fHistoPt[fCentralityBin][0]->Fill(pT);
    fHistoEta[fCentralityBin][0]->Fill(eta);
    fHistoPhi[fCentralityBin][0]->Fill(phi);
    fHistoNTPC[fCentralityBin][0]->Fill(NTPC);
    fHistoChiSquare[fCentralityBin][0]->Fill(chiSquare);
    fHistoNITS[fCentralityBin][0]->Fill(NITS);
    fHistoDCAxy[fCentralityBin][0]->Fill(DCAxy);
    fHistoDCAz[fCentralityBin][0]->Fill(DCAz);
    fHistoCharge[fCentralityBin][0]->Fill(charge);
  } // End: if (!fDoKineAnalysis && !fWriteMinimum).

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
    if (fFstFilter == 128) {Fatal(sMethod.Data(), "ERROR: TPConly tracks do not have ITS information.");}
    else if (NITS < fNITSMin) {return kFALSE;}
  }
  if (fCutOnDCAxy)
  { if (TMath::Abs(DCAxy) > fDCAxyMax) {return kFALSE;} }
  if (fCutOnDCAz)
  { if (TMath::Abs(DCAz) > fDCAzMax) {return kFALSE;} }
  if (fCutOnCharge)
  {
  // Keep only the positive tracks.
    if (fKeepPosCharges) { if (charge < 0) {return kFALSE;} }
  // Keep only the negative tracks.
    else if (!fKeepPosCharges) { if (charge > 0) {return kFALSE;} }
    else { Fatal(sMethod.Data(), "ERROR: Select the sign of the charges to keep."); }
  }

// 4. Fill the final QA histograms.
//    Fill them only if it is full analysis with full writing.
  if (!fDoKineAnalysis && !fWriteMinimum)
  {
    fHistoPt[fCentralityBin][1]->Fill(pT);
    fHistoEta[fCentralityBin][1]->Fill(eta);
    fHistoPhi[fCentralityBin][1]->Fill(phi);
    fHistoNTPC[fCentralityBin][1]->Fill(NTPC);
    fHistoChiSquare[fCentralityBin][1]->Fill(chiSquare);
    fHistoNITS[fCentralityBin][1]->Fill(NITS);
    fHistoDCAxy[fCentralityBin][1]->Fill(DCAxy);
    fHistoDCAz[fCentralityBin][1]->Fill(DCAz);
    fHistoCharge[fCentralityBin][1]->Fill(charge);
  } // End: if (!fDoKineAnalysis && !fWriteMinimum).

// 5. Reset the variables for the next track.
  pT = 0.;
  eta = 0.;
  phi = 0.;
  NTPC = 0;
  chiSquare = 0.;
  NITS = 0;
  //DCAx = 0.;
  //DCAy = 0.;
  DCAz = 0.;
  charge = 0;
  DCAxy = 0.;

  return kTRUE;
} // End: Bool_t ApplyTrackSelection().

/* ----------------------------------------------------------------------------------------- */
Bool_t AliAnalysisTaskTwoMultiCorrelations::ApplyTrackSelection(AliAODMCParticle *aMCtrack)
{
/* Apply the track selection criteria to the MC particle..................................... /
/ 1. Declare all the observables for the kine track.                                          /
/ 2. Fill the initial QA histograms.                                                          /
/ 3. Apply the track selection.                                                               /
/ 4. Reset the variables for the next track.                                                 */
  TString sMethodName = "Bool_t Ali-TwoMultiCorrelations::ApplyTrackSelection(AliAODMCParticle *aMCtrack)";

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
  fHistoPt[fCentralityBin][0]->Fill(pT);
  fHistoEta[fCentralityBin][0]->Fill(eta);
  fHistoPhi[fCentralityBin][0]->Fill(phi);
  fHistoCharge[fCentralityBin][0]->Fill(charge);

// 3. Apply the track selection.
  if (fKeepWeakSecondaries)   // Both the primaries and weak secondaries can be selected.
  { if (!isPrimary && !isWeakSecondary) {return kFALSE;} }
  else { if (!isPrimary) {return kFALSE;} } // Only the primaries can be selected.

  if (charge == 0) {return kFALSE;}
  if (fCutOnPt)
  {if ( (pT < fPtMin) || (pT > fPtMax) ) {return kFALSE;}}
  if (fCutOnEta)
  {if ( (eta < fEtaMin) || (eta > fEtaMax) ) {return kFALSE;}}
  if (fCutOnCharge)
  {
    if (fKeepPosCharges) {if (charge < 0) {return kFALSE;}} // Keep only the positive tracks.
    else if (!fKeepPosCharges) {if (charge > 0) {return kFALSE;}} // Keep only the negative tracks.
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
void AliAnalysisTaskTwoMultiCorrelations::CalculateWeight(Int_t runNumber, long long numberOfParticles, Float_t* pWeights, Float_t* angles, Float_t* pt, Float_t* eta)
{
/* Get the particle weights for all the tracks in the current event.......................... /
/ 1. Get the index corresponding to the current run.                                          /
/ 2. Get the phi, pT and eta weights.                                                         /
/ 3. Fill "pWeights".                                                                        */

// 1. Get the index corresponding to the current run.
  Int_t indexRun = GetRunIndex(runNumber);

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
      iBin = fHistoPhiWeight[indexRun][fCentralityBin]->FindBin(angles[i]);
      weight_phi = fHistoPhiWeight[indexRun][fCentralityBin]->GetBinContent(iBin);
    }
    if (fUsePtWeights)
    {
      iBin = fHistoPtWeight[indexRun][fCentralityBin]->FindBin(pt[i]);
      weight_pt = fHistoPtWeight[indexRun][fCentralityBin]->GetBinContent(iBin);
    }
    if (fUseEtaWeights)
    {
      iBin = fHistoEtaWeight[indexRun][fCentralityBin]->FindBin(eta[i]);
      weight_eta = fHistoEtaWeight[indexRun][fCentralityBin]->GetBinContent(eta[i]);
    }

  // 3. Fill "pWeights".
    pWeights[i] = weight_phi*weight_pt*weight_eta;
  } // End: for (long long i = 0; i < numberOfParticles; i++).

} // End: void CalculateWeight(Int_t, long long, Float_t*, Float_t*, Float_t*, Float_t*).

/* ----------------------------------------------------------------------------------------- */
Int_t AliAnalysisTaskTwoMultiCorrelations::GetRunIndex(Int_t runNumber)
{
/* Return for the given run the index in the run-by-run arrays.                              */
  TString sMethod = "Int_t Ali-TwoMultiCorrelations::GetRunIndex()";
  Int_t cRun = -1; // Current index in the loop.

// Find the position of the given run into the list of runs.
  for (Int_t iRun = 0; iRun < fNumberRuns; iRun++)
  {
    if (fListRuns[iRun] == runNumber)
    {
      cRun = iRun;
      break;
    } // End: for (Int_t iRun = 0; iRun < fNumberRuns; iRun++).
  } // End: iRun.

  return cRun;
} // End: Int_t GetRunIndex(Int_t).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Calculate all the Q-vectors for the given arrays of azimuthal angles and particle weights. */
  Float_t   iAngle            = 0.;   // Azimuthal angle of the current track.
  Float_t   iWeight           = 0.;   // Particle weight of the current track.
  Float_t   iWeightToPowerP   = 0.;   // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 11; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (long long iTrack = 0; iTrack < numberOfParticles; iTrack++)
  {
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < 11; iPower++)
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
    if (!fWriteMinimum) {fHistoReducedQvectors[fCentralityBin][i]->Fill(reducedQ);}
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
void AliAnalysisTaskTwoMultiCorrelations::ComputeSCsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Compute and save all the multiparticle correlators for the SCs for the current event. */
  Int_t twoZerosArray[2] = {0};     // 2d array with zeros for the denominator.
  Int_t twoPartArray[2] = {0};      // 2d array for the current harmonic.
  Int_t fourZerosArray[4] = {0};    // 4d array with zeros for the denominator.
  Int_t fourPartArray[4] = {0};     // 4d array for the current couple of harmonics.
  Int_t sixZerosArray[6] = {0};     // 6d array with zeros for the denominator.
  Int_t sixPartArray[6] = {0};      // 6d array for the current triplet of harmonics.
  Double_t eventWeight = 0.;        // Event weight and denominator of the correlator.
  TComplex complexCorrelator = TComplex(0., 0.);   // Complex form of the correlator.
  Double_t realCorrelator = 0.;     // Real part of the correlator.
  Int_t iBin = 1;                   // Bin for the current correlator.
  Double_t iMiddleBin = 0.;         // Middle of the bin for the current correlator.

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
    fProfileTwoPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
    fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(n, Form("%d", n));

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
      fProfileFourPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
      fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

    // Reset of the variables for the next pair of harmonics.
      iBin++;
      complexCorrelator = TComplex(0., 0.);
      realCorrelator    = 0.;
      iMiddleBin        = 0.;
    }   // End of the loop over the second harmonic for the 4-p correlators.
  }   // End of the loop over the first harmonic for the 4-p correlators.

  iBin                  = 1;
  if (!fWriteMinimum) // Compute the cross-check histograms only in full writing mode.
  {
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
        fProfileFourPartCorrelCheck[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
        fProfileFourPartCorrelCheck[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

      // Reset of the variables for the next pair of harmonics.
       iBin++;
       complexCorrelator  = TComplex(0., 0.);
       realCorrelator     = 0.;
       iMiddleBin         = 0.;
      }   // End of the loop over the second harmonic for cross-check.
    }   // End of the loop over the first harmonic for cross-check.
  } // End: if (!fWriteMinimum).
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
        fProfileSixPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
        fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d,%d)", l, m, n));

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

} // End of "void ComputeSCsCorrelators()".

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::ComputeACsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Compute and save all the multiparticle correlators for the ACs for the current event. */

// Compute the 2-particle correlations.
  Int_t twoZerosArray[2] = {0}; // Null 2d array for the event weight and denominator.
  Double_t eventWeight = (CalculateRecursion(2, twoZerosArray)).Re();
    // Event weight and denominator of the correlator.

/// <<2>>_{k,-k}, first bin.
  Int_t twoPartArray[2] = {fACHarmoOne, -1*fACHarmoOne};  // Array of harmonics.
  TComplex complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
    // Complex form of the multiparticle correlator.
  Double_t realCorrelator = complexCorrelator.Re();  // Real part of the multiparticle correlator.
  fProfileTwoPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("%d", twoPartArray[0]));

  realCorrelator = 0.;  // Partial reset for the next 2-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<2>>_{l,-l}, second bin.
  twoPartArray[0] = fACHarmoTwo; twoPartArray[1] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTwoPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("%d", twoPartArray[0]));

  realCorrelator = 0.;  // Partial reset for the next 2-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<2>>_{m,-m}, third bin.
  twoPartArray[0] = fACHarmoThree; twoPartArray[1] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTwoPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("%d", twoPartArray[0]));

// Full reset for the 4-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 4-particle correlations.
  Int_t fourZerosArray[4] = {0};  // Null 4d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(4, fourZerosArray)).Re();

/// <<4>>_{k,k,-k,-k}., first bin.
  Int_t fourPartArray[4] = {fACHarmoOne, fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{k,l,-k,-l}, second bin.
  fourPartArray[1] = fACHarmoTwo; fourPartArray[3] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{k,m,-k,-m}, third bin.
  fourPartArray[1] = fACHarmoThree; fourPartArray[3] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{l,m,-l,-m}, fourth bin.
  fourPartArray[0] = fACHarmoTwo; fourPartArray[2] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(3.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(4, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

// Full reset for the 6-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 6-particle correlations.
  Int_t sixZerosArray[6] = {0};   // Null 6d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(6, sixZerosArray)).Re();

/// <<6>>_{k,k,k,-k,-k,-k}, first bin.
  Int_t sixPartArray[6] = {fACHarmoOne, fACHarmoOne, fACHarmoOne,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,k,l,-k,-k,-l}, second bin.
  sixPartArray[2] = fACHarmoTwo; sixPartArray[5] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,k,m,-k,-k,-m}, third bin.
  sixPartArray[2] = fACHarmoThree; sixPartArray[5] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,l,m,-k,-l,-m}, fourth bin.
  sixPartArray[1] = fACHarmoTwo; sixPartArray[4] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(3.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(4, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

// Full reset for the 8-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 8-particle correlations.
  Int_t eightZerosArray[8] = {0};   // Null 8d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(8, eightZerosArray)).Re();

/// <<8>>_{k,k,k,k,-k,-k,-k,-k}, first bin.
  Int_t eightPartArray[8] = {fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoOne,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

  realCorrelator = 0.;  // Partial reset for the next 8-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<8>>_{k,k,k,l,-k,-k,-k,-l}, second bin.
  eightPartArray[3] = fACHarmoTwo; eightPartArray[7] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

  realCorrelator = 0.;  // Partial reset for the next 8-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<8>>_{k,k,l,m,-k,-k,-l,-m}, third bin.
  eightPartArray[2] = fACHarmoTwo; eightPartArray[6] = -1*fACHarmoTwo;
  eightPartArray[3] = fACHarmoThree; eightPartArray[7] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

/// Full reset for the 10-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 10-particle correlations.
  Int_t tenZerosArray[10] = {0};  // Null 10d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(10, tenZerosArray)).Re();

/// <<10>>_{k,k,k,k,l,-k,-k,-k,-k,-l}, first bin.
  Int_t tenPartArray[10] = {fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoTwo,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoTwo};
  complexCorrelator = (CalculateRecursion(10, tenPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTenPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileTenPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d,%d,%d)", tenPartArray[0], tenPartArray[1], tenPartArray[2], tenPartArray[3], tenPartArray[4]));

// Reset all variables for the next event.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

} // End of "void ComputeACsCorrelators()".

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
    if (fUseParticleWeights) {iWeightToP = iWeight;}   // All weights are multiplied to get the final one.

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
      fProfileTPCEta[fCentralityBin][iGap]->Fill(iHarmo + 0.5, realCorrel, Mminus[iGap][iHarmo]*Mplus[iGap][iHarmo]);

    // Reset the 2-particle correlator.
      complexCorrel = TComplex(0.,0.);
      realCorrel    = 0.;
    }   // End of the loop over the gaps.
  }   // End of the loop over the harmonics.

}   // End of "void ComputeTPCWithEtaGaps()".

/* ========================================================================================== /
/ Methods called in "Terminate".                                                              /
/ ========================================================================================== */
/* ----------------------------------------------------------------------------------------- */
/* ========================================================================================== /
/ Methods called in "UserCreateOutputObjects".                                                /
/ ========================================================================================== */
void AliAnalysisTaskTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists in fMainList........................................................... /
// 1. Check if the mother list exists.                                                        /
// 2. Book the daughter list for the minimum QA histograms.                                   /
// 3. Book the daughter list for the multiplicity histograms.                                 /
// 4. Book the daughter list for the event QA histograms.                                     /
// 5. Book the daughter list for the track QA histograms.                                     /
// 6. Book the daughter list for the MPC histograms.                                          /
// 7. Book the daughter list for the eta gaps profiles.                                      */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookAllLists()";

// 1. Check if the mother list exists.
  if (!fMainList) {Fatal(sMethodName.Data(), "ERROR: 'fMainList' does not exist.");}

// 2. Book the daughter list for the minimum QA histograms.
  fMinimumQAList = new TList();
  fMinimumQAList->SetName("fMinimumQAList");
  fMinimumQAList->SetOwner(kTRUE);
  fMainList->Add(fMinimumQAList);

// 3. Book the daughter list for the multiplicity histograms
  if (!fWriteMinimum)
  {
    fMultiplicityList = new TList();
    fMultiplicityList->SetName("fMultiplicityList");
    fMultiplicityList->SetOwner(kTRUE);
    fMainList->Add(fMultiplicityList);
  } // End: if (!fWriteMinimum).

// 4. Book the daughter list for the event QA histograms.
  if (!fWriteMinimum)
  {
    fEventQAList = new TList();
    fEventQAList->SetName("fEventQAList");
    fEventQAList->SetOwner(kTRUE);
    fMainList->Add(fEventQAList);
  } // End: if (!fWriteMinimum).

// 5. Book the daughter list for the track QA histograms.
  if (!fWriteMinimum && !fGetFiltersCorrel)
  {
    fTrackQAList = new TList();
    fTrackQAList->SetName("fTrackQAList");
    fTrackQAList->SetOwner(kTRUE);
    fMainList->Add(fTrackQAList);
  } // End: if (!fWriteMinimum).

// 6. Book the daughter list for the MPC histograms.
  if (!fGetFiltersCorrel)
  {
    fMPCList = new TList();
    fMPCList->SetName("fMPCList");
    fMPCList->SetOwner(kTRUE);
    fMainList->Add(fMPCList);
  } // End: if (!fGetFiltersCorrel)

// 7. Book the daughter list for the eta gaps profiles.
  if (fComputeEtaGaps && !fGetFiltersCorrel)
  {
    fTPCEtaList = new TList();
    fTPCEtaList->SetName("fTPCEtaList");
    fTPCEtaList->SetOwner(kTRUE);
    fMainList->Add(fTPCEtaList);
  } // End: if (fComputeEtaGaps).

} // End: void BookAllLists().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMinimumQAList()
{
/* Book the QA histograms still needed in minimum writing mode............................... /
// Number of events at each selection step. 4th bin is used in the offline bootstrap.        */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookMinimumQAList()";

// Summary of the cuts in the event selection.
  fProfileEventCuts = new TProfile("", "", 15, 0., 15.);
  fProfileEventCuts->SetName("fProfileEventCuts");
  fProfileEventCuts->SetTitle("Configuration of the event selection");
  fProfileEventCuts->SetStats(kFALSE);
  fProfileEventCuts->GetXaxis()->SetBinLabel(1, "1st centrality");
  fProfileEventCuts->GetXaxis()->SetBinLabel(2, "2nd centrality");
  fProfileEventCuts->GetXaxis()->SetBinLabel(3, "PV_{x} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(4, "PV_{x} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(5, "PV_{y} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(6, "PV_{y} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(7, "PV_{z} min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(8, "PV_{z} max");
  fProfileEventCuts->GetXaxis()->SetBinLabel(9, "Multiplicity min");
  fProfileEventCuts->GetXaxis()->SetBinLabel(10, "1st filter");
  fProfileEventCuts->GetXaxis()->SetBinLabel(11, "2nd filter");
  fProfileEventCuts->GetXaxis()->SetBinLabel(12, "HMO min A");
  fProfileEventCuts->GetXaxis()->SetBinLabel(13, "HMO min B");
  fProfileEventCuts->GetXaxis()->SetBinLabel(14, "HMO max A");
  fProfileEventCuts->GetXaxis()->SetBinLabel(15, "HMO max B");
  fMinimumQAList->Add(fProfileEventCuts);

// Summary of the cuts in the track selection.
  fProfileTrackCuts = new TProfile("", "", 15, 0., 15.);
  fProfileTrackCuts->SetName("fProfileTrackCuts");
  fProfileTrackCuts->SetTitle("Configuration of the track selection");
  fProfileTrackCuts->SetStats(kFALSE);
  fProfileTrackCuts->GetXaxis()->SetBinLabel(1, "p_{T} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(2, "p_{T} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(3, "#eta min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(4, "#eta max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(5, "N_{TPC} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(6, "#chi^{2} method");  
  fProfileTrackCuts->GetXaxis()->SetBinLabel(7, "#chi^{2} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(8, "#chi^{2} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(9, "N_{ITS} min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(10, "DCA_{xy} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(11, "DCA_{z} max");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(12, "Charge");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(13, "p_{T} weights?");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(14, "#eta weights?");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(15, "#phi weights?");
  fMinimumQAList->Add(fProfileTrackCuts);

// Number of events at each main selection step.
  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.

    fHistoNumberEvents[iCent] = new TH1I("", "", 4, 0., 4.);
    fHistoNumberEvents[iCent]->SetName(Form("fHistoNumberEvents_Bin%d", iCent));
    fHistoNumberEvents[iCent]->SetTitle(Form("Number of events, Centrality bin %d", iCent));
    fHistoNumberEvents[iCent]->SetStats(kTRUE);
    fHistoNumberEvents[iCent]->GetXaxis()->SetBinLabel(1, "After centrality selection");
    fHistoNumberEvents[iCent]->GetXaxis()->SetBinLabel(2, "After event selection");
    fHistoNumberEvents[iCent]->GetXaxis()->SetBinLabel(3, "After track selection");
    fHistoNumberEvents[iCent]->GetXaxis()->SetBinLabel(4, "Final number of events");
    fHistoNumberEvents[iCent]->GetYaxis()->SetTitle("Number of events");
    fMinimumQAList->Add(fHistoNumberEvents[iCent]);
  } // End: iCent.
} // End. void BookMinimumQAList().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMultiplicityList()
{
/* Book the histograms related to the multiplicity distributions............................. /
// 1. Initial correlations between the estimators.                                            /
// 2. Initial centrality distributions.                                                       /
// 3. Multiplicity distributions.                                                             /
// 4. Multiplicity distributions per filter bits. Only meaningful for AOD events.             /
// 5. Correlations between filters. Only to decide HMOs!                                      /
// 6. Final centrality distributions.                                                         /
// 7. Final correlations between the estimators.                                             */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookMultiplicityList()";

// 1. Initial correlations between the estimators.
  if (fGetEstimCorrel)  // Printed only for checking purposes.
  {
    fHistoInitCorrelEstim = new TH2F("", "", 101, 0., 101., 101, 0., 101.);
    fHistoInitCorrelEstim->SetName("fHistoInitCorrelEstim");
    fHistoInitCorrelEstim->SetTitle(Form("Initial correlations between %s and %s",
        fFstCentralityEstim.Data(), fSndCentralityEstim.Data()));
    fHistoInitCorrelEstim->SetStats(kFALSE);
    fHistoInitCorrelEstim->GetXaxis()->SetTitle(Form("Centrality from %s",
        fFstCentralityEstim.Data()));
    fHistoInitCorrelEstim->GetYaxis()->SetTitle(Form("Centrality from %s",
        fSndCentralityEstim.Data()));
    fMultiplicityList->Add(fHistoInitCorrelEstim);
  } // End: if (fGetEstimCorrel).

  for (Int_t iEst = 0; iEst < 2; iEst++)
  {
  // 2. Initial centrality distributions.
    TString iEstim = fSndCentralityEstim;
    if (iEst == 0) {iEstim = fFstCentralityEstim;}

    fHistoInitCentrality[iEst] = new TH1F("", "", 101, 0., 101.);
    fHistoInitCentrality[iEst]->SetName(Form("fHistoInitCentrality%s", iEstim.Data()));
    fHistoInitCentrality[iEst]->SetTitle(Form("Initial distribution of the centrality for %s",
        iEstim.Data()));
    fHistoInitCentrality[iEst]->SetStats(kTRUE);
    fHistoInitCentrality[iEst]->GetXaxis()->SetTitle("Centrality percentile");
    fHistoInitCentrality[iEst]->GetYaxis()->SetTitle("Number of events");
    fMultiplicityList->Add(fHistoInitCentrality[iEst]);

    TString iState = "Final";
    if (iEst == 0) {iState = "Initial";}

  // 3. Correlations between filters.
      if (fGetFiltersCorrel)  // Only filled to get the HMOs cuts.
      {
        fHistoCorrelFilters[iEst] = new TH2I("", "", 5000, 0., 5000., 5000, 0., 5000.);
        fHistoCorrelFilters[iEst]->SetName(Form("fHistoCorrelFilters%s",
            iState.Data()));
        fHistoCorrelFilters[iEst]->SetTitle(Form("%s distribution of correlated filters",
            iState.Data()));
        fHistoCorrelFilters[iEst]->SetStats(kTRUE);
        fHistoCorrelFilters[iEst]->GetXaxis()->SetTitle(Form("Multiplicity_{Filter %d}",
            fSndFilter));
        fHistoCorrelFilters[iEst]->GetYaxis()->SetTitle(Form("Multiplicity_{Filter %d}",
            fFstFilter));
        fMultiplicityList->Add(fHistoCorrelFilters[iEst]);
      } // End: if (fGetFiltersCorrel).
  } // End: for (Int_t iEst = 0; iEst < 2; iEst++).

  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.

  // 4. Multiplicity distributions.
    for (Int_t iSt = 0; iSt < 2; iSt++)
    {
      TString iState = "Final";
      if (iSt == 0) {iState = "Initial";}

      fHistoMultiplicity[iCent][iSt] = new TH1I("", "",
          fNumberBinsMulti, 0., fNumberBinsMulti);
      fHistoMultiplicity[iCent][iSt]->SetName(Form("fHistoMultiplicity%s_Bin%d",
          iState.Data(), iCent));
      fHistoMultiplicity[iCent][iSt]->SetTitle(Form("%s distribution of the number of tracks, Centrality bin %d", iState.Data(), iCent));
      fHistoMultiplicity[iCent][iSt]->SetStats(kTRUE);
      fHistoMultiplicity[iCent][iSt]->GetXaxis()->SetTitle("Number of tracks");
      fHistoMultiplicity[iCent][iSt]->GetYaxis()->SetTitle("Number of events");
      fMultiplicityList->Add(fHistoMultiplicity[iCent][iSt]);

  // 5. Multiplicity distributions per filter bits.
      if (!fDoKineAnalysis) // Only filled for AOD events.
      {
      // Main filter bit.
        fHistoMultiplicityMain[iCent][iSt] = new TH1I("", "",
            fNumberBinsMulti, 0., fNumberBinsMulti);
        fHistoMultiplicityMain[iCent][iSt]->SetName(Form("fHistoMultiplicityMain%s_Bin%d",
            iState.Data(), iCent));
        fHistoMultiplicityMain[iCent][iSt]->SetTitle(Form("%s distribution of the number of tracks, Centrality bin %d (Filter %d)", iState.Data(), iCent, fFstFilter));
        fHistoMultiplicityMain[iCent][iSt]->SetStats(kTRUE);
        fHistoMultiplicityMain[iCent][iSt]->GetXaxis()->SetTitle("Number of tracks");
        fHistoMultiplicityMain[iCent][iSt]->GetYaxis()->SetTitle("Number of events");
        fMultiplicityList->Add(fHistoMultiplicityMain[iCent][iSt]);

      // Multiplicity distributions for the global filter bit.
        fHistoMultiplicityGlobal[iCent][iSt] = new TH1I("", "",
            fNumberBinsMulti, 0., fNumberBinsMulti);
        fHistoMultiplicityGlobal[iCent][iSt]->SetName(Form("fHistoMultiplicityGlobal%s_Bin%d",
            iState.Data(), iCent));
        fHistoMultiplicityGlobal[iCent][iSt]->SetTitle(Form("%s distribution of the number of tracks, Centrality bin %d (Filter %d)", iState.Data(), iCent, fSndFilter));
        fHistoMultiplicityGlobal[iCent][iSt]->SetStats(kTRUE);
        fHistoMultiplicityGlobal[iCent][iSt]->GetXaxis()->SetTitle("Number of tracks");
        fHistoMultiplicityGlobal[iCent][iSt]->GetYaxis()->SetTitle("Number of events");
        fMultiplicityList->Add(fHistoMultiplicityGlobal[iCent][iSt]);
      } // End: if (!fDoKineAnalysis).
    } // End: iSt.

  // 6. Final centrality distributions.
    fHistoFstCentrality[iCent] = new TH1F("", "", 101, 0., 101.);
    fHistoFstCentrality[iCent]->SetName(Form("fHistoFstCentrality_Bin%d", iCent));
    fHistoFstCentrality[iCent]->SetTitle(Form("Final distribution of the centrality for %s, Centrality bin %d", fFstCentralityEstim.Data(), iCent));
    fHistoFstCentrality[iCent]->SetStats(kTRUE);
    fHistoFstCentrality[iCent]->GetXaxis()->SetTitle("Centrality percentile");
    fHistoFstCentrality[iCent]->GetYaxis()->SetTitle("Number of events");
    fMultiplicityList->Add(fHistoFstCentrality[iCent]);

    fHistoSndCentrality[iCent] = new TH1F("", "", 101, 0., 101.);
    fHistoSndCentrality[iCent]->SetName(Form("fHistoSndCentrality_Bin%d", iCent));
    fHistoSndCentrality[iCent]->SetTitle(Form("Final distribution of the centrality for %s, Centrality bin %d", fSndCentralityEstim.Data(), iCent));
    fHistoSndCentrality[iCent]->SetStats(kTRUE);
    fHistoSndCentrality[iCent]->GetXaxis()->SetTitle("Centrality percentile");
    fHistoSndCentrality[iCent]->GetYaxis()->SetTitle("Number of events");
    fMultiplicityList->Add(fHistoSndCentrality[iCent]);

  // 7. Final correlations between estimators.
    if (fGetEstimCorrel)  // Printed only for checking purposes.
    {
      fHistoFinCorrelEstim[iCent] = new TH2F("", "", 101, 0., 101., 101, 0., 101.);
      fHistoFinCorrelEstim[iCent]->SetName(Form("fHistoFinCorrelEstim_Bin%d", iCent));
      fHistoFinCorrelEstim[iCent]->SetTitle(Form("Final correlations between %s and %s, Centrality bin %d", fFstCentralityEstim.Data(), fSndCentralityEstim.Data(), iCent));
      fHistoFinCorrelEstim[iCent]->SetStats(kFALSE);
      fHistoFinCorrelEstim[iCent]->GetXaxis()->SetTitle(Form("Centrality from %s",
          fFstCentralityEstim.Data()));
      fHistoFinCorrelEstim[iCent]->GetYaxis()->SetTitle(Form("Centrality from %s",
          fSndCentralityEstim.Data()));
      fMultiplicityList->Add(fHistoFinCorrelEstim[iCent]);
    } // End: if (fGetEstimCorrel).
  } // End: iCent.
} // End: void BookMultiplicityList().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookEventQAList()
{
/* Book the control histograms for the event selection criteria.............................. /
// 1. Primary Vertex x, y, z.                                                                */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookEventQAList()";

  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.

    for (Int_t iSt = 0; iSt < 2; iSt++)
    {
      TString iState = "Final";
      if (iSt == 0) {iState = "Initial";}

    // 1.1 Primary Vertex x.
      fHistoPVx[iCent][iSt] = new TH1F("", "", fNumberBinsPV, -1.*fMaxHistoPV, fMaxHistoPV);
      fHistoPVx[iCent][iSt]->SetName(Form("fHistoPVx%s_Bin%d", iState.Data(), iCent));
      fHistoPVx[iCent][iSt]->SetTitle(Form("%s distribution of PV_{x}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoPVx[iCent][iSt]->SetStats(kTRUE);
      fHistoPVx[iCent][iSt]->GetXaxis()->SetTitle("PV_{x} [cm]");
      fHistoPVx[iCent][iSt]->GetYaxis()->SetTitle("dN/dPV_{x}");
      fEventQAList->Add(fHistoPVx[iCent][iSt]);

    // 1.2 Primary Vertex y.
      fHistoPVy[iCent][iSt] = new TH1F("", "", fNumberBinsPV, -1.*fMaxHistoPV, fMaxHistoPV);
      fHistoPVy[iCent][iSt]->SetName(Form("fHistoPVy%s_Bin%d", iState.Data(), iCent));
      fHistoPVy[iCent][iSt]->SetTitle(Form("%s distribution of PV_{y}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoPVy[iCent][iSt]->SetStats(kTRUE);
      fHistoPVy[iCent][iSt]->GetXaxis()->SetTitle("PV_{y} [cm]");
      fHistoPVy[iCent][iSt]->GetYaxis()->SetTitle("dN/dPV_{y}");
      fEventQAList->Add(fHistoPVy[iCent][iSt]);

    // 1.3 Primary Vertex z.
      fHistoPVz[iCent][iSt] = new TH1F("", "", fNumberBinsPV, -1.*fMaxHistoPV, fMaxHistoPV);
      fHistoPVz[iCent][iSt]->SetName(Form("fHistoPVz%s_Bin%d", iState.Data(), iCent));
      fHistoPVz[iCent][iSt]->SetTitle(Form("%s distribution of PV_{z}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoPVz[iCent][iSt]->SetStats(kTRUE);
      fHistoPVz[iCent][iSt]->GetXaxis()->SetTitle("PV_{z} [cm]");
      fHistoPVz[iCent][iSt]->GetYaxis()->SetTitle("dN/dPV_{z}");
      fEventQAList->Add(fHistoPVz[iCent][iSt]);
    } // End: iSt.
  } // End: iCent.

} // End: void BookEventQAList().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookTrackQAList()
{
/* Book the control histograms for the track selection criteria.............................. /
// 1. Transverse momentum.                                                                    /
// 2. Pseudorapidity.                                                                         /
// 3. Azimuthal angles.                                                                       /
// 4. Number of TPC clusters.                                                                 /
// 5. chi^2/NDF of p per TPC clusters.                                                        /
// 6. Number of ITS clusters.                                                                 /
// 7. DCA in the xy-plane.                                                                    /
// 8. DCA along the z-axis.                                                                   /
// 9. Electric charge of the tracks.                                                         */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookTrackQAList()";

  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.

    for (Int_t iSt = 0; iSt < 2; iSt++)
    {
      TString iState = "Final";
      if (iSt == 0) {iState = "Initial";}

    // 1. Transverse momentum.
      fHistoPt[iCent][iSt] = new TH1F("", "", fNumberBinsPt, 0., 20.);
      fHistoPt[iCent][iSt]->SetName(Form("fHistoPt%s_Bin%d", iState.Data(), iCent));
      fHistoPt[iCent][iSt]->SetTitle(Form("%s distribution of p_{T}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoPt[iCent][iSt]->SetStats(kTRUE);
      fHistoPt[iCent][iSt]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistoPt[iCent][iSt]->GetYaxis()->SetTitle("dN/dp_{T}");
      fTrackQAList->Add(fHistoPt[iCent][iSt]);

    // 2. Pseudorapidity.
      fHistoEta[iCent][iSt] = new TH1F("", "", fNumberBinsEta, -1., 1.);
      fHistoEta[iCent][iSt]->SetName(Form("fHistoEta%s_Bin%d", iState.Data(), iCent));
      fHistoEta[iCent][iSt]->SetTitle(Form("%s distribution of #eta, Centrality bin %d",
          iState.Data(), iCent));
      fHistoEta[iCent][iSt]->SetStats(kTRUE);
      fHistoEta[iCent][iSt]->GetXaxis()->SetTitle("#eta");
      fHistoEta[iCent][iSt]->GetYaxis()->SetTitle("dN/d#eta");
      fTrackQAList->Add(fHistoEta[iCent][iSt]);

    // 3. Azimuthal angles.
      fHistoPhi[iCent][iSt] = new TH1F("", "", fNumberBinsPhi, 0., 6.3);
      fHistoPhi[iCent][iSt]->SetName(Form("fHistoPhi%s_Bin%d", iState.Data(), iCent));
      fHistoPhi[iCent][iSt]->SetTitle(Form("%s distribution of #varphi, Centrality bin %d",
          iState.Data(), iCent));
      fHistoPhi[iCent][iSt]->SetStats(kTRUE);
      fHistoPhi[iCent][iSt]->GetXaxis()->SetTitle("#varphi");
      fHistoPhi[iCent][iSt]->GetYaxis()->SetTitle("dN/d#varphi");
      fTrackQAList->Add(fHistoPhi[iCent][iSt]);

    // 4. Number of TPC clusters.
      fHistoNTPC[iCent][iSt] = new TH1I("", "", 170, 0., 170.);
      fHistoNTPC[iCent][iSt]->SetName(Form("fHistoNTPC%s_Bin%d", iState.Data(), iCent));
      fHistoNTPC[iCent][iSt]->SetTitle(Form("%s distribution of N_{TPC}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoNTPC[iCent][iSt]->SetStats(kTRUE);
      fHistoNTPC[iCent][iSt]->GetXaxis()->SetTitle("N_{TPC}");
      fHistoNTPC[iCent][iSt]->GetYaxis()->SetTitle("dN/dN_{TPC}");
      fTrackQAList->Add(fHistoNTPC[iCent][iSt]);

    // 5. chi^2/NDF of p per TPC clusters.
      fHistoChiSquare[iCent][iSt] = new TH1F("", "", 1000, 0., 20.);
      fHistoChiSquare[iCent][iSt]->SetName(Form("fHistoChiSquare%s_Bin%d",
          iState.Data(), iCent));
      fHistoChiSquare[iCent][iSt]->SetTitle(Form("%s distribution of #chi^{2}/NDF, Centrality bin %d", iState.Data(), iCent));
      fHistoChiSquare[iCent][iSt]->SetStats(kTRUE);
      fHistoChiSquare[iCent][iSt]->GetXaxis()->SetTitle("#chi^{2}/NDF of p in TPC");
      fHistoChiSquare[iCent][iSt]->GetYaxis()->SetTitle("Number of tracks");
      fTrackQAList->Add(fHistoChiSquare[iCent][iSt]);

    // 6. Number of ITS clusters.
      fHistoNITS[iCent][iSt] = new TH1I("", "", 14, 0., 7.);
      fHistoNITS[iCent][iSt]->SetName(Form("fHistoNITS%s_Bin%d", iState.Data(), iCent));
      fHistoNITS[iCent][iSt]->SetTitle(Form("%s distribution of N_{ITS}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoNITS[iCent][iSt]->SetStats(kTRUE);
      fHistoNITS[iCent][iSt]->GetXaxis()->SetTitle("Number of ITS clusters");
      fHistoNITS[iCent][iSt]->GetYaxis()->SetTitle("dN/dN_{ITS}");
      fTrackQAList->Add(fHistoNITS[iCent][iSt]);

    // 7. DCA in the xy-plane.
      fHistoDCAxy[iCent][iSt] = new TH1F("", "", 1000, -10., 10.);
      fHistoDCAxy[iCent][iSt]->SetName(Form("fHistoDCAxy%s_Bin%d", iState.Data(), iCent));
      fHistoDCAxy[iCent][iSt]->SetTitle(Form("%s distribution of DCA_{xy}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoDCAxy[iCent][iSt]->SetStats(kTRUE);
      fHistoDCAxy[iCent][iSt]->GetXaxis()->SetTitle("DCA_{xy} [cm]");
      fHistoDCAxy[iCent][iSt]->GetYaxis()->SetTitle("dN/dDCA_{xy}");
      fTrackQAList->Add(fHistoDCAxy[iCent][iSt]);

    // 8. DCA along the z-axis.
      fHistoDCAz[iCent][iSt] = new TH1F("", "", 1000, -10., 10.);
      fHistoDCAz[iCent][iSt]->SetName(Form("fHistoDCAz%s_Bin%d", iState.Data(), iCent));
      fHistoDCAz[iCent][iSt]->SetTitle(Form("%s distribution of DCA_{z}, Centrality bin %d",
          iState.Data(), iCent));
      fHistoDCAz[iCent][iSt]->SetStats(kTRUE);
      fHistoDCAz[iCent][iSt]->GetXaxis()->SetTitle("DCA_{z} [cm]");
      fHistoDCAz[iCent][iSt]->GetYaxis()->SetTitle("dN/dDCA_{z}");
      fTrackQAList->Add(fHistoDCAz[iCent][iSt]);

    // 9. Electric charge of the tracks.
      fHistoCharge[iCent][iSt] = new TH1I("", "", 2, -2, 2);
      fHistoCharge[iCent][iSt]->SetName(Form("fHistoCharge%s_Bin%d", iState.Data(), iCent));
      fHistoCharge[iCent][iSt]->SetTitle(Form("%s distribution of the charges, Centrality bin %d", iState.Data(), iCent));
      fHistoCharge[iCent][iSt]->SetStats(kTRUE);
      fHistoCharge[iCent][iSt]->GetXaxis()->SetTitle("Charge");
      fHistoCharge[iCent][iSt]->GetYaxis()->SetTitle("Number of tracks");
      fTrackQAList->Add(fHistoCharge[iCent][iSt]);
    } // End: iSt.
  } // End: iCent.

  fHistoEfficiency = new TH1F("fHistoEfficiency",
    "Distribution of the efficiency correction", 1000, 0., 20.);
  fHistoEfficiency->SetStats(kTRUE);
  fHistoEfficiency->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoEfficiency->GetYaxis()->SetTitle("dN/dp_{T}");
  fTrackQAList->Add(fHistoEfficiency);

  fHistoEffInverse = new TH1F("fHistoEffInverse",
    "Distribution of the inverse of the efficiency correction", 1000, 0., 20.);
  fHistoEffInverse->SetStats(kTRUE);
  fHistoEffInverse->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistoEffInverse->GetYaxis()->SetTitle("dN/dp_{T}");
  fTrackQAList->Add(fHistoEffInverse);

} // End: void BookTrackQAList().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookMPCList()
{
/* Book the histograms and profiles for the MPC............................................. /
// 1. */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookMPCList()";

// Define the numbers of bins depending on if the code computes SCs or ACs.
// The default binnings are for SCs.
  Int_t numberBinsTwo = 8;
  Int_t numberBinsFour = 36;
  Int_t numberBinsSix = 56;
  if (fComputeACs) {numberBinsTwo = 3; numberBinsFour = 4; numberBinsSix = 4;}

// Declare all the profiles with their settings.
  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.
    if (!fWriteMinimum) // Reduced q-vectors part of QA.
    {
    // Modulus of reduced Q-vectors.
      for (Int_t i = 0; i < 8; i++)
      {
        fHistoReducedQvectors[iCent][i] = new TH1F("", "", 1000, 0., 10.);
        fHistoReducedQvectors[iCent][i]->SetName(Form("fHistoReducedQvectors%d_Bin%d",
            i, iCent));
        fHistoReducedQvectors[iCent][i]->SetTitle(Form("Distribution of q_{%d}, Centrality bin %d", i+1, iCent));
        fHistoReducedQvectors[iCent][i]->SetStats(kTRUE);
        fHistoReducedQvectors[iCent][i]->GetXaxis()->SetTitle(Form("q_{%d}", i+1));
        fHistoReducedQvectors[iCent][i]->GetYaxis()->SetTitle(Form("dN/dq_{%d}", i+1));
        fMPCList->Add(fHistoReducedQvectors[iCent][i]);
      } // End: for (Int_t i = 0; i < 8; i++).
    } // End: if (!fWriteMinimum).

  // Profile with the 2-particle correlators.
    fProfileTwoPartCorrel[iCent] = new TProfile(Form("fProfileTwoPartCorrel_Bin%d", iCent),
      Form("#LT#LT2#GT#GT_{n,-n}, Centrality bin %d (no #eta gap)", iCent),
      numberBinsTwo, 0., numberBinsTwo);
    fProfileTwoPartCorrel[iCent]->Sumw2();
    fProfileTwoPartCorrel[iCent]->SetStats(kTRUE);
    fProfileTwoPartCorrel[iCent]->GetXaxis()->SetTitle("n");
    fMPCList->Add(fProfileTwoPartCorrel[iCent]);

  // Profiles with the 4-particle correlators.
    fProfileFourPartCorrel[iCent] = new TProfile(Form("fProfileFourPartCorrel_Bin%d", iCent),
      Form("#LT#LT4#GT#GT_{m,n,-m,-n}, Centrality bin %d", iCent), numberBinsFour, 0., numberBinsFour);
    fProfileFourPartCorrel[iCent]->Sumw2();
    fProfileFourPartCorrel[iCent]->SetStats(kTRUE);
    fProfileFourPartCorrel[iCent]->GetXaxis()->SetTitle("(m,n)");
    fMPCList->Add(fProfileFourPartCorrel[iCent]);

    if (!fWriteMinimum) // Cross-check is part of QA.
    {
      fProfileFourPartCorrelCheck[iCent] = new TProfile(
        Form("fProfileFourPartCorrelCheck_Bin%d", iCent),
        Form("#LT#LT4#GT#GT_{m,n,-m,-n}, Centrality bin %d", iCent),
          28, 0., 28.);
      fProfileFourPartCorrelCheck[iCent]->Sumw2();
      fProfileFourPartCorrelCheck[iCent]->SetStats(kTRUE);
      fProfileFourPartCorrelCheck[iCent]->GetXaxis()->SetTitle("(m,n)");
      fMPCList->Add(fProfileFourPartCorrelCheck[iCent]);
    } // End: if (!fWriteMinimum).

  // Profiles with the needed 6-particle correlators.
    fProfileSixPartCorrel[iCent] = new TProfile(Form("fProfileSixPartCorrel_Bin%d" ,iCent),
      Form("#LT#LT6#GT#GT_{m,n,o,-m,-n,-o}, Centrality bin %d", iCent),
      numberBinsSix, 0., numberBinsSix);
    fProfileSixPartCorrel[iCent]->Sumw2();
    fProfileSixPartCorrel[iCent]->SetStats(kTRUE);
    fMPCList->Add(fProfileSixPartCorrel[iCent]);

    if (fComputeACs)
    {
    // Profiles with the needed 8-particle correlators for ACs.
      fProfileEightPartCorrel[iCent] = new TProfile(Form("fProfileEightPartCorrel_Bin%d" ,iCent),
        Form("#LT#LT8#GT#GT_{m,n,o,p,-m,-n,-o,-p}, Centrality bin %d", iCent),
        3, 0., 3.);
      fProfileEightPartCorrel[iCent]->Sumw2();
      fProfileEightPartCorrel[iCent]->SetStats(kTRUE);
      fMPCList->Add(fProfileEightPartCorrel[iCent]);

    // Profiles with the needed 10-particle correlators for ACs.
      fProfileTenPartCorrel[iCent] = new TProfile(Form("fProfileTenPartCorrel_Bin%d" ,iCent),
        Form("#LT#LT10#GT#GT_{m,n,o,p,q,-m,-n,-o,-p,-q}, Centrality bin %d", iCent),
        1, 0., 1.);
      fProfileTenPartCorrel[iCent]->Sumw2();
      fProfileTenPartCorrel[iCent]->SetStats(kTRUE);
      fMPCList->Add(fProfileTenPartCorrel[iCent]);
    } // End: if (fComputeACs).
  } // End: iCent.
} // End: void BookMPCList().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisTaskTwoMultiCorrelations::BookTPCEtaList()
{
/* Book the profiles for the 2-particle correlators with eta gaps........................... */
  Float_t etaGaps[11] = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};
    // Possible values for the eta gap.

  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.
    for (Int_t i = 0; i < 11; i++)
    {
      fProfileTPCEta[iCent][i] = new TProfile("", "", 8, 0., 8.);
      fProfileTPCEta[iCent][i]->SetName(Form("fProfileTPCEta_EtaGap%d_Bin%d", i, iCent));
      fProfileTPCEta[iCent][i]->SetTitle(Form("#LT#LT2#GT#GT_{n,-n}, Centrality bin %d (#eta gap: %.1f)", iCent, etaGaps[i]));
      fProfileTPCEta[iCent][i]->SetStats(kTRUE);
      fProfileTPCEta[iCent][i]->Sumw2();
      fProfileTPCEta[iCent][i]->GetXaxis()->SetTitle("n");
      fTPCEtaList->Add(fProfileTPCEta[iCent][i]);

    // Set bin labels.
      for (Int_t gap = 1; gap < 9; gap++)
      {
       fProfileTPCEta[iCent][i]->GetXaxis()->SetBinLabel(gap, Form("%d", gap));
      }   // End of the loop to set the labels.
    }   // End of the loop to setup the profiles.
  } // End: iCent.
}   // End of "void BookTPCEtaList()".


