/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/
//////////////////////////////////////////////////////////////
//
//    Flow Qn Vector correction options:
//
//
//    Calibration method: 0: Calibration of Q/sqrt(M)
//                        1: Calibration of Q/M
//                        2: Calibration of Q/|Q|
//                        3: Calibration of Q
//
//    Calibration step  : 0: Raw
//                        1: Equalization
//                        2: Recentering
//                        3: Twist
//                        4: Scaling
//
//    Equalization method : 0: M/<M>
//                          1: 1+(M-<M>)/sigma(M)
//
//    Channel list      : Array of channel numbers that are included in Q-vector calculation
//
//    Twist and Scaling method: to be implemented
//                              0: Double harmonic track wise (advisable for TPC)
//                              1: Double harmonic Q wise
//                              2: Correlations
//
///////////////////////////////////////////////////////////////

#ifdef __ECLIPSE_IDE

#include <TTree.h>
#include <TSystem.h>
#include <TMath.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TObject.h>
#include <TFile.h>
#include <AliLog.h>
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsHistos.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsCutWithin.h"
#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsDetector.h"
#include "AliQnCorrectionsDetectorConfigurationTracks.h"
#include "AliQnCorrectionsDetectorConfigurationChannels.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsInputGainEqualization.h"
#include "AliQnCorrectionsQnVectorRecentering.h"
#include "AliQnCorrectionsQnVectorAlignment.h"
#include "AliQnCorrectionsQnVectorTwistAndRescale.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

#ifdef __CLING__
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsHistos.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsCutWithin.h"
#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsDetector.h"
#include "AliQnCorrectionsDetectorConfigurationTracks.h"
#include "AliQnCorrectionsDetectorConfigurationChannels.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsProfileCorrelationComponents.h"
#include "AliQnCorrectionsProfile3DCorrelations.h"
#include "AliQnCorrectionsInputGainEqualization.h"
#include "AliQnCorrectionsQnVectorRecentering.h"
#include "AliQnCorrectionsQnVectorAlignment.h"
#include "AliQnCorrectionsQnVectorTwistAndRescale.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisTaskQnVectorAnalysis.h"
// to load external macro under ROOT6
#include "AliForwardMCMultiplicityTask.h"
#include "AliForwardMultiplicityTask.h"
#include <PWGLF/FORWARD/analysis2/AddTaskForwardMult.C>
#endif


#include "runAnalysis.H"

using std::cout;
using std::endl;

void DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass);

void AddVZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddTPC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddTZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddFMDTaskForESDanalysis();
void AddRawFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddZDC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);
void AddSPD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager);

Int_t varForEventMultiplicity;

AliAnalysisDataContainer* AddTaskFlowQnVectorCorrections() {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowQnVectorCorrections", "No analysis manager found.");
    return 0;
  }


  AliQnCorrectionsManager *QnManager = new AliQnCorrectionsManager();
  AliAnalysisTaskFlowVectorCorrections *taskQnCorrections = new AliAnalysisTaskFlowVectorCorrections("FlowQnVectorCorrections");

  /* let's establish the event cuts for event selection */
  AliQnCorrectionsCutsSet *eventCuts = new AliQnCorrectionsCutsSet();
  eventCuts->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kVtxZ,zvertexMin,zvertexMax));
  if (bUseMultiplicity) {
    varForEventMultiplicity = AliQnCorrectionsVarManagerTask::kVZEROMultPercentile;
  }
  else {
    varForEventMultiplicity = AliQnCorrectionsVarManagerTask::kCentVZERO;
  }
  eventCuts->Add(new AliQnCorrectionsCutWithin(varForEventMultiplicity,centralityMin,centralityMax));
  taskQnCorrections->SetEventCuts(eventCuts);
  taskQnCorrections->SetUseOnlyCentCalibEvents(bUseOnlyCentCalibEvents);

  /* and the physics selection also */
  if (!b2015DataSet) {
    taskQnCorrections->SelectCollisionCandidates(AliVEvent::kMB);  // Events passing trigger and physics selection for analysis
  }
  else
    taskQnCorrections->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);  // Events passing trigger and physics selection for analysis

  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  histClass+= "TrackQA_NoCuts;";

  /* add the selected detectors */
  if (bUseTPC) {
    AddTPC(taskQnCorrections, QnManager);
    histClass+= "TrackQA_TPC;";
  }
  if (bUseSPD) {
    AddSPD(taskQnCorrections, QnManager);
    histClass+= "TrackletQA_SPD;";
  }
  if (bUseVZERO) {
    AddVZERO(taskQnCorrections, QnManager);
  }
  if (bUseTZERO) {
    AddTZERO(taskQnCorrections, QnManager);
  }
  if (bUseFMD) {
    AddFMD(taskQnCorrections, QnManager);
  }
  if (bUseRawFMD) {
    AddRawFMD(taskQnCorrections, QnManager);
  }
  if (bUseZDC) {
    AddZDC(taskQnCorrections, QnManager);
  }

  QnManager->SetShouldFillQnVectorTree(kFALSE);
  QnManager->SetShouldFillQAHistograms(kTRUE);
  QnManager->SetShouldFillNveQAHistograms(kTRUE);
  QnManager->SetShouldFillOutputHistograms(kTRUE);

  taskQnCorrections->SetFillExchangeContainerWithQvectors(kTRUE);
  taskQnCorrections->SetFillEventQA(kTRUE);

  taskQnCorrections->SetAliQnCorrectionsManager(QnManager);
  taskQnCorrections->DefineInOutput();
  taskQnCorrections->SetRunsLabels(&listOfRuns);

  /* let's handle the calibration file */
  cout << "=================== CALIBRATION FILE =============================================" << endl;
  TString inputCalibrationFilename = Form("%s/%s", szCorrectionsFilePath.Data(), szCorrectionsFileName.Data());
  if (szCorrectionsSource.EqualTo("local")) {
    cout << "\t File " << inputCalibrationFilename << endl << "\t being taken locally when building the task object" << endl;
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_local, inputCalibrationFilename.Data());
  }
  else if (szCorrectionsSource.EqualTo("aliensingle")) {
    cout << "\t File " << inputCalibrationFilename << " being taken from alien in the execution nodes" << endl;
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_aliensingle, inputCalibrationFilename.Data());
  }
  else if (szCorrectionsSource.EqualTo("alienmultiple")) {
    cout << "\t File " << inputCalibrationFilename << " being taken from alien in the execution nodes on a per run basis " << endl;
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_alienmultiple, inputCalibrationFilename.Data());
  }
  else if (szCorrectionsSource.EqualTo("OADBsingle")) {
    cout << "\t File " << inputCalibrationFilename << " being taken from OADB in the execution nodes" << endl;
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_OADBsingle, inputCalibrationFilename.Data());
  }
  else if (szCorrectionsSource.EqualTo("OADBmultiple")) {
    cout << "\t File " << inputCalibrationFilename << " being taken from OADB in the execution nodes on a per run basis " << endl;
    taskQnCorrections->SetCalibrationHistogramsFile(AliAnalysisTaskFlowVectorCorrections::CALIBSRC_OADBmultiple, inputCalibrationFilename.Data());
  }
  else {
    Error("AddTaskFlowQnVectorCorrections", "\t CALIBRATION FILE SOURCE NOT SUPPORTED. ABORTING!!!");
    return NULL;
  }
  cout << "==================================================================================" << endl;

  AliQnCorrectionsHistos* hists = taskQnCorrections->GetEventHistograms();
  DefineHistograms(QnManager, hists, histClass);


  mgr->AddTask(taskQnCorrections);

  mgr->ConnectInput(taskQnCorrections,  0, mgr->GetCommonInputContainer());

  //create output containers
  if (QnManager->GetShouldFillOutputHistograms()) {
    AliAnalysisDataContainer *cOutputHist =
      mgr->CreateContainer(QnManager->GetCalibrationHistogramsContainerName(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "CalibrationHistograms.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQn(), cOutputHist );
  }

  if (QnManager->GetShouldFillQnVectorTree()) {
    AliAnalysisDataContainer *cOutputQvec =
      mgr->CreateContainer("CalibratedQvector",
          TTree::Class(),
          AliAnalysisManager::kOutputContainer,
          "QvectorsTree.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotTree(), cOutputQvec );
  }

  if (QnManager->GetShouldFillQAHistograms()) {
    AliAnalysisDataContainer *cOutputHistQA =
      mgr->CreateContainer(QnManager->GetCalibrationQAHistogramsContainerName(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "CalibrationQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQA(), cOutputHistQA );
  }

  if (QnManager->GetShouldFillNveQAHistograms()) {
    AliAnalysisDataContainer *cOutputHistNveQA =
      mgr->CreateContainer(QnManager->GetCalibrationNveQAHistogramsContainerName(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "CalibrationQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistNveQA(), cOutputHistNveQA );
  }

  if (taskQnCorrections->GetFillEventQA()) {
    AliAnalysisDataContainer *cOutputQnEventQA =
      mgr->CreateContainer("QnEventQA",
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "QnEventQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotEventQA(), cOutputQnEventQA );
  }

  AliAnalysisDataContainer *cOutputQvecList =
    mgr->CreateContainer("CalibratedQvectorList",
        TList::Class(),
        AliAnalysisManager::kExchangeContainer,
        "QvectorsList.root");

  if (taskQnCorrections->GetFillExchangeContainerWithQvectors())
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotGetListQnVectors(), cOutputQvecList );

  return cOutputQvecList;
}

void AddVZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){


  Bool_t VZEROchannels[4][64];
  for(Int_t iv0=0; iv0<4; iv0++) for(Int_t ich=0; ich<64; ich++) VZEROchannels[iv0][ich] = kFALSE;

  for(Int_t ich=32; ich<64; ich++) VZEROchannels[0][ich] = kTRUE;  // channel list: kTRUE if channel should be used
  for(Int_t ich=0; ich<32; ich++) VZEROchannels[1][ich] = kTRUE;


  Int_t channelGroups[64];
  for(Int_t ich=0; ich<64; ich++) channelGroups[ich] = Int_t(ich / 8);

  //-----------------------------------------------------------
  // Our event classes for V0
  //
  const Int_t nVZEROdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nVZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the VZERO detector */
  AliQnCorrectionsDetector *VZERO = new AliQnCorrectionsDetector("VZERO", AliQnCorrectionsVarManagerTask::kVZERO);

  /* the VZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *VZEROAconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "VZEROA",
          CorrEventClasses,
          64, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROAconf->SetChannelsScheme(VZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kTRUE);
  VZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  VZEROAconf->AddCorrectionOnQnVector(alignA);
  /* lets configrure the QA histograms */
  VZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROAconf->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPC","VZEROC");
  /* now we add it to the detector configuration */
  VZEROAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROAconf);

  /* the VZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *VZEROCconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "VZEROC",
          CorrEventClasses,
          64, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROCconf->SetChannelsScheme(VZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kTRUE);
  VZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  VZEROCconf->AddCorrectionOnQnVector(alignC);
  /* lets configrure the QA histograms */
  VZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROCconf->SetQAMultiplicityAxis(100, 0.0, 500.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPC","VZEROA");
  /* now we add it to the detector configuration */
  VZEROCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(VZERO);
}

void AddTPC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){

  /////////////// Add TPC subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for TPC
  //
  const Int_t nTPCdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nTPCdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TPC  detector */
  AliQnCorrectionsDetector *TPC = new AliQnCorrectionsDetector("TPC", AliQnCorrectionsVarManagerTask::kTPC);

  /* the TPC detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks *TPCconf =
      new AliQnCorrectionsDetectorConfigurationTracks(
          "TPC",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScale = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScale->SetApplyTwist(kTRUE);
  twScale->SetApplyRescale(kFALSE);
  twScale->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  TPCconf->AddCorrectionOnQnVector(twScale);

  /* define the cuts to apply */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  AliQnCorrectionsCutsSet *cutsTPC = new AliQnCorrectionsCutsSet();
  if(!isESD){
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFilterBitMask768,0.5,1.5));
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta,-0.8,0.8));
    cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt,0.2,5.));
  }
  else {
    Bool_t UseTPConlyTracks=kFALSE;   // Use of TPC standalone tracks or Global tracks (only for ESD analysis)
    task->SetUseTPCStandaloneTracks(UseTPConlyTracks);
    if(UseTPConlyTracks){
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY,-3.0,3.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ,-3.0,3.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta,-0.8,0.8));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt,0.2,5.));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCnclsIter1,70.0,161.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2Iter1,0.2,4.0));
    }
    else{
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaXY,-0.3,0.3));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kDcaZ,-0.3,0.3));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kEta,-0.8,0.8));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kPt,0.2,5.));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCncls,70.0,161.0));
      cutsTPC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kTPCchi2,0.2,4.0));
    }
  }
  TPCconf->SetCuts(cutsTPC);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TPC);
}


void AddSPD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){

  /////////////// Add SPD subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for SPD
  //
  const Int_t nSPDdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nSPDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the SPD detector */
  AliQnCorrectionsDetector *SPD = new AliQnCorrectionsDetector("SPD", AliQnCorrectionsVarManagerTask::kSPD);

  /* the SPD detector configuration */
  AliQnCorrectionsDetectorConfigurationTracks *SPDconf =
      new AliQnCorrectionsDetectorConfigurationTracks(
          "SPD",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  SPDconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  SPDconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector twist correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScale = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScale->SetApplyTwist(kTRUE);
  twScale->SetApplyRescale(kFALSE);
  twScale->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
  SPDconf->AddCorrectionOnQnVector(twScale);

  /* add the configuration to the detector */
  SPD->AddDetectorConfiguration(SPDconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(SPD);
}

void AddTZERO(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){

  /////////////// Add TZERO subdetectors ///////////////////

  Bool_t TZEROchannels[2][24];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<24; ich++) TZEROchannels[iv0][ich] = kFALSE;

  for(Int_t ich=12; ich<24; ich++) TZEROchannels[0][ich] = kTRUE;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<12; ich++) TZEROchannels[1][ich] = kTRUE;

  Int_t channelGroups[24];
  for(Int_t ich=0; ich<24; ich++) channelGroups[ich] = Int_t(ich / 12);

  //-----------------------------------------------------------
  // Our event classes for TZERO
  //
  const Int_t nTZEROdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nTZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TZERO detector */
  AliQnCorrectionsDetector *TZERO = new AliQnCorrectionsDetector("TZERO", AliQnCorrectionsVarManagerTask::kTZERO);

  /* the TZEROA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *TZEROAconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "TZEROA",
          CorrEventClasses,
          24, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROAconf->SetChannelsScheme(TZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kFALSE);
  TZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  TZEROAconf->AddCorrectionOnQnVector(alignA);
  /* let's configure the QA histograms */
  TZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROAconf->SetQAMultiplicityAxis(100, 0.0, 150.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPC","TZEROC");
  /* now we add it to the detector configuration */
  TZEROAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROAconf);

  /* the TZEROC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *TZEROCconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "TZEROC",
          CorrEventClasses,
          24, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROCconf->SetChannelsScheme(TZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kFALSE);
  TZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  TZEROCconf->AddCorrectionOnQnVector(alignC);
  /* let's configure the QA histograms */
  TZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROCconf->SetQAMultiplicityAxis(100, 0.0, 150.0);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPC","TZEROA");
  /* now we add it to the detector configuration */
  TZEROCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TZERO);
}


void AddZDC(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){
  /////////////// Add ZDC subdetectors ///////////////////

  Bool_t ZDCchannels[2][10];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = kFALSE;

  for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = kTRUE;
  for(Int_t ich=1; ich<5; ich++)  ZDCchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for ZDC
  //
  const Int_t nZDCdim = 3;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nZDCdim);
  Double_t VtxXbinning[][2] = {{ -0.3, 2}, {0.3, 10 }};
  Double_t VtxYbinning[][2] = {{ -0.3, 2}, {0.3, 10 }};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxX,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxX), VtxXbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxY,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxY), VtxYbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the ZDC detector */
  AliQnCorrectionsDetector *ZDC = new AliQnCorrectionsDetector("ZDC", AliQnCorrectionsVarManagerTask::kZDC);

  /* the ZDCA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *ZDCAconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "ZDCA",
          CorrEventClasses,
          10, /* number of channels */
          3); /* number of harmonics: 1, 2 and 3 */
  ZDCAconf->SetChannelsScheme(ZDCchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCAconf);

  /* the ZDCC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *ZDCCconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "ZDCC",
          CorrEventClasses,
          10, /* number of channels */
          3); /* number of harmonics: 1, 2 and 3 */
  ZDCCconf->SetChannelsScheme(ZDCchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(ZDC);
}

void AddFMDTaskForESDanalysis(){

  gSystem->Load("libPWGLFforward2");  // for FMD

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // Create the FMD task and add it to the manager
  //===========================================================================
  //--- AOD output handler -----------------------------------------
  AliAODHandler* ret = new AliAODHandler();

  ret->SetOutputFileName("AliAOD.pass2.root");
  mgr->SetOutputEventHandler(ret);

#ifndef __CLING__
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");
#endif  

  ULong_t run = 0; // 0: get from data???
  UShort_t sys = 0; // 0: get from data, 1: pp, 2: AA
  UShort_t sNN = 0; // 0: get from data, otherwise center of mass energy (per nucleon pair)
  Short_t  fld = 0; // 0: get from data, otherwise L3 field in kG



  const Char_t* config = "$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/ForwardAODConfig2.C";
  AliAnalysisTask *taskFmd  = AddTaskForwardMult(bMC, run, sys, sNN, fld, config,0,0);

  // --- Make the output container and connect it --------------------
  AliAnalysisDataContainer* histOut =
    mgr->CreateContainer("Forward", TList::Class(),
        AliAnalysisManager::kExchangeContainer);

  AliAnalysisDataContainer *output =
    mgr->CreateContainer("ForwardResultsP", TList::Class(),
        AliAnalysisManager::kExchangeContainer);

  mgr->ConnectInput(taskFmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskFmd, 1, histOut);
  mgr->ConnectOutput(taskFmd, 2, output);

}

void AddFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  if(isESD) AddFMDTaskForESDanalysis();


  Bool_t FMDchannels[2][4000];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<4000; ich++) FMDchannels[iv0][ich] = kFALSE;

  for(Int_t ich=2000; ich<4000; ich++) FMDchannels[0][ich] = kTRUE;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<2000; ich++) FMDchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the FMD detector */
  AliQnCorrectionsDetector *FMD = new AliQnCorrectionsDetector("FMD", AliQnCorrectionsVarManagerTask::kFMD);

  /* the FMDA detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *FMDAconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "FMDA",
          CorrEventClasses,
          4000, /* number of channels */
          4); /* number of harmonics: 1, 2 and 3 */
  FMDAconf->SetChannelsScheme(FMDchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDAconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDAconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  FMDAconf->AddCorrectionOnQnVector(alignA);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleA = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleA->SetApplyTwist(kTRUE);
  twScaleA->SetApplyRescale(kTRUE);
  twScaleA->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleA->SetReferenceConfigurationsForTwistAndRescale("TPC","FMDC");
  /* now we add it to the detector configuration */
  FMDAconf->AddCorrectionOnQnVector(twScaleA);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDAconf);

  /* the FMDC detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *FMDCconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "FMDC",
          CorrEventClasses,
          4000, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDCconf->SetChannelsScheme(FMDchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDCconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDCconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  FMDCconf->AddCorrectionOnQnVector(alignC);
  /* let's configure the twist and rescale correction step */
  AliQnCorrectionsQnVectorTwistAndRescale *twScaleC = new AliQnCorrectionsQnVectorTwistAndRescale();
  twScaleC->SetApplyTwist(kTRUE);
  twScaleC->SetApplyRescale(kTRUE);
  twScaleC->SetTwistAndRescaleMethod(AliQnCorrectionsQnVectorTwistAndRescale::TWRESCALE_correlations);
  twScaleC->SetReferenceConfigurationsForTwistAndRescale("TPC","FMDA");
  /* now we add it to the detector configuration */
  FMDCconf->AddCorrectionOnQnVector(twScaleC);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMD);
}


void AddRawFMD(AliAnalysisTaskFlowVectorCorrections *task, AliQnCorrectionsManager* QnManager){

  /////////////// Add FMD subdetectors ///////////////////
  /* FMD1 and FMD2 make FMDA and FMD3 make FMDC */

  const Int_t nNoOfDetectors       = 3;           ///< the number of FMD detectors
  const Int_t detectorNumber[]     = {1,2,3};     ///< the number of the FMD detector
  const Int_t nNoOfRings[]         = {1,2,2};     ///< the number of rings for each detector
  const Int_t ringNoOfSectors[]    = {20,40};   ///< ring number of sectors
  const Int_t nTotalNoOfChannels   = ringNoOfSectors[0] * 3 + ringNoOfSectors[1] * 2; ///< three inner sectors plus two outer ones
  const Int_t nTotalNoOfGroups     = nNoOfRings[0] + nNoOfRings[1] + nNoOfRings[2]; ///< each ring one channel group
  const Int_t FMDCdetectorNumber   = 3;           ///< the number of the detector associated to FMDC
  Int_t nSectorId = 0;                            ///< the sector id used as channel number
  Int_t nRingId = 0;                              ///< the ring id (0..4) used as group number

  Bool_t FMDchannels[2][nTotalNoOfChannels];      ///< the assignment of channels to each subdetector
  Int_t FMDchannelGroups[nTotalNoOfChannels];     ///< the group associated to each channel
  for (Int_t i = 0; i < 2; i++) for (Int_t c = 0; c < nTotalNoOfChannels; c++) FMDchannels[i][c] = kFALSE;

  for(Int_t detector = 0; detector < nNoOfDetectors; detector++) {
    for(Int_t ring = 0; ring < nNoOfRings[detector]; ring++) {
      for(Int_t sector = 0; sector < ringNoOfSectors[ring]; sector++) {
        FMDchannels[0][nSectorId] = ((detectorNumber[detector] != FMDCdetectorNumber) ? kTRUE : kFALSE);
        FMDchannels[1][nSectorId] = ((detectorNumber[detector] != FMDCdetectorNumber) ? kFALSE : kTRUE);
        FMDchannelGroups[nSectorId] = nRingId;
        nSectorId++;
      }
      nRingId++;
    }
  }

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  AliQnCorrectionsEventClassVariablesSet *CorrEventClasses = new AliQnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(AliQnCorrectionsVarManagerTask::kVtxZ,
      task->VarName(AliQnCorrectionsVarManagerTask::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new AliQnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  AliQnCorrectionsCutsSet *cutFMDA = new AliQnCorrectionsCutsSet();
  cutFMDA->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFMDEta,0.0,6.0));

  AliQnCorrectionsCutsSet *cutFMDC = new AliQnCorrectionsCutsSet();
  cutFMDC->Add(new AliQnCorrectionsCutWithin(AliQnCorrectionsVarManagerTask::kFMDEta,-6.0,0.0));

  /* the FMD detector */
  AliQnCorrectionsDetector *FMDraw = new AliQnCorrectionsDetector("FMDraw", AliQnCorrectionsVarManagerTask::kFMDraw);

  /* the FMDAraw detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *FMDArawconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "FMDAraw",
          CorrEventClasses,
          nTotalNoOfChannels,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDArawconf->SetChannelsScheme(FMDchannels[0], FMDchannelGroups);
  /* let's configure the Q vector calibration */
  FMDArawconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqA = new AliQnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kTRUE);
  FMDArawconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDArawconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignA = new AliQnCorrectionsQnVectorAlignment();
  alignA->SetHarmonicNumberForAlignment(2);
  alignA->SetReferenceConfigurationForAlignment("TPC");
  FMDArawconf->AddCorrectionOnQnVector(alignA);
  /* and add the cuts */
  FMDArawconf->SetCuts(cutFMDA);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDArawconf);

  /* the FMDCraw detector configuration */
  AliQnCorrectionsDetectorConfigurationChannels *FMDCrawconf =
      new AliQnCorrectionsDetectorConfigurationChannels(
          "FMDCraw",
          CorrEventClasses,
          nTotalNoOfChannels,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDCrawconf->SetChannelsScheme(FMDchannels[1], FMDchannelGroups);
  /* let's configure the Q vector calibration */
  FMDCrawconf->SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QVNORM_QoverM);
  /* lets configure the equalization of input data */
  AliQnCorrectionsInputGainEqualization *eqC = new AliQnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(AliQnCorrectionsInputGainEqualization::GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kTRUE);
  FMDCrawconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDCrawconf->AddCorrectionOnQnVector(new AliQnCorrectionsQnVectorRecentering());
  /* let's add the Q vector alignment correction step */
  AliQnCorrectionsQnVectorAlignment *alignC = new AliQnCorrectionsQnVectorAlignment();
  alignC->SetHarmonicNumberForAlignment(2);
  alignC->SetReferenceConfigurationForAlignment("TPC");
  FMDCrawconf->AddCorrectionOnQnVector(alignC);
  /* and add the cuts */
  FMDCrawconf->SetCuts(cutFMDC);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDCrawconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMDraw);
}

//__________________________________________________________________
void DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  const Char_t* histClasses = histClass.Data();

  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClass<< endl;

  //fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};

  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;

    // Event wise histograms
    if(classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,AliQnCorrectionsVarManagerTask::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,AliQnCorrectionsVarManagerTask::kIsPhysicsSelection, 0,0.0,0.0,AliQnCorrectionsVarManagerTask::kNothing, 0,0.0,0.0,AliQnCorrectionsVarManagerTask::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-30.0,30.0,AliQnCorrectionsVarManagerTask::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,AliQnCorrectionsVarManagerTask::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManagerTask::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManagerTask::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsMultPVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kVZEROMultPercentile, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100,0.,100., AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets, 100,0.,100., AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100,0.,100., AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy, 100,0.,100., AliQnCorrectionsVarManagerTask::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManagerTask::kNtracksSelected, 100, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult, 100, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kSPDntracklets, 100, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, AliQnCorrectionsVarManagerTask::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 9500.0, AliQnCorrectionsVarManagerTask::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 16000.0, AliQnCorrectionsVarManagerTask::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, AliQnCorrectionsVarManagerTask::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManagerTask::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManagerTask::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManagerTask::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManagerTask::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManagerTask::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManagerTask::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManagerTask::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManagerTask::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"MultPercentVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, AliQnCorrectionsVarManagerTask::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManagerTask::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,AliQnCorrectionsVarManagerTask::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManagerTask::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManagerTask::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManagerTask::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManagerTask::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManagerTask::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManagerTask::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManagerTask::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManagerTask::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC);
      continue;
    }  // end if className contains "Event"


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",AliQnCorrectionsVarManagerTask::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTrigger, 2, -0.5, 1.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentTPC, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentSPD, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentZDC, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManagerTask::kOfflineTriggerFired2, 200, -20.0, 20.0, AliQnCorrectionsVarManagerTask::kVtxZ, 0, 0.0, 0.0, AliQnCorrectionsVarManagerTask::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManagerTask::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, AliQnCorrectionsVarManagerTask::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, AliQnCorrectionsVarManagerTask::kTPCncls);
      histos->AddHistogram(classStr.Data(), "TPCsa_TPCncls", "TPC standalone TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, AliQnCorrectionsVarManagerTask::kTPCnclsIter1);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, 0.0, 50.0, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -1.5, 1.5, AliQnCorrectionsVarManagerTask::kEta);
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManagerTask::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManagerTask::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 100, 0.0, 2.2, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 100, 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManagerTask::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManagerTask::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManagerTask::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManagerTask::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManagerTask::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManagerTask::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManagerTask::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManagerTask::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManagerTask::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManagerTask::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManagerTask::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManagerTask::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManagerTask::kSinNPhi+ih);
      }
    }

    // Tracklet histograms
    if(classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -3.0, 3.0, AliQnCorrectionsVarManagerTask::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManagerTask::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -3.0, +3.0, AliQnCorrectionsVarManagerTask::kSPDtrackletEta, 100, 0.0, 6.3, AliQnCorrectionsVarManagerTask::kSPDtrackletPhi);
    }

  }

  cout << " done" << endl;
}

