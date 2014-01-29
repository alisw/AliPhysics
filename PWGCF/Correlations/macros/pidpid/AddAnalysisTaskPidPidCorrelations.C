// -*- c++ -*-
// $Id: AddAnalysisTaskPidPidCorrelations.C

const Int_t nBinCent = 1;
Double_t centLimits[nBinCent+1] = {0.,5.};
// Double_t centLimits[nBinCent+1] = {0.,1.,2.,3.,4.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.1};
const Int_t nBinZvtx = 10;
// Double_t zvtxLimits[nBinZvtx+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
Double_t zvtxLimits[nBinZvtx+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
const Int_t nBinPt = 8;
Double_t ptLimits[nBinPt+1] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0,, 5.0, 6.0, 7.0};
// Double_t ptLimits[nBinPt+1] = {0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0};
const Int_t nBinEta = 20;
Double_t etaLimits[nBinEta+1] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};


AliAnalysisTaskPidPidCorrelations*
AddAnalysisTaskPidPidCorrelations(TString centralityEstimator = "V0M"
				, Bool_t eventMixing = kTRUE
				, Double_t ptMin = 0.2
				, Double_t ptMax = 6.0
				, Double_t etaMin = -0.8
				, Double_t etaMax = 0.8
				, Double_t centrMin = 0.
				, Double_t centrMax = 5.
				, Double_t vertexZ = 10.
				, Int_t triggerPID = 1
				, Int_t assocPID = 1
// 				, Int_t AODfilterBit = 768 // PbPb LHC11h , AOD 142 
				, TString fileNameBase = "AnalysisResults"
				)
{
  TString centralityName = "";
  centralityName += Form("%s_cent%.0f%.0f_zvtx%.0f_trig%1d_assoc%1d",centralityEstimator.Data(),centrMin,centrMax,vertexZ,triggerPID,assocPID);
  TString outputFileName = fileNameBase + "_" + centralityName + "_" + (eventMixing ? "mix" : "nomix") + "_";
  
  //===========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr == NULL) { ::Error("AddTaskBalancePsiEfficiency", "No analysis manager to connect to."); return NULL; }
  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (mgr->GetInputEventHandler() == NULL) { ::Error("AddAnalysisTaskPidPidCorrelations", "This task requires an input event handler"); return NULL; }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType();
//   if (analysisType != "AOD") { ::Error("AddAnalysisTaskPidPidCorrelations", "This task runs only on AOD data"); return NULL; }
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
  
  outputFileName += analysisType.Data();
  outputFileName.Append(".root");
  
  
  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskPidPidCorrelations *taskPIDCorr = new AliAnalysisTaskPidPidCorrelations(Form("AliAnalysisTaskPidPidCorrelations_%s",centralityEstimator.Data()));

  taskPIDCorr -> SetTriggerMask(AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral);
//     taskPIDCorr->SelectCollisionCandidates(AliVEvent::kMB);//MB  //now inside task
//   taskPIDCorr -> SetFilterBit(AODfilterBit);
  
  //______ DATA or MC
  taskPIDCorr -> SetMC(kFALSE);

  //______ kinematics cut
  taskPIDCorr -> SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
//   taskPIDCorr -> SetEtaRange(-0.9, 0.9);
  //_____ vertex
  taskPIDCorr -> SetVertexDiamond(.3,.3,vertexZ);  
  //____ centrality
  taskPIDCorr -> SetCentralityEstimator(centralityEstimator);
  taskPIDCorr -> SetCentralityRange(centrMin,centrMax);
  taskPIDCorr -> SetEventMixing(eventMixing);
//     taskPIDCorr -> SetTriggerRestrictEta(0.5);
  taskPIDCorr -> SetEtaOrdering(kFALSE);
  taskPIDCorr -> SetPairCuts(kFALSE, kFALSE);
  taskPIDCorr -> SetTwoTrackEfficiencyCut(0.02,0.8);
  taskPIDCorr -> SetFillpT(kFALSE);
  taskPIDCorr -> SetMixingTracks(50000,1000);  
  taskPIDCorr -> SetWeightPerEvent(kFALSE);
  taskPIDCorr -> SetRejectResonanceDaughters(0);
  taskPIDCorr -> SetSelectCharge(0);
  taskPIDCorr -> SetSelectTriggerCharge(0);
  taskPIDCorr -> SetSelectAssociatedCharge(0);
  taskPIDCorr -> SetOnlyOneEtaSide(0);
  taskPIDCorr -> SetPtOrder(kFALSE);
  taskPIDCorr -> UseMomentumDifferenceCut(kFALSE,0.01);
  taskPIDCorr -> SetCentBinning(nBinCent, centLimits);
  taskPIDCorr -> SetZvtxBinning(nBinZvtx, zvtxLimits);
  taskPIDCorr -> SetPtBinning(nBinPt, ptLimits);
  taskPIDCorr -> SetEtaBinning(nBinEta, etaLimits);
  taskPIDCorr -> SetPIDsToCorrelate(triggerPID,assocPID);
  
  // ADD the task
  //===========================================================================
  mgr -> AddTask(taskPIDCorr);  
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
//   TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":outputPIDCorrAnalysis";
  AliAnalysisDataContainer *listPIDCorr1 = mgr -> CreateContainer( taskPIDCorr->GetOutputListName()/*Form("listPIDCorr_%s",centralityName.Data())*/, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
  AliAnalysisDataContainer *listPIDCorr2 = mgr -> CreateContainer( Form("OutputCFCont_%s",centralityName.Data()), AliCFContainer::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
//  AliAnalysisDataContainer *listPIDCorr3 = mgr -> CreateContainer( Form("MyCorrTTree_%s",centralityName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
  
  mgr -> ConnectInput( taskPIDCorr, 0, mgr -> GetCommonInputContainer());
  mgr -> ConnectOutput( taskPIDCorr, 1, listPIDCorr1);
  mgr -> ConnectOutput( taskPIDCorr, 2, listPIDCorr2);
//  mgr -> ConnectOutput( taskPIDCorr, 3, listPIDCorr3);
  
  return taskPIDCorr;
}
