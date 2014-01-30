//=========================================================================//
//                                                                         //
//           Analysis AddTask for Net-Charge Higher Moment Analysis        //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//                                                                         //
//=========================================================================//
TString fileNameBase="AnalysisResults.root";
//_________________________________________________________//
AliAnalysisTask*  AddAliEbyEHigherMomentsTaskCentrality(Double_t vx,
							Double_t vy,
							Double_t vz,
							Double_t ptl,
							Double_t pth,
							Int_t nptbins,
							Double_t eta,
							Int_t AODfilterBit = 128,
							const char* centralityEstimator,
							Bool_t trigger = kFALSE,
							TString  analysis,
							const char* taskss) {
  


   
  Bool_t IsMC = kFALSE;  
  if(analysis=="MCAOD") IsMC = kTRUE;
  
  TString taskname = "HM";
  taskname.Append(taskss);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();
  
  AliEbyEHigherMomentsTask *taskHM = new AliEbyEHigherMomentsTask("HigherMomentsTask");
  taskHM->SetVertexDiamond(vx,vy,vz);
  taskHM->SetCentralityEstimator(centralityEstimator);
  taskHM->SetAnalysisType(analysis);
  taskHM->SetAODtrackCutBit(AODfilterBit);
  taskHM->SetKinematicsCutsAOD(ptl,pth,eta);
  taskHM->SetNumberOfPtBins(nptbins);
  
  if(trigger) taskHM->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else taskHM->SelectCollisionCandidates(AliVEvent::kMB);
  // cout << " Check analysis type " << analysisType << endl;
  
  mgr->AddTask(taskHM);
 
  AliAnalysisDataContainer *coutFA = mgr->CreateContainer(taskname.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,fileNameBase.Data());
  mgr->ConnectInput(taskHM, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHM, 1, coutFA);
  
  return taskHM;
}
