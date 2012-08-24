//=========================================================================//
//                                                                         //
//           Analysis AddTask for Net-Charge Higher Moment Analysis        //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//                               V0.0 23/08/2012                           //
//=========================================================================//

TString fileNameBase="AnalysisResults.root";

//  const char* analysisType        = "AOD"; // MC, ESD, AOD
// const char* centralityEstimator = "V0M"; // V0M, TRK, FMD, ....

//_________________________________________________________//

AliAnalysisTask* AddAliEbyEHigherMomentsTask(Double_t vx, 
					     Double_t vy, 
					     Double_t vz,
					     Double_t dcaxy, 
					     Double_t dcaz,
					     Double_t ptl, 
					     Double_t pth,
					     Double_t eta,
					     Int_t    nclus,
					     Int_t AODfilterBit, 
					     const char* centralityEstimator,
					     const char* analysisType,
					     const char* taskss) {
  
  
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
  taskHM->SetKinematicsCutsAOD(ptl,pth,eta);
  taskHM->SetDCA(dcaxy, dcaz);
  
  taskHM->SetCentralityEstimator(centralityEstimator);
  taskHM->SetAnalysisType(analysisType);
  taskHM->SetPtRange(ptl,pth);
  taskHM->SetEta(eta);
  taskHM->SetTPCNclus(nclus);
  taskHM->SetAODtrackCutBit(AODfilterBit);
  
  mgr->AddTask(taskHM);
  
  AliAnalysisDataContainer *coutFA = mgr->CreateContainer(taskname.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  fileNameBase.Data());
  mgr->ConnectInput(taskHM, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHM, 1, coutFA);
  
  return taskHM;
}
