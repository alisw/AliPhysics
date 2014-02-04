//=================================================================//
//       Analysis AddTask for Multiplicity Fluctuation  Analysis 
//=================================================================//

TString fileNameBase="AnalysisResults.root";
//=======================================================//

AliAnalysisTask* AddAliEbyEMultFluctuationTask(Double_t vx, 
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
					     const char* taskstring) {

TString taskname = "MF";
  taskname.Append(taskstring);

AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMultFluctuations", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMultFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();

 AliEbyEMultFluctuationTask *taskMF = new AliEbyEMultFluctuationTask(taskstring);

 //taskMF->SetVertexDiamond(vx,vy,vz);
 /// taskMF->SetKinematicsCutsAOD(ptl,pth,eta);
 // taskMF->SetDCA(dcaxy, dcaz);
 // taskMF->SetCentralityEstimator(centralityEstimator);
 // taskMF->SetAnalysisType(analysisType);
 // taskMF->SetPtRange(ptl,pth);
 // taskMF->SetEta(eta);
 // taskMF->SetTPCNclus(nclus);
 // taskMF->SetAODtrackCutBit(AODfilterBit);
  
  taskMF->SelectCollisionCandidates(AliVEvent::kMB);

  mgr->AddTask(taskMF);
  

AliAnalysisDataContainer *coutFA = mgr->CreateContainer(taskname.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  fileNameBase.Data());
  mgr->ConnectInput(taskMF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMF, 1, coutFA);
  
  return taskMF;
}

  
