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
				  Double_t dcaxy,
				  Double_t dcaz,
				  Double_t ptl,
				  Double_t pth,
				  Double_t eta,
				  Double_t rapidity,
				  Int_t    nclus,
				  Double_t chi2ndf,
				  TString particle,
				  Double_t nsigma,
				  Int_t AODfilterBit = 128,
				  const char* centralityEstimator,
				  Bool_t trigger = kFALSE,
				  Bool_t usepid,
				  Bool_t checkEff,
				  TString  analysis,
				  const char* taskss) {
  
  
  TString taskname = "HMQA";
  taskname.Append(taskss);
  TString taskname1 = "HM";
  taskname1.Append(taskss);
  
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
  taskHM->SetDCA(dcaxy, dcaz);
  taskHM->SetPtRange(ptl,pth);
  taskHM->SetEta(eta);
  taskHM->SetTPCNclus(nclus);
  taskHM->SetChi2PerNDF(chi2ndf);
  taskHM->SetAODtrackCutBit(AODfilterBit);
  taskHM->SetKinematicsCutsAOD(ptl,pth,eta);
  taskHM->SetUsePid(usepid);
  taskHM->SetEfficencyJob(checkEff);

  if(trigger) taskHM->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else taskHM->SelectCollisionCandidates(AliVEvent::kMB);
  
  if( usepid ){
    taskHM->SetNSigmaCut(nsigma);
    taskHM->SetRapidityCut(rapidity);
    if( particle == "Proton" ){  taskHM->SetParticleSpecies(AliPID::kProton); }
    else if( particle == "Kaon") { taskHM->SetParticleSpecies(AliPID::kKaon); }
    else if( particle == "Pion" ){ taskHM->SetParticleSpecies(AliPID::kPion); }
  }
  // cout << " Check analysis type " << analysisType << endl;
  
  mgr->AddTask(taskHM);
  
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(taskname.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,fileNameBase.Data());
  AliAnalysisDataContainer *coutFA = mgr->CreateContainer(taskname1.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,fileNameBase.Data());
  mgr->ConnectInput(taskHM, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHM, 1, coutQA);
  mgr->ConnectOutput(taskHM, 2, coutFA);
  
  return taskHM;
}
