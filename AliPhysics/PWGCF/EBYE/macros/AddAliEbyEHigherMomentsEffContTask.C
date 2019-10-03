//=========================================================================//
//                                                                         //
//           Analysis AddTask for Net-Charge Higher Moment Analysis        //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//                                                                         //
//=========================================================================//
TString fileNameBase="AnalysisResults.root";

//const char* analysisType        = "MCAOD"; // MC, ESD, AOD
//const char* centralityEstimator = "V0M"; // V0M, TRK, FMD, ....

//_________________________________________________________//

AliAnalysisTask* AddAliEbyEHigherMomentsEffContTask(Double_t vx,
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
						    Int_t AODfilterBit = 128,
						    Bool_t usepid,
						    TString particle,
						    Double_t nsigma,
						    Bool_t CheckCont,
						    Bool_t efficiency,
						    const char* centralityEstimator,
						    Bool_t trigger = kFALSE,
						    const char* analysisType,
						    const char* taskss) {
  
  
  TString taskname = "HMQA";
  TString taskname1 = "HM";
  taskname.Append(taskss);
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
  
  AliEbyEHigherMomentsEffContTask *taskHM = new AliEbyEHigherMomentsEffContTask("HigherMomentsTask");
  taskHM->SetVertexDiamond(vx,vy,vz);
  taskHM->SetKinematicsCutsAOD(ptl,pth,eta);
  taskHM->SetDCA(dcaxy, dcaz);
  taskHM->SetCentralityEstimator(centralityEstimator);
  taskHM->SetAnalysisType(analysisType);
  taskHM->SetPtRange(ptl,pth);
  taskHM->SetEta(eta);
  taskHM->SetTPCNclus(nclus);
  taskHM->SetChi2PerNDF(chi2ndf);
  taskHM->SetAODtrackCutBit(AODfilterBit);
  taskHM->SetUsePid(usepid);
  taskHM->SetContaMinationCheck(CheckCont);
  taskHM->SetEfficencyJob(efficiency);
  
  if( usepid ){
    taskHM->SetNSigmaCut(nsigma);
    taskHM->SetRapidityCut(rapidity);
    if( particle == "Proton" ){  taskHM->SetParticleSpecies(AliPID::kProton); }
    else if( particle == "Kaon") { taskHM->SetParticleSpecies(AliPID::kKaon); }
    else if( particle == "Pion" ){ taskHM->SetParticleSpecies(AliPID::kPion); }
  }
  
  if(trigger) taskHM->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else taskHM->SelectCollisionCandidates(AliVEvent::kMB);
  
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
