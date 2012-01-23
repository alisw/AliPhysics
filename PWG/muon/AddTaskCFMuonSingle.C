AliCFMuonSingleTask1* AddTaskSingleMuonCF(Bool_t isMC=kTRUE)
{
  //
  // just a template of the AddTask macros
  // suppose we apply the cuts for same variables
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuonHFCorrAna", "No analysis manager to connect to.");
    return;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskMuonHFCorrAna", "MuonsHF task needs the manager to have an ESD or AOD input handler.");
    return;
  }
  if (type.Contains("ESD") && isMC && !mgr->GetMCtruthEventHandler()) {
    AliMCEventHandler *mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcHandler);
  }

  enum             { kEta,  kY, kPhi, kPt, kP3, kHit, kChi2Fit, kTrM, kChi2TrM,  kContrN,  kVt,  kVz, kTrig, kDCA, kZcoor, kRabs, kCharge, kTheta,kNVars };
  Int_t  nBins[] = {   5 , 5  ,  45 , 60 ,150 ,  20 ,      20 ,  4  ,      20 ,     202 , 100 ,  100 ,  10 , 500 ,   1000,    7 ,     3 , 100      };
  Double_t min[] = {  -4.,-4. ,-180.,  0.,  0.,   0.,       0., -0.5,       0.,     -2.5,   0., -100.,   0.,   0., -3000.,  171.,  -1.5 , 2.95      };
  Double_t max[] = {-2.5.,-2.5, 180., 30.,150.,  20.,      20.,  3.5,      10.,    199.5, 200.,  100.,  10., 500.,  1000.,  178.,   1.5 , 3.15     };

  Double_t *binLimits = 0;
  Int_t nSteps=1; if (isMC) nSteps=2;
  AliCFContainer* contCF = new AliCFContainer("container", "", nSteps, kNVars, nBins);
  for (Int_t var=kNVars; var--;) {
    binLimits = new Double_t[nBins[var]+1];
    for (Int_t i=0; i<=nBins[var]; i++) binLimits[i]=min[var]+i*(max[var]-min[var])/nBins[var];
    contCF->SetBinLimits(var, binLimits);
    delete [] binLimits; binLimits=0;
  }

  TList *qaList = new TList();
  TObjArray *genList = new TObjArray(0);
  AliCFTrackKineCuts *genKineCuts = new AliCFTrackKineCuts("genKineCuts", "GenKineCuts");
  genKineCuts->SetPtRange(min[kPt], max[kPt]);
  genKineCuts->SetRapidityRange(min[kY], max[kY]);
  genKineCuts->SetQAOn(qaList);
  genList->AddLast(genKineCuts);

  TObjArray *recList = new TObjArray(0);
  AliCFTrackKineCuts *recKineCuts = new AliCFTrackKineCuts("recKineCuts", "RecKineCuts");
  recKineCuts->SetPtRange(min[kPt], max[kPt]);
  recKineCuts->SetRapidityRange(min[kY], max[kY]);
  recKineCuts->SetQAOn(qaList);
  recList->AddLast(recKineCuts);

  AliCFManager* managerCF = new AliCFManager() ;
  managerCF->SetParticleContainer(contCF);
  managerCF->SetParticleCutsList(AliCFManager::kPartGenCuts, genList);
  managerCF->SetParticleCutsList(AliCFManager::kPartAccCuts, recList);

  AliCFMuonSingleTask1 *taskMuonCF = new AliCFMuonSingleTask1("AliMuonSingleTask1");
  taskMuonCF->SetCFManager(managerCF);
  taskMuonCF->SetQAList(qaList);
  taskMuonCF->SetUseMC(isMC);
  mgr->AddTask(taskMuonCF);

  mgr->ConnectInput(taskMuonCF, 0, mgr->GetCommonInputContainer());

  char *fileName = "muonCF.root";
  mgr->ConnectOutput(taskMuonCF,1,mgr->CreateContainer("chist",TH1I::Class(),AliAnalysisManager::kOutputContainer,fileName));
  mgr->ConnectOutput(taskMuonCF,2,mgr->CreateContainer("ccont",AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,fileName));

  return taskMuonCF;
}
