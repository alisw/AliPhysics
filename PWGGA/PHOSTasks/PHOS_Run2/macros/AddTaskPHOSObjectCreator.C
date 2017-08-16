AliAnalysisTaskPHOSObjectCreator* AddTaskPHOSObjectCreator(
    const char* name     = "PHOSObjectCreator",
    const UInt_t trigger = AliVEvent::kINT7|AliVEvent::kPHI7,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t isMC = kFALSE,
    const Double_t BunchSpace = 25,
    const Double_t distanceBC = 0,
    const TString period = "LHC15n"
    )
{
  //Add a task AliAnalysisTaskPHOSObjectCreator to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSObjectCreator", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSObjectCreator", "This task requires an input event handler");
    return NULL;
  }


  TString taskname = Form("%s_BS%dns_DBC%dcm",name,(Int_t)BunchSpace,(Int_t)(distanceBC*10));

  AliAnalysisTaskPHOSObjectCreator* task = new AliAnalysisTaskPHOSObjectCreator(taskname);
  task->SelectCollisionCandidates(trigger);
  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  task->SetBunchSpace(BunchSpace);//in unit of ns
  task->SetMinimumDistanceFromBC(distanceBC);

  if(isMC){  
    TF1 *f1nonlin = new TF1("f1nonlin","[2]*(1.+[0]*TMath::Exp(-x*x/2./[1]/[1]))",0,100);
    f1nonlin->SetNpx(1000);
    f1nonlin->SetParNames("a","b");
    //f1nonlin->FixParameter(0,0);//for core E
    //f1nonlin->FixParameter(1,1);//for core E
    //f1nonlin->FixParameter(0,-0.03);//for core E at ZS5
    //f1nonlin->FixParameter(1,  0.6);//for core E at ZS5
    //f1nonlin->FixParameter(0,-0.06); //for full E at ZS5 ;
    //f1nonlin->FixParameter(1,  0.4); //for full E at ZS5 ;
    if(period.Contains("LHC15n")){
      cout << "pp LHC15n is called." << endl;
      f1nonlin->FixParameter(0,-0.05); //for full E, ZS = 20MeV;
      f1nonlin->FixParameter(1,  0.6); //for full E, ZS = 20MeV;
      f1nonlin->FixParameter(2,1.006); //for full E, ZS = 20MeV;//decalib 3% on M123
      //f1nonlin->FixParameter(0,-0.03);//for core E at ZS 20 MeV
      //f1nonlin->FixParameter(1,  0.9);//for core E at ZS 20 MeV
    }
    else if(period.Contains("LHC15o")){
      cout << "PbPb LHC15o is called." << endl;
      f1nonlin->FixParameter(0,-0.04);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(1,  0.8);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(2, 0.99);//for core E at ZS 20 MeV with only MIP cut
    }
    task->SetUserDefinedNonlinearity(f1nonlin);
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());//event characerization
  //mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

