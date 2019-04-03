AliAnalysisTaskPHOSObjectCreator* AddTaskPHOSObjectCreator(
    const char* name           = "PHOSObjectCreator",
    const UInt_t trigger       = AliVEvent::kINT7|AliVEvent::kPHI7,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t isMC          = kFALSE,
    const Double_t BunchSpace  = 25.,
    const Bool_t NonLinCorr    = kTRUE,
    const Bool_t excludeM4     = kTRUE,
    const TString period       = "LHC15n",
    const Bool_t isSingleSim   = kFALSE
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

  TString taskname = Form("%s_BS%dns",name,(Int_t)BunchSpace);

  AliAnalysisTaskPHOSObjectCreator* task = new AliAnalysisTaskPHOSObjectCreator(taskname);
  if(!isSingleSim) task->SelectCollisionCandidates(trigger);
  task->SetTenderFlag(usePHOSTender);
  task->SetMCFlag(isMC);
  task->SetBunchSpace(BunchSpace);//in unit of ns
  task->ExcludeM4(excludeM4);

  if(isMC && NonLinCorr){
    TF1 *f1nonlin = new TF1("f1nonlin","[2]*(1.+[0]/(1. + TMath::Power(x/[1],2)))",0,100);
    f1nonlin->SetNpx(1000);
    f1nonlin->SetParNames("a","b (GeV)","c");

    if(period.Contains("LHC15n")){
      f1nonlin->FixParameter(0,-0.06); //for full E, ZS = 20MeV;
      f1nonlin->FixParameter(1,  0.7); //for full E, ZS = 20MeV;
      f1nonlin->FixParameter(2,1.012); //for full E, ZS = 20MeV;//decalib 3% on M123
    }
    else if(period.Contains("LHC15o")){
      f1nonlin->FixParameter(0,-0.06);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(1,  0.7);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(2,1.002);//for core E at ZS 20 MeV with only MIP cut
      //f1nonlin->FixParameter(0,-0.04);//for core E at ZS 20 MeV with only MIP cut
      //f1nonlin->FixParameter(1,  0.8);//for core E at ZS 20 MeV with only MIP cut
      //f1nonlin->FixParameter(2,0.993);//for core E at ZS 20 MeV with only MIP cut
    }
    else if(period.Contains("LHC17p") || period.Contains("LHC17q")){

      if(isSingleSim){
        f1nonlin->FixParameter(0,-0.06);//for core E at ZS 20 MeV with only MIP cut
        f1nonlin->FixParameter(1,  0.7);//for core E at ZS 20 MeV with only MIP cut
        f1nonlin->FixParameter(2,1.013);//for core E at ZS 20 MeV with only MIP cut
      }
      else{
        f1nonlin->FixParameter(0,-0.06);//for core E at ZS 20 MeV with only MIP cut
        f1nonlin->FixParameter(1,  0.7);//for core E at ZS 20 MeV with only MIP cut
        f1nonlin->FixParameter(2,1.013);//for core E at ZS 20 MeV with only MIP cut

      }

    }
    else{
      f1nonlin->FixParameter(0,-0.06);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(1,  0.7);//for core E at ZS 20 MeV with only MIP cut
      f1nonlin->FixParameter(2,1.000);//for core E at ZS 20 MeV with only MIP cut
    }
    task->SetUserDefinedNonlinearity(f1nonlin);

    printf("Non-linearity correction in M.C. is ON\n");
    printf("function shape : %s , where p0 = %4.3f , p1 = %4.3f , p2 = %4.3f\n",f1nonlin->GetTitle(),f1nonlin->GetParameter(0),f1nonlin->GetParameter(1),f1nonlin->GetParameter(2));

  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  //const TString listname = Form("hist_%s",name);
  const TString listname = Form("hist_%s",taskname.Data());

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",listname.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

