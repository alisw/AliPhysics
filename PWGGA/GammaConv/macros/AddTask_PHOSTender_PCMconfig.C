AliPHOSTenderTask* AddTask_PHOSTender_PCMconfig(
                    const char* taskName        = "PHOSTenderTask",
                    const char* tenderName      = "PHOStender",
                    const char* options         = "",
                    Int_t pass                  = 1,
                    Bool_t isMC                 = kFALSE,
                    Int_t forceBadChannelMap    = 0, //0: no forced map, 1: forced OADB map, 2: forced single map file
                    TString specificBCMap       = "",
                    TString nonLinName          = "Default",                     // "Default", "Run2", "MC", "NoCorrection"
                    Bool_t isRun2               = kFALSE
)
{
  //Add a task with PHOS tender which works with AOD to the analysis train
  //Author: D.Peressounko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddPHOSTender_PCMConfig", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddPHOSTender_PCMConfig", "This task requires an input event handler");
    return NULL;
  }

  // create and add task
  AliPHOSTenderTask * tenderTask = new AliPHOSTenderTask(taskName) ;
  AliPHOSTenderSupply *PHOSSupply=new AliPHOSTenderSupply(tenderName) ;
  PHOSSupply->SetReconstructionPass(pass) ;
  PHOSSupply->SetNonlinearityVersion(nonLinName) ;

  tenderTask->SetPHOSTenderSupply(PHOSSupply) ;
  if(isMC) //handle MC data
    PHOSSupply->SetMCProduction(options) ;
  if (isRun2 && !isMC)
    PHOSSupply->ApplyZeroSuppression(0.020)

  if (forceBadChannelMap==1){
    std::cout << "=============================================================" << std::endl;
    std::cout << "INFO: AddPHOSTender_PCMConfig: "<< "You are setting a specific bad channel map using a full OADB file: " <<  specificBCMap.Data() << std::endl;
    std::cout << "=============================================================" << std::endl;
    PHOSSupply->SetPrivateOADBBadMap(specificBCMap.Data());
  }
  if (forceBadChannelMap==2){
    std::cout << "=============================================================" << std::endl;
    std::cout << "INFO: AddPHOSTender_PCMConfig: "<< "You are setting a specific bad channel independent of the run: " <<  specificBCMap.Data() << std::endl;
    std::cout << "=============================================================" << std::endl;
    PHOSSupply->ForceUsingBadMap(specificBCMap.Data());

  }
  //Need MagFeild
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(tenderTask);

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, mgr->GetCommonInputContainer());

  return tenderTask;
}
