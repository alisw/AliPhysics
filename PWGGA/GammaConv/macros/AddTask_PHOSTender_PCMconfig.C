AliPHOSTenderTask* AddTask_PHOSTender_PCMconfig(
                    const char* taskName        = "PHOSTenderTask",
                    const char* tenderName      = "PHOStender",
                    const char* options         = "",
                    Int_t pass                  = 1,
                    Bool_t isMC                 = kFALSE,
                    Int_t forceBadChannelMap    = 0, //0: no forced map, 1: forced OADB map, 2: forced single map file
                    TString specificBCMap       = "",
                    Bool_t useStandardPHOSNL    = kTRUE
)
{
  //Add a task with PHOS tender which works with AOD to the analysis train
  //Author: D.Peressounko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAODPHOSTender", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddAODPHOSTender", "This task requires an input event handler");
    return NULL;
  }

//  // input must be AOD
//  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
//  if( "AOD" != inputDataType )
//    ::Error("AddAODPHOSTender", Form("AOD input data required, input data is of type: %s", inputDataType.Data()));

  // create and add task
  AliPHOSTenderTask * tenderTask = new AliPHOSTenderTask(taskName) ;
  AliPHOSTenderSupply *PHOSSupply=new AliPHOSTenderSupply(tenderName) ;
  PHOSSupply->SetReconstructionPass(pass) ;
  if(!useStandardPHOSNL)
    PHOSSupply->SetNonlinearityVersion("NoCorrection") ;
  tenderTask->SetPHOSTenderSupply(PHOSSupply) ;
  if(isMC) //handle MC data
    PHOSSupply->SetMCProduction(options) ;
  if (forceBadChannelMap==1)
      PHOSSupply->SetPrivateOADBBadMap(specificBCMap.Data());
  if (forceBadChannelMap==2)
      PHOSSupply->ForceUsingBadMap(specificBCMap.Data());

  //Need MagFeild
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(tenderTask);

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, mgr->GetCommonInputContainer());

  return tenderTask;
}
