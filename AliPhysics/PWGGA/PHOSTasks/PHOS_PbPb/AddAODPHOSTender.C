AliPHOSTenderTask* AddAODPHOSTender(const char* taskName = "PHOSTenderTask",
				    const char* tenderName = "PHOStender",
				    const char* options = "",
				    Int_t pass = 1,
				    Bool_t isMC = kFALSE,
                                    const char* nonLinType = "",
                                    Double_t zsSimulation = 0
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
  tenderTask->SetPHOSTenderSupply(PHOSSupply) ;
  if(isMC){ //handle MC data
    PHOSSupply->SetMCProduction(options) ;
    PHOSSupply->ApplyZeroSuppression(zsSimulation);
  }
  PHOSSupply->SetNonlinearityVersion(nonLinType) ;
  
  //Need MagFeild
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(tenderTask);

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, mgr->GetCommonInputContainer());

  return tenderTask;
}
