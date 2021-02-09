AliPHOSTenderTask* AddTask_PHOSTender_PCMconfig(
                    TString taskName            = "PHOSTenderTask",
                    TString tenderName          = "PHOStender",
                    TString options             = "",
                    Int_t pass                  = 1,
                    Bool_t isMC                 = kFALSE,
                    Int_t forceBadChannelMap    = 0, //0: no forced map, 1: forced OADB map, 2: forced single map file
                    TString specificBCMap       = "",
                    TString nonLinName          = "Default",                     // "Default", "Run2", "MC", "NoCorrection"
                    Bool_t isRun2               = kFALSE,
                    Bool_t enableSpecialNL      = kFALSE
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

  std::cout << "=============================================================" << std::endl;
  std::cout << "======INFO: AddPHOSTender_PCMConfig:        =================" << std::endl;
  std::cout << "options: " << options.Data() << std::endl;
  std::cout << "pass: " << pass << std::endl;
  std::cout << "isMC: " << isMC << std::endl;
  std::cout << "forceBadChannelMap: " << forceBadChannelMap << std::endl;
  std::cout << "specificBCMap: " << specificBCMap.Data() << std::endl;
  std::cout << "nonLinName: " << nonLinName.Data() << std::endl;
  std::cout << "isRun2: " << isRun2 << std::endl;
  std::cout << "enableSpecialNL: " << enableSpecialNL << std::endl;
  std::cout << "=============================================================" << std::endl;


  // create and add task
  AliPHOSTenderTask * tenderTask = new AliPHOSTenderTask(taskName.Data()) ;
  AliPHOSTenderSupply *PHOSSupply=new AliPHOSTenderSupply(tenderName.Data()) ;
  PHOSSupply->SetReconstructionPass(pass) ;
  PHOSSupply->SetNonlinearityVersion(nonLinName) ;

  if(isMC){ //handle MC data
    if (enableSpecialNL){
      std::cout << "=============================================================" << std::endl;
      std::cout << "=========== setting specific NL =============================" << std::endl;
      std::cout << "=========== Period:" << options.Data() << "==================" << std::endl;
      std::cout << "=============================================================" << std::endl;
      if (options.CompareTo("LHC13b2") == 0){
          PHOSSupply->SetMCProduction(options.Data()) ;
          PHOSSupply->SetNonlinearityVersion("MC");
          Int_t n = 3;
          Double_t par[3]={ 0.993, 0.02, 1.6 } ;
          std::cout << "============== NL version: MC with parameters ======================" << std::endl;
          for (Int_t i  = 0; i< 3; i++){
            std::cout << "parameter " << i << ": \t" << par[i] << endl;
          }
          std::cout << "=============================================================" << std::endl;
          PHOSSupply->SetNonlinearityParams(n,par);
      } else if (options.CompareTo("MBMCLHC15o") == 0){
        PHOSSupply->SetMCProduction("Run2Default") ;
        PHOSSupply->SetNonlinearityVersion("Run2MC");//only for MC
        Int_t n = 3;
        Double_t par[3] = {1.002,-0.06,0.7};//only for MC
        std::cout << "============== NL version: MC with parameters ======================" << std::endl;
        for (Int_t i  = 0; i< 3; i++){
          std::cout << "parameter " << i << ": \t" << par[i] << endl;
        }
        PHOSSupply->SetNonlinearityParams(3,par);//only for MC
      }
    } else {
      PHOSSupply->SetMCProduction(options.Data()) ;
    }
  }
  if (isRun2 && isMC)
    PHOSSupply->ApplyZeroSuppression(0.020);


  if (forceBadChannelMap==1){
    std::cout << "=============================================================" << std::endl;
    std::cout << "INFO: AddPHOSTender_PCMConfig: "<< "You are setting a specific bad channel map using a full OADB file: " <<  specificBCMap.Data() << std::endl;
    std::cout << "=============================================================" << std::endl;
    PHOSSupply->SetPrivateOADBBadMap((char*)specificBCMap.Data());
  }
  if (forceBadChannelMap==2){
    std::cout << "=============================================================" << std::endl;
    std::cout << "INFO: AddPHOSTender_PCMConfig: "<< "You are setting a specific bad channel independent of the run: " <<  specificBCMap.Data() << std::endl;
    std::cout << "=============================================================" << std::endl;
    PHOSSupply->ForceUsingBadMap((char*)specificBCMap.Data());

  }
  tenderTask->SetPHOSTenderSupply(PHOSSupply) ;

  //Need MagFeild
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(tenderTask);

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, mgr->GetCommonInputContainer());

  return tenderTask;
}
