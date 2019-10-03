void AddTask_GammaCocktailMC(Bool_t runLightOutput = kFALSE, Double_t maxpTset = 50., Double_t binWidthPt = 0.05, TString maxyetaset = "80"/*80=0.80*/) {

  TObjArray *rConfigRapandEta = maxyetaset.Tokenize("_");
  if(rConfigRapandEta->GetEntries()<1){cout << "ERROR: AddTask_GammaCocktailMC during parsing of maxyetaset String '" << maxyetaset.Data() << "'" << endl; return;}
  TObjString* rmaxyset;
  TObjString* rmaxetaset;
  Bool_t fEtaSet = kFALSE;
  for(Int_t i = 0; i<rConfigRapandEta->GetEntries() ; i++){
    if(i==0)
      rmaxyset                = (TObjString*) rConfigRapandEta->At(i);
    else{
      rmaxetaset              = (TObjString*) rConfigRapandEta->At(i);
      fEtaSet                 = kTRUE;
    }
  }
  TString maxyset             = rmaxyset->GetString();
  Double_t maxy               = maxyset.Atof();
  maxy                        /= 100;  // needed to enable subwagon feature on grid
  cout << "running with max y cut of: " << maxy << endl;

  TString maxetaset;
  Double_t maxeta;
  if(fEtaSet){
    maxetaset                 = rmaxetaset->GetString();
    maxeta                    = maxetaset.Atof();
    maxeta                    /= 100;  // needed to enable subwagon feature on grid
    cout << "running with max eta cut of: " << maxeta << endl;
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_GammaCocktailMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaCocktailMC *task=NULL;
  task= new AliAnalysisTaskGammaCocktailMC(fEtaSet ? Form("GammaCocktailMC_%1.2f_%1.2f",maxy,maxeta) : Form("GammaCocktailMC_%1.2f",maxy));
  task->SetMaxY(maxy);
  if(fEtaSet)
    task->SetMaxEta(maxeta);
  task->SetLightOutput(runLightOutput);
  task->SetMaxPt(maxpTset);
  task->SetPtBinWidth(binWidthPt);
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(fEtaSet ? Form("GammaCocktailMC_%1.2f_%1.2f",maxy,maxeta) : Form("GammaCocktailMC_%1.2f",maxy), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaCocktailMC",AliAnalysisManager::GetCommonFileName()));
    
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
  
}
