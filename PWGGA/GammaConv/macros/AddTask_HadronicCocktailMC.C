void AddTask_HadronicCocktailMC(Int_t particleFlag = 0, Bool_t runLightOutput = kFALSE, Double_t maxpTset = 50., Double_t binWidthPt = 0.05, TString maxyetaset = "0.80"/*80=0.80*/) {

    TObjArray *rConfigRapandEta = maxyetaset.Tokenize("_");
    if(rConfigRapandEta->GetEntries()<1){cout << "ERROR: AddTask_HadronicCocktailMC during parsing of maxyetaset String '" << maxyetaset.Data() << "'" << endl; return;}
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
    Error("AddTask_HadronicCocktailMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskHadronicCocktailMC *task=NULL;
  task = new AliAnalysisTaskHadronicCocktailMC(fEtaSet ? Form("HadronicCocktailMC_%1.2f_%1.2f",maxy,maxeta) : Form("HadronicCocktailMC_%1.2f",maxy));
  task->SetMaxY(maxy);
  if(fEtaSet)
    task->SetMaxEta(maxeta);
  task->SetLightOutput(runLightOutput);
  task->SetAnalyzedParticle(particleFlag);          // switch to run: 0 - pi0, 1 - eta, 2 - pi+-
  task->SetMaxPt(maxpTset);
  task->SetPtBinWidth(binWidthPt);

  TString                   analyzedParticle = "";
  if (particleFlag==0)      analyzedParticle = "pi0";
  else if (particleFlag==1) analyzedParticle = "eta";
  else if (particleFlag==2) analyzedParticle = "pi+-";

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(fEtaSet ? Form("HadronicCocktailMC_%s_%1.2f_%1.2f",analyzedParticle.Data(),maxy,maxeta) : Form("HadronicCocktailMC_%s_%1.2f",analyzedParticle.Data(),maxy), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:HadronicCocktailMC",AliAnalysisManager::GetCommonFileName()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
