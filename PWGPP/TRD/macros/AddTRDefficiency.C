void AddTRDefficiency(AliAnalysisDataContainer **ci, Int_t effmc, Int_t mult)
{
  Info("AddTRDefficiency",  "[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\" [4]=\"%s\"", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName(), ci[4]->GetName()) ;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) return;

  //AliLog::SetClassDebugLevel("AliTRDefficiency", 5);
  AliAnalysisDataContainer *evInfoContainer = ci[3];
  AliTRDrecoTask *eff(NULL);
//        trackStatus = 0; // barrel tracks
//                    = 1; // ITS tracks
//                    = 2; // Kink tracks
  const Char_t *suffix[]={"", "ITS", "K"};
  for(Int_t its(0); its<1; its++){
    mgr->AddTask(eff = new AliTRDefficiency(Form("TRDefficiency%s", suffix[its])));
    eff->SetMCdata((Bool_t)mgr->GetMCtruthEventHandler());
    eff->SetDebugLevel(0);
    mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  // connect main (ESD) container
    mgr->ConnectInput(eff, 1, ci[its]);                         // conect track info container
    mgr->ConnectInput(eff, 2, evInfoContainer);                 // conect event info container
    mgr->ConnectInput(eff, 3, ci[4]);                           // conect onl.tracklets container
    mgr->ConnectInput(eff, 4, ci[5]);                           // conect clusters container
    mgr->ConnectOutput(eff,1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
  }

  // TRD combined tracking efficiency
  if(mgr->GetMCtruthEventHandler() && effmc) {
    mgr->AddTask(eff = new AliTRDefficiencyMC((char*)"TRDefficiencyMC"));
    eff->SetDebugLevel(0);
    //AliLog::SetClassDebugLevel("AliTRDefficiencyMC", 5);  

    // Create containers for input/output
    mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(eff, 1, ci[0]);
    mgr->ConnectInput(eff, 2, evInfoContainer);
    mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
  }

  // TRD single track selection
  if(!mult) return;

  mgr->AddTask(eff = new AliTRDmultiplicity((char*)"TRDmultiplicity"));
  eff->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDmultiplicity", 5);  
  mgr->ConnectInput(eff, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(eff, 1, ci[0]);
  mgr->ConnectOutput(eff, 1, mgr->CreateContainer(eff->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
}

