
// Test Macro, shows how to execute the Trigger, and prints the results. 
// Author: Gustavo Conesa

void TestEMCALTrigger(){

  //Loader  
  AliRunLoader* rl=0x0;
  
  cout<<"TestEMCALTrigger: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "read");
  if (rl == 0x0)
    {
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    }

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));

  //Get Maximum number of events

  Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"n events "<<maxevent<<endl;

  //Load Digits  
  rl->LoadDigits("EMCAL");
  
  //Create trigger pointer and set thresholds if you want
  //Defautl threshold values need to be fixed
  AliEMCALTrigger *tr = new AliEMCALTrigger();  
  tr->SetL0MBPbPbThreshold(500);
  tr->SetL0MBppThreshold(100);
  tr->SetL1JetLowPtThreshold(2000);
  tr->SetL1JetMediumPtThreshold(10000);
  tr->SetL1JetHighPtThreshold(20000);

  //event loop
  for(Int_t iEvent = 0; iEvent < maxevent ; iEvent++){
    rl->GetEvent(iEvent);
    cout<<">>>>>>>>>>> Event >>> "<<iEvent<<endl;
    
    //Select trigger for each event
    tr->Trigger();//Do the trigger algorithm
    
    tr->Print("");//Print results. 
    
  }
}
