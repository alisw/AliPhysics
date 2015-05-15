
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
//  maxevent=5;
  //Load Digits  
  rl->LoadDigits("EMCAL");

  //event loop
  for(Int_t iEvent = 0; iEvent < maxevent ; iEvent++){
    rl->GetEvent(iEvent);
    cout<<">>>>>>>>>>> Event >>> "<<iEvent<<endl;
    AliEMCALTrigger *tr = new AliEMCALTrigger();
    //Create trigger pointer and set thresholds if you want
    //Default threshold values need to be fixed  
    //     tr->SetL0Threshold(10000);
    //     tr->SetL1JetLowPtThreshold(10000);
    //     tr->SetL1JetMediumPtThreshold(10000);
    //     tr->SetL1JetHighPtThreshold(10000);
    tr->SetPatchSize(4);//0 means 2x2, 1->4x4, 2->8x8, 3->16x16 ...
    //Select trigger for each event
    tr->Trigger();//Do the trigger algorithm
    //     cout<<"Patch "<<tr->GetPatchSize()<<endl;
    //     cout<<"Trigger patch "<< tr->GetPatchSize()
    // 	<<" 2x2 maximum amplitude sum "<<tr->Get2x2MaxAmplitude()
    //   	<<" nxn max amp sum "<<tr->GetnxnMaxAmplitude()<<endl;
    tr->Print("");//Print results. 
    
  }
}
