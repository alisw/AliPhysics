
// Test Macro, shows how to execute the Trigger, and prints the results. 
// Author: Gustavo Conesa

void TestPHOSTrigger(){

  //Loader  
  AliRunLoader* rl=0x0;
  
//   cout<<"TestPHOSTrigger: Creating Run Loader ..."<<endl;
//   rl = AliRunLoader::Open("galice.root",
// 			  AliConfig::GetDefaultEventFolderName(),
// 			  "read");
//   if (rl == 0x0)
//     {
//       gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
//       return;
//     }

//   AliPHOSLoader *emcalLoader = dynamic_cast<AliPHOSLoader*>
//     (rl->GetDetectorLoader("PHOS"));

//   //Load Digits  
//   rl->LoadDigits("PHOS");

  AliPHOSGetter * gime = AliPHOSGetter::Instance("./galice.root");
  //Get Maximum number of events
  Int_t maxevent = gime->MaxEvent();
  //Int_t maxevent =  rl->GetNumberOfEvents();
  cout<<"n events "<<maxevent<<endl;
//  maxevent=5;


  //event loop
  for(Int_t iEvent = 0; iEvent < maxevent ; iEvent++){
    //    rl->GetEvent(iEvent);
    gime->Event(iEvent,"D"); //Only Digits
    cout<<">>>>>>>>>>> Event >>> "<<iEvent<<endl;
    AliPHOSTrigger *tr = new AliPHOSTrigger();
    //Create trigger pointer and set thresholds if you want
    //Default threshold values need to be fixed  
    //     tr->SetL0Threshold(10000);
    //     tr->SetL1JetLowPtThreshold(10000);
    //     tr->SetL1JetHighPtThreshold(10000);
    tr->SetPatchSize(1);//0 means 2x2, 1->4x4, 2->8x8, 3->16x16 ...
    //Select trigger for each event
    tr->Trigger();//Do the trigger algorithm
    //cout<<"Patch "<<tr->GetPatchSize()<<endl;
    cout<<"Trigger patch "<< tr->GetPatchSize()
     	<<" 2x2 maximum amplitude sum "<<tr->Get2x2MaxAmplitude()
       	<<" nxn max amp sum "<<tr->GetnxnMaxAmplitude()<<endl;
    //tr->Print("");//Print results. 
    
  }
}
