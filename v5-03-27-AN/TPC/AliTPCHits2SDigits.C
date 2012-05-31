#ifndef __CINT__
  #include "alles.h"
  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliTPCtracker.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSRecPoint.h"
  #include "AliITSclusterV2.h"
  #include "AliITSsimulationFastPoints.h"
  #include "AliITStrackerV2.h"

#endif
Int_t AliTPCHits2SDigits(Int_t nevent=1)
{

  // new version by J.Belikov

  // Connect the Root Galice file containing Geometry, Kine and Hits

  //it assures full cleaning of prevous session
   if (gAlice)
    {
      delete AliRunLoader::Instance();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }

  AliRunLoader *rl = AliRunLoader::Open("galice.root","Event","update");
  if (!rl) 
   {
    cerr<<"Can't load RunLoader from "<<inFile_new<<" !\n";
    return 1;
   }

  // Get AliRun object from file or create it if not on file
  
  rl->LoadgAlice();
 
  gAlice = rl->GetAliRun();
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    delete rl;
    return 2;
  }

  // gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
  AliLoader * tpcl = rl->GetLoader("TPCLoader");
  if ((TPC == 0x0) || (tpcl == 0x0))
   {
    cerr<<"AliTPCHits2Digits.C : Can not find TPC or TPCLoader\n";
    delete rl;
    return 3;
   }
  
  tpcl->LoadHits("READ");
  tpcl->LoadSDigits("RECREATE");
  
  TStopwatch timer;
  timer.Start();

  // uncomment lines below to set active sectors 
  //Int_t sec[10]={4,5,6,7,8,4+36,5+36,6+36,7+36,8+36};
  //TPC->SetActiveSectors(sec,10);

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d \n",eventn);
    rl->GetEvent(eventn);
    TPC->SetTreeAddress();
    TPC->SetActiveSectors(); // all sectors set active
    printf("\nActive sectors\n");
    for (Int_t i=0;i<72;i++) if (TPC->IsSectorActive(i)) printf("%d\t",i);
    
    TPC->Hits2SDigits2(eventn);
  } 

  delete rl;
  timer.Stop();
  timer.Print();

  return 0;
};




