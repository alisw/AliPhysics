/// \file AliTPCHits2Digits.C

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliTPC.h"

  #include "TStopwatch.h"
#endif

extern AliRun *gAlice;

Int_t AliTPCHits2Digits(Int_t nev=5) {
  /// Connect the Root Galice file containing Geometry, Kine and Hits

  if (gAlice) { 
     delete AliRunLoader::Instance();
     delete gAlice;//if everything was OK here it is already NULL
     gAlice = 0x0;
  }

  AliRunLoader *rl = AliRunLoader::Open("galice.root","Event","update");
  if (!rl) {
    cerr<<"Can't load RunLoader from "<<endl;
    return 1;
  }

  // Get AliRun object from file or create it if not on file

  rl->LoadgAlice();
 
  gAlice = rl->GetAliRun();
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }

  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      
  AliLoader * tpcl = rl->GetLoader("TPCLoader");
  if ((TPC == 0x0) || (tpcl == 0x0)) {
    cerr<<"AliTPCHits2Digits.C : Can not find TPC or TPCLoader\n";
    delete rl;
    return 3;
  }
  tpcl->LoadHits("READ");
  tpcl->LoadDigits("recreate");

  TStopwatch timer;
  timer.Start();

 // uncomment below lines to set sectors active
 // Int_t sec[10]={0,1,2,3,4,5,6,7,8,9};
 // TPC->SetActiveSectors(sec,10);

  for (Int_t i=0; i<nev; i++){
    printf("Processing event %d \n",i);
    if(rl->GetEvent(i)) break;
    TPC->SetActiveSectors(); // all sectors set active
    TPC->Hits2Digits(i);
  }

  delete rl;

  timer.Stop();
  timer.Print();

  return 0;
}

