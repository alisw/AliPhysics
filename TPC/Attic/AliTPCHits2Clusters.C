/*
  .L AliTPCHits2Clusters.C+
  Hits2ExactClusters();

*/

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliTPCFast.h"


AliRunLoader *Init();

void Hits2ExactClusters(){
  AliRunLoader * runLoader = Init();
  AliTPCFast fast;
  
  for (Int_t ievent =0; ievent<runLoader->GetNumberOfEvents(); ievent++){
    runLoader->SetEventNumber(ievent);
    printf("Event\t%d\n",ievent);
    fast.Hits2ExactClusters(runLoader);
  }
}





AliRunLoader* Init(){
  //
  // initialization
  //
  if (gAlice) {
    delete AliRunLoader::Instance();
    delete gAlice;//if everything was OK here it is already NULL
    gAlice = 0x0;
  }
  
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return 0;
  }
  
  if (rl->LoadgAlice()) {
    cerr<<"Error occured while l"<<endl;
    return 0;
  }
   
  gAlice=rl->GetAliRun();
  if (!gAlice) {
    cerr<<"Can't get gAlice !\n";
    return 0;
  }
  return rl;
}



