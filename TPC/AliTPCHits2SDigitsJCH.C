#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "TFile.h"
#include "TTree.h"
#include "AliRun.h"
#include "AliTPC.h"
#endif


Int_t AliTPCHits2SDigitsJCH(TString fileNameSDigits="sdigits.root", TString fileNameHits="rfio:galice.root", Int_t nEvents = 1, Int_t firstEvent = 0)
{

  TString fileMode = "read";
  if (fileNameSDigits == fileNameHits || fileNameSDigits == "") fileMode = "update";
  TFile *fileHits =  TFile::Open(fileNameHits.Data(),fileMode.Data());
  if (!fileHits->IsOpen()) {
    cerr<<"Can't open "<<fileNameHits.Data()<<" !\n";
    return 1;
  }

  // Get AliRun object from file or create it if not on file
  if (gAlice) delete gAlice; gAlice = 0;
  gAlice = (AliRun*)fileHits->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }

  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

  TStopwatch timer;
  timer.Start();

  if (fileNameSDigits != fileNameHits && fileNameSDigits != "")  gAlice->InitTreeFile("S",fileNameSDigits.Data());

  // uncomment lines below to set active sectors 
  //Int_t sec[10]={4,5,6,7,8,4+36,5+36,6+36,7+36,8+36};
  //TPC->SetActiveSectors(sec,10);

  for(Int_t iEvent = firstEvent;iEvent<firstEvent+nEvents;iEvent++){
    printf("Processing event %d \n",iEvent);
    gAlice->GetEvent(iEvent);
//    TPC->SetActiveSectors(1); // all sectors set active
    printf("\nActive sectors\n");
    for (Int_t i=0;i<72;i++) {
      if (i%10 == 0) printf("\n");
      if (TPC->IsSectorActive(i)) printf("%2d ",i);
    }
    printf("\n");
    TPC->Hits2SDigits2(iEvent);
  } 

  timer.Stop();
  timer.Print();
  delete gAlice; gAlice=0;
  fileHits->Close(); delete fileHits;

  return 0;
};

