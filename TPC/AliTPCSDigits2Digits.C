#ifndef __CINT__
  #include "alles.h"
  #include "AliTPCtracker.h"
#endif
Int_t AliTPCSDigits2Digits(Int_t nevent=1)
{

  // new version by J.Belikov

  // Connect the Root Galice file containing Geometry, Kine and Hits

  const char * inFile_old = "galice.root"; 
  const char * inFile_new = "galice.root";
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile_old);
  if (file) {file->Close(); delete file;}
  file =  TFile::Open(inFile_new,"UPDATE");
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<inFile_new<<" !\n";
    return 1;
  }

  // Get AliRun object from file or create it if not on file
  //  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }



  // gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

  TStopwatch timer;
  timer.Start();

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d\n",eventn);
    gAlice->GetEvent(eventn);

    TPC->SDigits2Digits2(eventn);
  } 

  delete gAlice; gAlice=0;
  file->Close(); delete file;
  timer.Stop();
  timer.Print();

  return 0;
};












