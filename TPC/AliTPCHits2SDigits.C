Int_t AliTPCHits2SDigits(Int_t nevent=1)
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
  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }



  // gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

  TStopwatch timer;
  timer.Start();

  // uncomment lines below to set active sectors 
  //Int_t sec[10]={4,5,6,7,8,4+36,5+36,6+36,7+36,8+36};
  //TPC->SetActiveSectors(sec,10);

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d \n",eventn);
    gAlice->GetEvent(eventn);
    TPC->SetActiveSectors(); // all sectors set active
    printf("\nActive sectors\n");
    for (Int_t i=0;i<72;i++) if (TPC->IsSectorActive(i)) printf("%d\t",i);
    TPC->Hits2SDigits2(eventn);
  } 

  delete gAlice; gAlice=0;
  file->Close(); delete file;
  timer.Stop();
  timer.Print();

  return 0;
};




