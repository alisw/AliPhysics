Int_t AliTPCHits2Digits(Int_t nevent=1)
{

  // new version by J.Belikov

  // Connect the Root Galice file containing Geometry, Kine and Hits

  const char * inFile_old = "galice.root"; 
  const char * inFile_new = "rfio:galice.root";
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

  // uncomment below lines to set sectors active
  //Int_t sec[10]={0,1,2,3,4,5,6,7,8,9};
  //TPC->SetActiveSectors(sec,10);

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d \n",eventn);
    gAlice->GetEvent(eventn);
    TPC->SetActiveSectors(); // all sectors set active
    TPC->Hits2Digits(eventn);
  }

  delete gAlice; gAlice=0;
  file->Close(); delete file;
  timer.Stop();
  timer.Print();

  return 0;
};

