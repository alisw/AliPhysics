#ifndef __CINT__
  #include "alles.h"
#endif
Int_t AliTPCTestMerge(Int_t nevent=1)
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
  //if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }

  // gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

  TStopwatch timer;
  timer.Start();
  TTree  tree[100];   //
  Int_t mask[100];
  for(Int_t eventn =0;eventn<nevent;eventn++){
  
    printf("Processing event %d",eventn);
    char  name[100];
    sprintf(name,"TreeS_75x40_100x60_%d",eventn);
    tree[eventn].Read(name);
    tree[eventn].Dump();
    mask[eventn] = 2<<(20+eventn);
    
  } 
  TPC->Merge(tree,mask,nevent,nevent);
  
  file->Close(); delete file;
  timer.Stop();
  timer.Print();
  delete gAlice; gAlice=0;
  return 0;
};




