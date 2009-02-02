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

   if (gAlice)
    {
      delete AliRunLoader::Instance();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }

  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile_old);
  if (file) {file->Close(); delete file;}

  AliRunLoader *rl = AliRunLoader::Open(inFile_new,"Event","update");

  if (!rl) {
    cerr<<"Can't open "<<inFile_new<<" !\n";
    return 1;
  }

  // Get AliRun object from file or create it if not on file
  //  if (gAlice) delete gAlice;

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
//    delete rl;
    return 3;
   }

  tpcl->LoadSDigits("READ");
  tpcl->LoadDigits("RECREATE");
  rl->CdGAFile();
  AliTPCParamSR *dig=(AliTPCParamSR *)gDirectory->Get("75x40_100x60");
   if(dig){
     cerr<<"2 pad-length geom hits with 3 pad-lengths geom digits\n";
     delete dig;
     dig = new AliTPCParamSR();
   }
   else
   {
     dig=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
   }
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 4;}
   TPC->SetParam(dig);
       

  TStopwatch timer;
  timer.Start();

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d\n",eventn);
    cout<<"rl->GetEvent(eventn);\n";
    rl->GetEvent(eventn);
    cout<<"TPC->SDigits2Digits2(eventn);\n";
    TPC->SDigits2Digits2(eventn);
  } 

  timer.Stop();
  timer.Print();

  delete rl;
  return 0;
};












