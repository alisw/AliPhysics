const char* filename = "corr.root";

void
WriteFCM()
{
  AliForwardCorrectionManagerOADB& fcm = 
    AliForwardCorrectionManagerOADB::Instance();
  fcm.Init(0, 1, 900, 5);
  fcm.Print("R");

  TFile* file = TFile::Open(filename, "RECREATE");
  fcm.Write();
  delete file;
}

void 
ReadFCM()
{
  TFile* file = TFile::Open(filename, "READ");
  file->Get("forwardCorrections");

  AliForwardCorrectionManagerOADB& fcm = 
    AliForwardCorrectionManagerOADB::Instance();
  // fcm.Init(0, 1, 900, 5);
  fcm.Print("R");
}


void 
TestFCM(bool write=true)
{
  gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");

  if (write) WriteFCM();
  else       ReadFCM();
}

    
