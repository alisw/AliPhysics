void FMDReconstructor (Int_t evNumber=1) 
{
  if (gClassTable->GetID("AliRun") < 0) 
    {
      gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
      loadlibs();
    }
  if (gSystem->Getenv("CONFIG_SPLIT_FILE"))
    cout << "SPLIT" << endl;
  else
    cout << "NO SPLIT" << endl ;
    
  
  
  //TFile * f = new TFile("galice.root","UPDATE");
  //gAlice = (AliRun*) f->Get("gAlice") ;
  
  AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();

  AliFMD* FMD  = (AliFMD *)gAlice->GetDetector("FMD");

  gAlice->RunReco("FMD");
}
  
