void FMDReconstruction (Int_t evNumber=1) 
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
  TFile * f = new TFile("mgalice.root","UPDATE");
  gAlice = (AliRun*) f->Get("gAlice") ;
  cout<<"gAlice="<<gAlice<<endl;
  AliFMD* FMD  = (AliFMD *)gAlice->GetDetector("FMD");
  cout<<"\nFMD="<<FMD<<endl;
  gAlice->RunReco("FMD") ;
}
  
