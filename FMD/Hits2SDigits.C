void Hits2SDigits(){  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
    loadlibs();
  }
   gSystem->Setenv("CONFIG_SPLIT_FILE","1") ;
   if (gSystem->Getenv("CONFIG_SPLIT_FILE"))
    cout << "SPLIT" << endl;
   else
    cout << "NO SPLIT" << endl ;
  TFile * f = new TFile("galice.root","UPDATE");
  gAlice = (AliRun*) f->Get("gAlice") ;
  cout<<"gAlice="<<gAlice<<endl;
  AliFMD* FMD  = (AliFMD *)gAlice->GetDetector("FMD") ;
  gAlice->Hits2SDigits("FMD") ;
 }
