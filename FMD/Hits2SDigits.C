void Hits2SDigits(){  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
    loadlibs();
  }
  //    gSystem->Setenv("CONFIG_SPLIT_FILE","1") ;
   if (gSystem->Getenv("CONFIG_SPLIT_FILE"))
    cout << "SPLIT" << endl;
   else
    cout << "NO SPLIT" << endl ;
  TFile * f = new TFile("mgalice.root","UPDATE");
  gAlice = (AliRun*) f->Get("gAlice") ;
  cout<<"gAlice="<<gAlice<<endl;
  AliFMD* FMD  = (AliFMD *)gAlice->GetDetector("FMD") ;
  cout<<" FMD "<<FMD<<endl;
  gAlice->GetEvent(0);
  gAlice->Hits2SDigits() ;
  gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
  }
