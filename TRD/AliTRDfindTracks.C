void AliTRDfindTracks() {

//////////////////////////////////////////////////////////////////////////
//
// Reads TRD clusters from file, finds tracks, and stores them in the file 
//
//////////////////////////////////////////////////////////////////////////


  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }       

  Char_t *alifile = "AliTRDclusters.root"; 

  cerr<<"got this far"<<endl;

  AliTRDtracker *Tracker =
    new AliTRDtracker("TheOnlyTRDtrackerSoFar","UniqueVariationOfIt");


  Tracker->GetEvent(alifile);

  Int_t inner, outer, delta=60;
  for(Int_t i=0; i<1; i++) {
    outer=179-i; inner=outer-delta;
    Tracker->MakeSeeds(inner,outer);
  }

  Tracker->FindTracks();

  Tracker->WriteTracks();


}
