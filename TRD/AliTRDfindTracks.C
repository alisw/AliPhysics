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

  Char_t *hitfile     = "galice.root"; 
  Char_t *clusterfile = "AliTRDclusters.root"; 
  Char_t *trackfile   = "AliTRDtracks.root"; 

  AliTRDtracker *Tracker =
    new AliTRDtracker("TheOnlyTRDtrackerSoFar","UniqueVariationOfIt");

  Tracker->GetEvent(hitfile,clusterfile);

  Tracker->Clusters2Tracks();

  Tracker->WriteTracks(trackfile);


}
