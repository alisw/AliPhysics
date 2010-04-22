// $Id$
/*
 * Example macro to run the HLT tracker embedded into AliRoot reconstruction.
 * The reconstruction is done from the TPC raw data. The TPCTrackHisto 
 * component fills histograms for the diagnostics of the tracking.
 *
 * Input is taken from TPC-clusters and TPC-globalmerger for the moment.
 *
 * Usage:
 * <pre>
 *
 * aliroot -b -q -l trackhisto.C'("file", "cdbURI", minEvent, maxEvent)' 
 * aliroot -b -q -l trackhisto.C'("raw://run115322", 0, 100)'
 * aliroot -b -q -l trackhisto.C'("raw.root", "local://$ALICE_ROOT/OCDB/", 0, 100)'
 * </pre>
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void trackhisto(const char *filename, const char *cdbURI, int minEvent=-1, int maxEvent=-1){

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Remove galice.root or run in a different folder." << endl;
    return;
  }

  if (!filename) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  
  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();

  AliHLTConfiguration trhiconf("trhi", "TPCTrackHisto", "TPC-clusters TPC-globalmerger", "");
  AliHLTConfiguration rwfconf("rfw", "ROOTFileWriter", "trhi", "-datafile TrackHisto -concatenate-events -overwrite");


  AliReconstruction rec;

  rec.SetRunQA(":");
  rec.SetInput(filename);
  
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri  = cdbURI;
  TString strfile = filename;
  if(struri.BeginsWith("raw://") ||  strfile.Contains("://") && !strfile.Contains("local://")){
     TGrid::Connect("alien");
  }
  
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);

  if(struri.BeginsWith("local://")) {
     // set specific storage for GRP entry
     // search in the working directory and one level above, the latter
     // follows the standard simulation setup like e.g. in test/ppbench
     if(!gSystem->AccessPathName("GRP/GRP/Data")) {
       man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
     } 
     else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
          man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");
     }
  }

  if(minEvent>=0 || maxEvent>minEvent){
     if(minEvent<0) minEvent=0;
     if(maxEvent<minEvent) maxEvent=minEvent;
     rec.SetEventRange(minEvent,maxEvent);
  }
 
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT"); 
  rec.SetOption("HLT", "loglevel=0x7c chains=rfw");

  rec.Run();

}

void trackhisto(const char *filename, int minEvent=-1, int maxEvent=-1){

  trackhisto(filename, "raw://", minEvent, maxEvent);
}

void trackhisto(){

  cout << "trackhisto: Run HLT TPC tracking and fill histograms" << endl;
  cout << " Usage: aliroot -b -q -l trackhisto.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
}

