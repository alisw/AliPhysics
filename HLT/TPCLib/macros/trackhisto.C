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
 * aliroot -b -q -l trackhisto.C'("file", minEvent, maxEvent)'
 * </pre>
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void trackhisto(const char *filename, int minEvent=-1, int maxEvent=-1){

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
  AliHLTConfiguration rwfconf("rfw", "ROOTFileWriter", "trhi", "-datafile trackhisto -concatenate-events -overwrite");


  AliReconstruction rec;

  rec.SetRunQA(":");
  rec.SetInput(filename);

  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB/");  
  rec.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));

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

void trackhisto(){

  cout << "trackhisto: Run HLT TPC tracking and fill histograms" << endl;
  cout << " Usage: aliroot -b -q -l trackhisto.C'(\"file\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
}

