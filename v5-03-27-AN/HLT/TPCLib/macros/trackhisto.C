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
 * aliroot -b -q -l trackhisto.C'("raw.root", "local://$ALICE_ROOT/OCDB/", 0, 100, kTRUE)'
 * for when we want to have an output file with tracklet properties.
 * </pre>
 *
 * In the latter case an extra converter component is attached to the output
 * of CATracker and produces the exact output datatype as the global merger.
 * For this reason it cannot be run in parallel to the normal chain
 * with the merged tracks. We need to direct these special tracklet data 
 * blocks differently and run the macro with kTRUE as the last argument.
 *
 * This last argument is implemented only for debugging and QA reasons. 
 * The standard user does not need it.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void trackhisto(const char *filename, const char *cdbURI, int minEvent=-1, int maxEvent=-1, Bool_t bTracklets=kFALSE){

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Remove galice.root or run in a different folder." << endl;
    return;
  }

  if (!filename) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  

  ////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  ////////////////////////////////////////////////////////////////////
 
  
  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();
 
  int iMinSlice =  0;
  int iMaxSlice = 35;
  int iMinPart  =  0;
  int iMaxPart  =  5;
  
  TString tracklets;
  
  for(int slice=iMinSlice; slice<=iMaxSlice; slice++){
      
      TString trackerInput;      
      for(int part=iMinPart; part<=iMaxPart; part++){
          
	  TString arg, publisher, cf;
          
          // raw data publisher components
          int ddlno=768;
          if (part>1) ddlno+=72+4*slice+(part-2);
          else ddlno+=2*slice+part;
      
          arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
          publisher.Form("DP_%02d_%d", slice, part);
          AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

          // cluster finder components
          cf.Form("CF_%02d_%d", slice, part);
          AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinder32Bit", publisher.Data(), "");
     
          if (trackerInput.Length()>0) trackerInput+=" ";
          trackerInput+=cf;     
    }
    
    TString tracker, converted_tracklet;
    // tracker components
    tracker.Form("TR_%02d", slice);
    converted_tracklet.Form("CONVTR_%02d", slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "");
    AliHLTConfiguration convertedconf(converted_tracklet.Data(), "TPCCATrackerOutputConverter", tracker.Data(), "");    
       
    if (tracklets.Length()>0) tracklets+=" ";
    tracklets+=converted_tracklet;
  }
   
  TString histoInput;
   
  if(bTracklets==kTRUE){
     if(histoInput.Length()>0) 
     histoInput+=" ";
     histoInput+=tracklets;
     
     AliHLTConfiguration tracklethiconf("tracklethi", "TPCTrackHisto", histoInput.Data(), "");
     AliHLTConfiguration trackletrwfconf("trackletrfw", "ROOTFileWriter", "tracklethi", "-datafile TrackletHisto -concatenate-events -overwrite");
  }
  else {
     if(histoInput.Length()>0) 
     histoInput+=" ";
     histoInput+="TPC-clusters";
     histoInput+=" ";
     histoInput+="TPC-globalmerger";
  
     AliHLTConfiguration trhiconf("trhi", "TPCTrackHisto", histoInput.Data(), "");
     AliHLTConfiguration rwfconf("rfw", "ROOTFileWriter", "trhi", "-datafile TrackHisto -concatenate-events -overwrite");
  }


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
 
  if(bTracklets==kTRUE){
     rec.SetOption("HLT", "loglevel=0x7c chains=trackletrfw");
  }
  else {
     rec.SetOption("HLT", "loglevel=0x7c chains=rfw");
  }

  rec.Run();

}

void trackhisto(const char *filename, int minEvent=-1, int maxEvent=-1, Bool_t bTracklets=kFALSE){

  trackhisto(filename, "raw://", minEvent, maxEvent, bTracklets);
}

void trackhisto(){

  cout << "trackhisto: Run HLT TPC tracking and fill histograms" << endl;
  cout << " Usage: aliroot -b -q -l trackhisto.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "   OR " << endl;
  cout << " Usage: aliroot -b -q -l trackhisto.C'(\"file\", \"cdb\", minEvent, maxEvent, kTRUE)'" << endl;
  cout << " if you want to fill a separate histogram file (\TrackletHisto.root)\ with properties of the CA tracker output, no merger called" << endl;
}

