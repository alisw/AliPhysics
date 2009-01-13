// $Id$
/*
 * Example macro to run the HLT global track merger embedded
 * into AliRoot reconstruction. The reconstruction is done from the raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-hlt-global-merger.C | tee rec-hlt-global-merger.log
 * </pre>
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-hlt-global-merger.C'("input.root")'
 * </pre>
 *
 * By the second parameter the digit reader can be chosen, default is
 * AliHLTTPCDigitReaderPacked (=false). Set to true to use
 * AliHLTTPCDigitReaderDecoder
 *
 * In the first section, an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_global
 * @author 
 */
void rec_hlt_global_merger(const char* input="./", bool bUseClusterFinderDecoder=true)
{
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  gSystem->Exec("rm galice.root");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  gSystem->Load("libHLTrec.so");
  AliHLTSystem* gHLT=AliHLTReconstructorBase::GetInstance();
  gHLT->SetGlobalLoggingLevel(0x4); 
  /*  enum AliHLTComponentLogSeverity {       
      kHLTLogNone      = 0,
      kHLTLogBenchmark = 0x1,
      kHLTLogDebug     = 0x2,
      kHLTLogInfo      = 0x4,
      kHLTLogWarning   = 0x8,
      kHLTLogError     = 0x10,
      kHLTLogFatal     = 0x20,      
      few important messages not to be filtered out.
      redirected to kHLTLogInfo in AliRoot
      kHLTLogImportant = 0x40,
      special value to enable all messages 
      kHLTLogAll       = 0x7f,
      the default logging filter 
      kHLTLogDefault   = 0x79
      useful           = 0x45
  */

  // set TPC debug stream level
  AliTPCReconstructor::SetStreamLevel(1);

  //gHLT->LoadComponentLibraries("libAliHLTUtil.so libAliHLTRCU.so libAliHLTTRD.so libAliHLTTPC.so libAliHLTGlobal.so");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //

  //  Global merger
  TString globalmergerInput;

  // TPC
  int iMinSlice=0;
  int iMaxSlice=35;
  //int iMaxSlice=17;
  int iMinPart=0;
  int iMaxPart=5;

  TString writerInput;
  TString mergerInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    TString trackerInput;
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, cf;

      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d",slice,part);
      if (bUseClusterFinderDecoder) {
	AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 446");
      } else {
	AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderPacked", publisher.Data(), "-timebins 446 -sorted");
      }
      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      if (writerInput.Length()>0) writerInput+=" ";
      writerInput+=cf;
    }
    TString tracker;
    // tracker finder components
    tracker.Form("TR_%02d",slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "-solenoidBz 5");
    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;
  }
  // TPC GlobalMerger component
  AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");

  //
  // TRD
  //
  int iMinSlice=0;
  int iMaxSlice=17;
  //int iMaxSlice=8;

  TString writerInput1;
  TString mergerInput1;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    //TString trackerInput1;
    TString arg1, publisher1, cf1;

    // raw data publisher components
    int ddlno=1024;
    ddlno += slice;

    arg1.Form("-minid %d -datatype 'DDL_RAW ' 'TRD ' -verbose", ddlno);
    publisher1.Form("DP1_%02d", slice);
    AliHLTConfiguration pubconf(publisher1.Data(), "AliRawReaderPublisher", NULL , arg1.Data());

    //arg1.Form("-datatype 'DDL_RAW ' 'TRD ' -datafile raw0/TRD_%d.ddl", ddlno);
    //publisher1.Form("DP1_%02d",slice);
    //AliHLTConfiguration pubconf(publisher1.Data(), "FilePublisher", NULL , arg.Data());
    //cout << arg.Data() << endl; 

    // cluster finder components
    cf1.Form("CF1_%02d",slice);
    AliHLTConfiguration cfconf(cf1.Data(), "TRDClusterizer", publisher1.Data(), "output_percentage 1000 -geometry geometry.root -lowflux -simulation");

    // tracker finder components
    TString tracker1;
    tracker1.Form("TR1_%02d",slice);
    AliHLTConfiguration trackerconf(tracker1.Data(), "TRDTrackerV1", cf1.Data(), "output_percentage 300 -NTimeBins 24 -magnetic_field_ON -lowflux -geometry geometry.root");

    if (writerInput1.Length()>0) writerInput1+=" ";
    writerInput1+=tracker1;

    if (mergerInput1.Length()>0) mergerInput1+=" ";
    mergerInput1 += tracker1;
  }

  if (globalmergerInput.Length()>0) globalmergerInput+=" ";
  globalmergerInput += "globalmerger";

  if (globalmergerInput.Length()>0) globalmergerInput+=" ";
  globalmergerInput += mergerInput1;

  //TString datatype = "sim";
  //AliHLTConfiguration HWriterTR( "HWriterTR", "FileWriter", writerInput.Data(), "-directory . -datafile tr_" + datatype + ".out");
  
  //GlobalTrackMerger component
  AliHLTConfiguration gtrackmergerconf("gtrackmerger","GlobalTrackMerger",globalmergerInput.Data(),"");

  // the root file writer configuration
  AliHLTConfiguration sink("sink1", "EsdCollector" , "gtrackmerger" , "-directory hlt-global-track");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field,kTRUE);
  rec.SetFillESD("HLT");
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTRD.so libAliHLTTPC.so libAliHLTGlobal.so chains=sink1");
  rec.Run();
}
