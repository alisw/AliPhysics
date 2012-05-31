// $Id$
/*
 * Example macro to run the HLT Conformal mapping tracker embedded into
 * AliRoot reconstruction. The reconstruction is done from the TPC raw
 * data. The last component in the chain is the global merger and the
 * output of the chain is of type {'TRACKS  ':'TPC '}. The output is
 * added to the HLTOUT and processed during the FillESD step.
 * Alternativly the merger is skipped and all {'TRAKSEGS':'TPC '} blocks
 * are forwarded to the output. The merger step is included in the
 * HLTOUT processing in that case (see bUseTrackSegs variable below).
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-hlt-tpc-tracks.C | tee rec-hlt-tpc-tracks.log
 * </pre>
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-hlt-tpc-tracks.C'("input.root")'
 * </pre>
 *
 * In the first section, an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 */
void rec_hlt_tpc_tracks(const char* input="./")
{
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int iMinSlice=0;
  int iMaxSlice=17;
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
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 446");

      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      if (writerInput.Length()>0) writerInput+=" ";
      writerInput+=cf;
    }
    TString tracker;
    // tracker finder components
    tracker.Form("TR_%02d", slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "-pp-run -solenoidBz 0.5");
    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;
  }

  bool bUseTrackSegs=false;

  if (!bUseTrackSegs) {
    // GlobalMerger component
    AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");
  } else {
    // collect all tracksegs blocks
    AliHLTConfiguration mergerconf("globalmerger","BlockFilter",mergerInput.Data(),"");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetFillESD("HLT");
  rec.SetFillTriggerESD(false);
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=globalmerger");
  rec.Run();
}
