// $Id$
/*
 * Example macro to run the HLT TPC Cellular Automaton tracker embedded
 * into AliRoot reconstruction. The reconstruction is done from the TPC
 * raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-hlt-tpc-ca.C | tee rec-hlt-tpc-ca.log
 * </pre>
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-hlt-tpc-ca.C'("input.root")'
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
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 */
void rec_hlt_tpc_ca(const char* input="./", bool bUseClusterFinderDecoder=true)
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
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int iMinSlice=0;
  int iMaxSlice=35;
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
      cf.Form("CF_%02d_%d", slice, part);
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
    tracker.Form("TR_%02d", slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "-solenoidBz 5");
    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;
  }

  // GlobalMerger component
  AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");

  // the esd converter configuration
  AliHLTConfiguration esdcconf("esd-converter", "TPCEsdConverter"   , "globalmerger", "-tree");
  
  // the root file writer configuration
  AliHLTConfiguration sink("sink1", "EsdCollector"   , "esd-converter", "-directory hlt-tpc-ca");

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
  rec.SetFillESD("HLT");
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so chains=sink1");
  rec.Run();
}
