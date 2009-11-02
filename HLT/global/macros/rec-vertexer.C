/*
 * Example macro to run the HLT vertexer embedded
 * into AliRoot reconstruction. The reconstruction is done from the raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-vertexer.C | tee rec-vertexer.log
 * </pre>
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-vertexer.C'("input.root")'
 * </pre>
 *
 * AliHLTTPCDigitReader32Bit is used to read the data
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
void rec_vertexer(const char* input="./")
{
  
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }

  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
 
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTGlobal.so loglevel=0x7c chains=";
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
  TString histoInput;
  TString histogramHandlerInputClusterFinder;
  TString cdumpInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    TString trackerInput;
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, cf;
      TString clusterHistoOutput;
      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      new AliHLTConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      new AliHLTConfiguration(cf.Data(), "TPCClusterFinder32Bit", publisher.Data(), "-solenoidBz -5");

      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      if (writerInput.Length()>0) writerInput+=" ";
      writerInput+=cf;
      if (histoInput.Length()>0) histoInput+=" ";
      histoInput+=cf;
      if (cdumpInput.Length()>0) cdumpInput+=" ";
      cdumpInput+=cf;
    }
    TString tracker;
    // tracker components
    tracker.Form("TR_%02d", slice);
    new AliHLTConfiguration(tracker.Data(), "TPCCATracker", trackerInput.Data(), "");

    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;
    //add all slice tracks to histo input
    //if (histoInput.Length()>0) histoInput+=" ";
    //histoInput+=tracker;
  }

  // GlobalMerger component
    new AliHLTConfiguration("globalmerger","TPCCAGlobalMerger",mergerInput.Data(),"");
  
    //add all global tracks to histo input
  if (histoInput.Length()>0) histoInput+=" ";
  histoInput+="globalmerger";
  
  // specify whether to write all blocks separately or merge the tracks
  // and convert to ESD
  bool writeBlocks=false;

  // the esd converter configuration
  new AliHLTConfiguration("esd-converter", "GlobalEsdConverter"   , "globalmerger", "-fitTracksToVertex 1");
  
  new AliHLTConfiguration("global-vertexer", "GlobalVertexer"   , "esd-converter", "");
  
  new AliHLTConfiguration("v0HistoOut", "V0Histo"   , "global-vertexer", "");

  new AliHLTConfiguration("GVhistorootfile", "ROOTFileWriter", "global-vertexer" , "-datafile primaryVertexHistograms -concatenate-events -overwrite");
  
  new AliHLTConfiguration("v0historootfile", "ROOTFileWriter", "v0HistoOut" , "-datafile secondaryVertexHistograms -concatenate-events -overwrite");

  option+="v0historootfile";

  option+=",GVhistorootfile";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;


  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // rec.SetDefaultStorage("local://$HOME/HCDB"); 
  //rec.SetSpecificStorage("GRP/GRP/Data",
  //			 Form("local://%s",gSystem->pwd()));

  //  rec.SetSpecificStorage("GRP/GRP/Data","local:///opt/HLT-public/OCDB/LHC09c");
  //  rec.SetSpecificStorage("GRP/CTP/Config","local:///opt/HLT-public/OCDB/LHC09c");
  //  rec.SetSpecificStorage("GRP/CTP/Config", "local:///opt/HLT-public/rec/LHC09c/");

  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");  //add TPC for comparison
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetOption("HLT",option);
  // switch off cleanESD
  rec.SetCleanESD(kFALSE);
  //  rec.SetEventRange(0, 100);
  rec.Run();
}
