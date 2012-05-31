// $Id: rec-hlt-tpc.C 35524 2009-10-14 14:39:46Z kkanaki $
/*
 * Example macro to run the HLT Conformal mapping tracker embedded into
 * AliRoot reconstruction. The reconstruction is done from the TPC raw
 * data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-hlt-tpc.C | tee rec-hlt-tpc.log
 *   aliroot -b -q rec-hlt-tpc.C'("./","decoder All")' | tee rec-hlt-tpc.log
 * </pre>
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-hlt-tpc.C'("input.root")'
 * </pre>
 *
 * The ClusterFinder uses the offline altro decoder v3 by default. The former
 * CF components can be selected by means of the second parameter:
 *    - decoder, uses TPCClusterFinderDecoder. This is default.
 *    - packed, uses TPCClusterFinderPacked
 *
 * Also in the second parameter you can set which output you would like to have:
 *    - ESD, gives you an ESD file. This is default.
 *    - TrackHistogram, will run the TrackHistogram component, and give 
 *      root files with histograms.
 *    - TrackDump, dumps the track struct to a text file.
 *    - ClusterHisto, gives you histograms of the cluster.
 *    - ClusterDump, dumps the cluster struct to text flie.
 *    - All, gives you all 3 output.
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
void rec_hlt_tpc_dEdx(const char* input="./", char* opt="ESD")
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
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Setting up which output to give
  //
  int clusterFinderType=0; // 0 = v3; 1 = decoder; 2 = packed (offline v1)
  bool bUseCA=true;   // use the CA tracker and merger
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=";
  Bool_t esdout=kFALSE, dumpout=kFALSE, histout=kFALSE, chout=kFALSE, cdout=kFALSE;
  TString allArgs=opt;
  TString argument;
  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries(); i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
 
      if (argument.Contains("-statistics=")) {
	argument.ReplaceAll("-statistics=", "");
	if (argument.Contains("root")) {
	  if (option.Length()>0) option+=",";
	  option+="statroot";
	}
	if (argument.Contains("raw")) {
	  if (option.Length()>0) option+=",";
	  option+="statraw";
	}
	continue;
      }
      if (argument.CompareTo("ca", TString::kIgnoreCase)==0) {
	bUseCA=true;
	continue;
      } 
      if (argument.CompareTo("cm", TString::kIgnoreCase)==0) {
	bUseCA=false;
	continue;
      }
      if (argument.CompareTo("decoder",TString::kIgnoreCase)==0) {
	clusterFinderType = 1;
	continue;
      }
      if (argument.CompareTo("packed",TString::kIgnoreCase)==0) {
	clusterFinderType = 2;
	continue;
      }
      if (argument.CompareTo("trackhistogram",TString::kIgnoreCase)==0) {
	histout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="histFile";
	continue;
      }
      if (argument.CompareTo("trackdump",TString::kIgnoreCase)==0) {
	dumpout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="dump";
	continue;
      }
      if (argument.CompareTo("esd",TString::kIgnoreCase)==0) {
	esdout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="esd-converter";
	continue;
      }
      if (argument.CompareTo("esdwrite",TString::kIgnoreCase)==0) {
	esdout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="esdfile";
	continue;
      }
      if (argument.CompareTo("clusterdump",TString::kIgnoreCase)==0) {
	cdout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="cdump";
	continue;
      }      
      if (argument.CompareTo("clusterhisto",TString::kIgnoreCase)==0) {
	chout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="chhisto";
	continue;
      }      
      if (argument.CompareTo("all",TString::kIgnoreCase)==0) {
	histout = kTRUE;
	dumpout = kTRUE;
	esdout = kTRUE;
	chout = kTRUE;
	cdout = kTRUE;
	if (option.Length()>0) option+=",";
	option+="esd-converter,histFile,dump,cdump,chhisto";
	continue;
      }
      else {
	break;
      }
    }
  }
  if (!histout && !dumpout && !esdout && !chout && !cdout) {
    cout << "you need to specify an output option, on or more out of: ESD, TrackHistogram, TrackDump, ClusterHisto, ClusterDump" << endl;
    return;
  }
  if ((option.Contains("statroot") || option.Contains("statraw")) && !esdout) {
    cout << "!!!!!!!! Warning: HLT statistics are only collected for output type ESD !!!!!!!!" << endl;
  }

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
  TString dEdxInput;
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
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      if (clusterFinderType==1) {
	AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 1001");
      } else if (clusterFinderType==2) {
	AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderPacked", publisher.Data(), "-timebins 1001 -sorted");
      } else {
	AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinder32Bit", publisher.Data(), "-timebins 1001");
      }

      if (dEdxInput.Length()>0) dEdxInput+=" ";
      dEdxInput+=cf;
      
      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      if (writerInput.Length()>0) writerInput+=" ";
      writerInput+=cf;
      if (histoInput.Length()>0) histoInput+=" ";
      histoInput+=cf;
      if (cdumpInput.Length()>0) cdumpInput+=" ";
      cdumpInput+=cf;

      if(chout){
	clusterHistoOutput.Form("CH_%02d_%d", slice, part);
	AliHLTConfiguration cfconf(clusterHistoOutput.Data(), "TPCClusterHisto", cf.Data(), "");
	if (histogramHandlerInputClusterFinder.Length()>0) histogramHandlerInputClusterFinder+=" ";
	histogramHandlerInputClusterFinder+=clusterHistoOutput;
      }
    }
    TString tracker;
    // tracker components
    tracker.Form("TR_%02d", slice);
    if (bUseCA) {
      AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "");
    } else {
      AliHLTConfiguration trackerconf(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "-pp-run");
    }

    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=tracker;
    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;
    //add all slice tracks to histo input
    //if (histoInput.Length()>0) histoInput+=" ";
    //histoInput+=tracker;
  }

  // GlobalMerger component
  if (bUseCA) {
    AliHLTConfiguration mergerconf("globalmerger","TPCCAGlobalMerger",mergerInput.Data(),"");
  } else {
    AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");
  }


  if (dEdxInput.Length()>0) dEdxInput+=" ";
  dEdxInput+="globalmerger";

  AliHLTConfiguration dEdX("dEdx", "TPCdEdx", dEdxInput.Data(), "");
  
  //add all global tracks to histo input
  if (histoInput.Length()>0) histoInput+=" ";
  histoInput+="globalmerger";
  
  // specify whether to write all blocks separately or merge the tracks
  // and convert to ESD
  bool writeBlocks=false;

  if(esdout){
    if (writeBlocks) {
      // the writer configuration
      AliHLTConfiguration fwconf("esdfile", "FileWriter"   , writerInput.Data(), "-specfmt=_%d -subdir=out_%d -blcknofmt=_0x%x -idfmt=_0x%08x");
    } else {
           
      // the esd converter configuration
      AliHLTConfiguration esdcconf("esd-converter", "GlobalEsdConverter"   , "globalmerger dEdx", "");
      
      // the root file writer configuration
      AliHLTConfiguration sink("esdfile", "EsdCollector"   , "esd-converter", "-directory hlt-tpc-esd");

      // optional component statistics
      AliHLTConfiguration statroot("statroot", "StatisticsCollector"   , "esd-converter", "-file HLT.statistics.root -publish 0");
      AliHLTConfiguration statraw("statraw", "FileWriter"   , "esd-converter", "-datatype COMPSTAT PRIV -datafile HLT.statistics.raw -concatenate-blocks -concatenate-events");
    }
  }
  //Chain with Track Histogram
  if(histout){
    AliHLTConfiguration histoconf("histo","TPCTrackHisto",histoInput.Data(),"");  
    AliHLTConfiguration fwconf("histFile", "ROOTFileWriter"   , "histo", "-datafile TrackHisto -concatenate-events -overwrite");
  }
  //Chain with Track Dump
  if(dumpout){
    AliHLTConfiguration dumpconf("dump","TPCTrackDump","globalmerger","-directory TrackDump");
  }
  if(chout){
    AliHLTConfiguration cfconf("HHCF", "TPCHistogramHandler", histogramHandlerInputClusterFinder.Data(), "-use-general");
    AliHLTConfiguration rootFileWriterClusters("chhisto", "ROOTFileWriter", "HHCF" , "-datafile histogramHandlerClusterFinder -concatenate-events -overwrite");
  }
  if(cdout){
    AliHLTConfiguration cdumpconf("cdump","TPCClusterDump",cdumpInput.Data(),"-directory ClusterDump");
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");   
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  rec.SetOption("HLT", option);
  rec.Run();
}
