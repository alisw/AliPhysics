// $Id$
// 
// @file run-compression.C
// @brief Define and run custom chains for the TPC data compression
// @author Matthias.Richter@ift.uib.no
//
// The macro can be used to either define chains to be run in the
// AliRoot reconstruction or to run a chain standalone. Some of the
// configurations can only be run embedded into AliRoot reconstruction
// in order to couple to the input.
// 1) AliRoot reconstruction examples, just use the macro in front of
//    the reconstruction macro, e.g.
// aliroot -b -q -l runmonitor.C(2) $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root", "local://OCDB", 0, 5, "HLT", "loglevel=0x7c chains=compressor-input-writer")'
// 
// 2) standalone chains can be used if all input files are available
// on disk and the configuration file for the publisher is existing
// aliroot -b -q -l runmonitor.C'(2, 2, 1, "writer", "compressor compressor-publisher", 167808, "local://OCDB")'
//
// Chains:
// 'hltout-cluster-writer' writes clusters for every event in directory
//              hltout-compressed-cluster and a publisher configuration
//              file hltout-compressed-cluster.txt 
//
// 'tpc-raw-writer' writes TPC raw DDL files per event in directory tpc-raw
//              and a publisher configuration file tpc-raw.txt
//
// 'compressor-input-writer' writes input to compression component (clusters
//              from HWCF and reconstructed tracks) to directory compressor-input
//              and publisher configuration file compressor-input.txt
//
// 'compressor' the compressor component is configured to write a
//              statistics file HLT.TPCDataCompression-histograms*.root
//
// 'huffmantrainer' huffman table trainer
//
// 'writer'     cluster monitor component, requires to define input, e.g
//              "compressor compressor-publisher"
//              "compressed-cluster-publisher"
//              "monitor-publisher"
//              "hltout-publisher"
//              "TPC-hwcfdata TPC-compression" (note: defined by AliHLTTPCAgent)
const int defaultMode=0;
const int defaultDeflaterMode=2;
const char* defaultMonitorInput="compressor compressor-publisher";
const char* defaultCDBUri="local://OCDB";
void run_compression(int mode=defaultMode, int deflaterMode=defaultDeflaterMode, int events=1,
		const char* chain=NULL,
		const char* monitorInput=defaultMonitorInput,
		int runno=-1,
		const char* cdbURI=defaultCDBUri)
{
  // setup the OCDB access
  // required to load the GRP entry in order to initialize the magnetic field
  if (runno>=0) {
    AliCDBManager::Instance()->SetDefaultStorage(cdbURI);
    AliCDBManager::Instance()->SetRun(runno);
    AliGRPManager grpman;
    grpman.ReadGRPEntry();
    grpman.SetMagField();
  }

  // init the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // list of configurations
  //
  ///////////////////////////////////////////////////////////////////////////////////////////

  // handling of compressed clusters from HLTOUT
  // 'hltout-cluster-writer' writes clusters for every event in directory hltout-compressed-cluster
  // and a publisher configuration file hltout-compressed-cluster.txt 
  // 'compressed-cluster-publisher' publishes according to list hltout-compressed-cluster.txt
  AliHLTConfiguration hltoutpublisher("hltout-publisher", "AliHLTOUTPublisher", "", "-datatype 'REMCLSCM' 'TPC '");
  AliHLTConfiguration hltoutclusterwriter("hltout-cluster-writer", "FileWriter", "hltout-publisher", "-directory hltout-compressed-cluster -subdir -specfmt=_0x%08x -blocknofmt= -publisher-conf hltout-compressed-cluster.txt");
  AliHLTConfiguration hltoutclusterpublisher("compressed-cluster-publisher", "FilePublisher", "", "-datafilelist hltout-compressed-cluster.txt");

  // processing of TPC raw data
  // 'tpc-raw-writer' writes TPC raw DDL files per event in directory tpc-raw and a publisher
  // configuration file tpc-raw.txt
  // 'compressor-input-writer' writes input to compression component (clusters from HWCF and reconstructed tracks)
  // to directory compressor-input and publisher configuration file compressor-input.txt
  // 'compressor-publisher' publishes files according to configuration file
  AliHLTConfiguration tpcrawwriter("tpc-raw-writer", "FileWriter","TPC-raw-data", "-directory tpc-raw -subdir -specfmt=_0x%08x -blocknofmt= -publisher-conf tpc-raw.txt");
  AliHLTConfiguration emulatorhwclust1writer("compressor-input-writer", "FileWriter","TPC-hwcfdata TPC-globalmerger", "-directory compressor-input -subdir -specfmt=_0x%08x -blocknofmt= -publisher-conf compressor-input.txt");
  AliHLTConfiguration compressorpublisher("compressor-publisher", "FilePublisher", "", "-datafilelist compressor-input.txt");

  // compressor configuration
  // input from file list created by 'compressor-input-writer'
  TString compressorArgument;
  if (mode>0) compressorArgument+=Form(" -mode %d", mode); // take default from configuration object if mode==0
  if (deflaterMode>0)   compressorArgument+=Form(" -deflater-mode %d", deflaterMode);  // take default from configuration object if deflaterMode==0
  compressorArgument+=Form(" -histogram-file HLT.TPCDataCompression-histograms-mode%d-%s.root -cluster-verification 0", mode, (deflaterMode==2?"huffman":"simple"));
  AliHLTConfiguration compressor("compressor", "TPCDataCompressor", "compressor-publisher",  compressorArgument.Data());

  // huffman trainer configuration
  // input from file list created by 'compressor-input-writer'
  TString trainerArgument;
  trainerArgument.Form("-deflater-mode 3 -mode %d", mode);
  AliHLTConfiguration trainer("huffmantrainer", "TPCDataCompressor", "compressor-publisher", trainerArgument.Data());

  // writer component for the compressor output
  AliHLTConfiguration compressoroutputwriter("compressor-output-writer", "FileWriter","compressor", "-directory compressor-output -subdir -specfmt=_0x%08x -blocknofmt= -publisher-conf compressor-output.txt");
  AliHLTConfiguration compressoroutputpublisher("compressor-data", "FilePublisher","", "-datafilelist compressor-output.txt");

  // specifc configuration to publish some data blocks for the
  // monitoring component
  AliHLTConfiguration monitorpublisher("monitor-publisher", "FilePublisher", "", "-datafilelist monitor-input.txt");

  // the monitoring ccomponent setup
  TString writerArguments(Form("-concatenate-events -overwrite -datafile HLT.TPCcluster-compression-mode%d-histograms.root",mode));
  AliHLTConfiguration monitor("monitor", "TPCDataCompressorMonitor", monitorInput, "");
  AliHLTConfiguration writer("writer", "ROOTFileWriter", "monitor", writerArguments);

  if (chain) {
    pHLT->ScanOptions("loglevel=0x7c");
    pHLT->BuildTaskList(chain);
    pHLT->Run(events);
  }
}

// an abbreviated version setting default deflaterMode to 2
void run_compression(int mode, int events,
		const char* chain,
		const char* monitorInput=defaultMonitorInput,
		int runno=-1,
		const char* cdbURI=defaultCDBUri)
{
  run_compression(mode, defaultDeflaterMode, events, chain, monitorInput, runno, cdbURI);
}

// an abbreviated version setting default deflaterMode to 2 and
// omitting monitorInput
// example:
//  aliroot -b -q -l runmonitor.C'(1, 5, "compressor", "", 167808)'
void run_compression(int mode, int events,
		const char* chain,
		int runno,
		const char* cdbURI=defaultCDBUri)
{
  run_compression(mode, defaultDeflaterMode, events, chain, defaultMonitorInput, runno, cdbURI);
}
