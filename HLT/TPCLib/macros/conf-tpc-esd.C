// $Id$
/*
 * Configuration macro defining a chain running the HLT Conformal mapping
 * tracker embedded into AliRoot simulation.
 * The reconstruction is done from the TPC digits.
 *
 * Several output options are possible:
 * - sink-writeall FileWriter writing all blocks
 * - sink-esd-file write ESD using the TPCEsdWriter
 * - sink-esd      write ESD using the TPCEsdConverter and EsdCollector
 * - sink-clusters relay cluster blocks to HLTOUT
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tpc
 */
{
  // load TPCParam from OCDB
  const char* cdbEntry="TPC/Calib/Parameters";
  AliCDBManager* pMan=AliCDBManager::Instance();
  AliTPCParam* pTPCParam=NULL;
  if (pMan) {
    AliCDBEntry *pEntry = pMan->Get(cdbEntry);
    if (pEntry && 
	pEntry->GetObject() &&
	(pTPCParam=dynamic_cast<AliTPCParam*>(pEntry->GetObject()))) {
    } else {
      HLTWarning("can not load AliTPCParam from CDB entry %s", cdbEntry);
    }
  }

  int iMinSlice=0; 
  int iMaxSlice=35;
  int iMinPart=0;
  int iMaxPart=5;
  TString mergerInput;
  TString writerInput;
  TString sinkClusterInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    TString trackerInput;
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, cf;

      // digit publisher components
      arg.Form("-slice %d -partition %d", slice, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "TPCDigitPublisher", NULL , arg.Data());

      // cluster finder components
      arg="-timebins ";
      if (pTPCParam) arg+=pTPCParam->GetMaxTBin()+1;
      else arg+=446; // old simulated data
      arg+=" -sorted ";
      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderUnpacked", publisher.Data(), arg.Data());
      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      if (sinkClusterInput.Length()>0) sinkClusterInput+=" ";
      sinkClusterInput+=cf;
    }
    TString tracker;
    // tracker finder components
    tracker.Form("TR_%02d", slice);
    AliHLTConfiguration trackerconf(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "-pp-run -solenoidBz 0.5");
    //AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "-solenoidBz 0.5");

    //add all trackers to writer input. Include if you would like all slice tracks written.
    //if (writerInput.Length()>0) writerInput+=" ";
    //writerInput+=tracker;

    // add all clusterfinders to the writer input
    if (writerInput.Length()>0) writerInput+=" ";
    writerInput+=trackerInput;

    if (mergerInput.Length()>0) mergerInput+=" ";
    mergerInput+=tracker;

  }

  // GlobalMerger component
  AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");

  if (writerInput.Length()>0) writerInput+=" ";
  writerInput+="globalmerger";

  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // output section
  //

  //////////////////////////////////////////////////////////////////////////////////////////
  // sink1: id=sink-writeall write all blocks
  AliHLTConfiguration sink1("sink-writeall", "FileWriter"   , writerInput.Data(), "-specfmt -subdir=event_%d -blcknofmt=_0x%x -idfmt=_0x%08x");


  //////////////////////////////////////////////////////////////////////////////////////////
  // sink2: id=sink-esd-file write ESD using the TPCEsdWriter
  AliHLTConfiguration sink2("sink-esd-file", "TPCEsdWriter"   , "globalmerger", "-datafile AliHLTESDs.root");


  //////////////////////////////////////////////////////////////////////////////////////////
  // sink3: id=sink-esd add ESD to HLTOUT using the TPCEsdConverter

  // the esd converter configuration
  //AliHLTConfiguration sink3("sink-esd", "TPCEsdConverter"   , "globalmerger", "");
  AliHLTConfiguration sink3("sink-esd", "TPCEsdConverter"   , mergerInput.Data(), "");


  //////////////////////////////////////////////////////////////////////////////////////////
  // sink4: id=sink-clusters add cluster data blocks to HLTOUT
  AliHLTConfiguration sink4("sink-clusters", "BlockFilter"   , sinkClusterInput.Data(), "-datatype 'CLUSTERS' 'TPC '");
}
