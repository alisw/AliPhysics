// $Id$
/**
 * @file hltout-collect-esd.C
 * @brief Example for the AliHLTEsdCollectorComponent
 *
 * Example macro to run a small chain with the AliHLTOUTPublisherComponent
 * and the AliHLTEsdCollectorComponent. The AliHLTOUTPublisherComponent
 * is configured to publish all ESD objects from the HLTOUT data, the
 * AliHLTEsdCollectorComponent writes the files using the AliHLTEsdManager.
 *
 * Usage: aliroot -b -q \
 *             hltout-collect-esd.C'("raw.root","local://$ALICE_ROOT/OCDB",0,5)' | tee hltout-collect-esd.log
 * or 
 *  
 *    aliroot -b -l -q \
 *                hltout-collect-esd.C'("alien:///alice/data/2010/LHC10b/000115322/raw/10000115322040.80.root","raw://",0,5)'  | tee hltout-collect-esd.log
 * 
 * 
 * The macro asumes HLTOUT raw data ddl files in order to open an
 * AliRawReaderFile, e.g. simulated using the default AliSimulation with
 * at least SetWriteRawData("HLT") enabled.
 *
 * \b Note: The example disables a few steps in the AliReconstruction,
 * mainly because of crashes in various parts of AliRoot. This does not
 * have any impact to the HLT features to be presented.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void hltout_collect_esd(const char *filename,
                        const char *cdbURI,
                        int minEvent=-1,
                        int maxEvent=-1)
{
  
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=filename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");
    }
  }
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  if (!pHLT) {
    cerr << "fatal error: can not get HLT instance" << endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the configuration chain
  TString arg;

  // the publisher configuration
  arg.Form("-typeid ALIESDV0");
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, arg.Data());

  // the writer configuration
  arg.Form("");
  AliHLTConfiguration collector("sink1", "EsdCollector"   , "hltout-publisher", arg.Data());

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstruction
  AliReconstruction rec; 
  
  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  rec.SetInput(filename);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetFillESD("");
  rec.SetRunQA(":");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetOption("HLT", "libAliHLTUtil.so loglevel=0x7c chains=sink1");
  rec.Run();
}
