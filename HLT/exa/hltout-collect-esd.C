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
 * Usage: aliroot -b -q hltout-collect-esd.C | tee hltout-collect-esd.log
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
void hltout_collect_esd()
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase!::GetInstance();
  if (!pHLT) {
    cerr << "fatal error: can not get HLT instance" << endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the configuration chain
  TString arg;

  // the publisher configuration
  arg.Form("-typeid ESD_TREE -typeid ALIESDV0");
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, arg.Data());

  // the writer configuration
  arg.Form("");
  AliHLTConfiguration collector("sink1", "EsdCollector"   , "hltout-publisher", arg.Data());

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstruction
  AliReconstruction rec;
  rec.SetInput("./");
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetFillESD("");
  rec.SetRunQA(":");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetOption("HLT", "libAliHLTUtil.so loglevel=0x7c chains=sink1");
  rec.Run();
}
