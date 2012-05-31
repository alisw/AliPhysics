// $Id$
/**
 * @file split-hltout.C
 * @brief Example for the AliHLTOUTPublisherComponent
 *
 * Example macro to run a small chain with the AliHLTOUTPublisherComponent.
 * The AliHLTOUTPublisherComponent is a tool to publish data blocks from
 * the HLTOUT data stream into an AliHLT reconstruction chain.
 *
 * Usage: aliroot -b -q split-hltout.C | tee split-hltout.log
 *
 * The macro asumes HLTOUT raw data ddl files in order to open an
 * AliRawReaderFile, e.g. simulated by sim-hlt-rawddl.C. A small chain with
 * just a FileWriter connected is run embedded into AliReconstruction.
 *
 * \b Note: The example disables a few steps in the AliReconstruction,
 * mainly because of crashes in various parts of AliRoot. This does not
 * have any impact to the HLT features to be presented.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void split_hltout()
{
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
  arg.Form("");
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, arg.Data());

  // the writer configuration
  arg.Form("-subdir=out%%d -datafile hltout.dat -specfmt=_0x%%x");
  // see AliHLTFileWriter
  AliHLTConfiguration fwconf("sink1", "FileWriter"   , "hltout-publisher", arg.Data());

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
