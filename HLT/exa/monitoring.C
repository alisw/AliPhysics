// $Id$
/**
 * @file monitoring.C
 * @brief Sample macro for the a monitoring component;
 *
 * Usage:
 * <pre>
 *   aliroot -b -q monitoring.C | tee monitoring.log
 * </pre>
 *
 * This macro illustrates the creation of an HLT monitoring component.
 * The \b Sample-MonitoringComponent component ignores all input data
 * and just fakes some histograms. The histograms can be sent out via
 * different ways. A real component has to have input data and a parent
 * component providing the data in the chain.
 *
 * The ROOTFileWriter component (AliHLTRootFileWriterComponent) provides
 * a simple means to write ROOT objects to a ROOT file.
 *
 * See AliHLTSampleMonitoringComponent for detailed description.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // global initialization of the HLT

  // this is just a tool to switch the logging systems
  AliHLTLogging log;
  //log.SwitchAliLog(0);

  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x7c);

  // load the component library
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTSample.so");

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // now we build up a small chain

  // publisher for the reconfigure event
  TString arg;
  AliHLTConfiguration monitoring("monitoring", "Sample-MonitoringComponent", NULL , "-push-histograms");

  AliHLTConfiguration writer("writer", "ROOTFileWriter", "monitoring" , "");

  // run the chain
  gHLT.BuildTaskList("writer");
  gHLT.Run(1);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // cleanup

}
