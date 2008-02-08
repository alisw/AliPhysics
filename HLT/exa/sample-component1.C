// $Id:  1.1 2007/11/08 12:17:24 richterm Exp $
/**
 * @file sample-component1.C
 * @brief Sample macro for the component initialization and configuration.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q sample-component1.C | tee sample-component1.log
 * </pre>
 *
 * This macro illustrates the creation of an HLT component and it's
 * initialization and configuration.
 *
 * A component can be initialized by command line arguments. The scan
 * and interpretation of those arguments must be implemented in the
 * DoInit function of the component. The command line arguments are
 * specified in the chain configuration. In this example, it is an
 * AliRoot HLT chain. The same applies for PubSub online HLT chains.
 *
 * Configuration of components is done from configuration objects in
 * the CDB. It's the responsibility of the component to retrieve the
 * CDB entry from the CDB and interprete it. The CDB initialization
 * is done by the framework.
 *
 * The component can also decide if it wants to configure already during
 * initialization (DoInit) from configuration objects.
 *
 * The Sample-component1 (AliHLTSampleComponent1) implements configuration
 * via a string of arguments like e.g. '-config1 config-param -config2'.
 * Two different ways of configuration are implemented:
 * - configuration arguments can be part of the initialization arguments.
 *   All arguments not known to the argument scan in DoInit are treated
 *   as configuration arguments. Scanning of those remaining arguments
 *   is done at the end of the DoInit
 * - if there are no configuration arguments, the configuration is done
 *   from the default object in the CDB
 * - The implemented Reconfigure method retrieves either the object
 *   specified in the reconfiguration event or the default object.
 *
 * The macro defines the CDB in the /tmp folder and creates an object
 * in the CDB. Then it defines a very simple chain, which models the
 * respons of the component to the reconfiguration event.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // some parameters for this test macro

  // the path of a temporary file needed to send the reconfigure event
  const char* tmpFile="/tmp/samplecomponent1-dummy.dat";

  // path of the cdb entry
  const char* cdbEntryPath="HLT/ConfigSample/SampleComponent1";

  // path of the CDB to be created
  const char* cdbUri="local:///tmp/OCDB";

  // initialization arguments for the component
  const char* componentInit="-mandatory1 testarg -mandatory2";

  // configuration arguments for the component
  const char* componentConfig="-config1 config-param -config2";


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // global initialization of the HLT

  // this is just a tool to switch the logging systems
  AliHLTLogging log;
  //log.SwitchAliLog(0);

  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x3c);

  // load the component library
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTSample.so");

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // CDB initialization

  // this is a tool to send the reconfiguration event and the
  // path of the cdb entry
  FILE* fp = fopen(tmpFile, "w");
  if (fp) {
    fprintf(fp, cdbEntryPath);
    fclose(fp);
  }

  // now we create the actual entry in the CDB
  // the CDB is created in /tmp
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet())
  {
    man->SetDefaultStorage(cdbUri);
    man->SetRun(0);
    
    // here is the actual content of the configuration object
    TObjString obj=componentConfig;
    AliCDBPath cdbPath(cdbEntryPath);
    AliCDBId cdbId(cdbPath, 0, 0);
    AliCDBMetaData cdbMetaData;
    man->Put(&obj, cdbId, &cdbMetaData);
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // now we build up a small chain

  // publisher for the reconfigure event
  TString arg;
  arg.Form("-datatype COM_CONF PRIV -datafile %s", tmpFile);
  AliHLTConfiguration reconfevent("reconfevent", "FilePublisher", NULL , arg.Data());

  AliHLTConfiguration sc1("sc1", "Sample-component1", "reconfevent" , componentInit);

  // run the chain
  gHLT.BuildTaskList("sc1");
  gHLT.Run();

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // cleanup

  // delete temporary file
  TString cmd;
  cmd.Form("rm -r %s", tmpFile, cdbUri);
  gSystem->Exec(cmd.Data());
}
