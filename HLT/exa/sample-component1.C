// $Id$
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
 * initialization and configuration, including update of values from
 * DCS.
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
  const char* tmpFile1="/tmp/samplecomponent1-comconf.dat";
  const char* tmpFile2="/tmp/samplecomponent1-updtdcs.dat";

  // path of the cdb entry
  const char* cdbEntryPath="HLT/ConfigSample/SampleComponent1";

  // path of detectors with 'active' preprocessors
  const char* prepDetectors="TPC PHOS";

  // path of the CDB to be created
  const char* cdbLocation="/tmp/OCDB";
  TString cdbUri; cdbUri.Form("local://%s", cdbLocation);

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
  gHLT.SetGlobalLoggingLevel(0x7c);

  // load the component library
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTSample.so");

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // CDB initialization

  // this is a tool to send the reconfiguration event and the
  // path of the cdb entry
  FILE* fp = fopen(tmpFile1, "w");
  if (fp) {
    fprintf(fp, cdbEntryPath);
    fclose(fp);
  }

  // this is a tool to send the update DCS event and the
  // path of the cdb entry
  FILE* fp = fopen(tmpFile2, "w");
  if (fp) {
    fprintf(fp, prepDetectors);
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
  arg.Form("-datatype COM_CONF PRIV -datafile %s -nextevent "
	   "-datatype UPDT_DCS PRIV -datafile %s", tmpFile1, tmpFile2);
  AliHLTConfiguration sep("steeringevents", "FilePublisher", NULL , arg.Data());

  AliHLTConfiguration sc1("sc1", "Sample-component1", "steeringevents" , componentInit);

  // build the chain
  gHLT.BuildTaskList("sc1");

  // run two events, in the 1st event the component reconfiguration is emulated
  // in the 2nd one the update of Preprocessor values
  gHLT.Run(2);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  // cleanup

  // delete temporary file
  TString cmd;
  cmd.Form("rm -r %s %s %s", tmpFile1, tmpFile2, cdbLocation);
  gSystem->Exec(cmd.Data());
}
