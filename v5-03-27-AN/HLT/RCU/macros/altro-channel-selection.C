// $Id$
/**
 * Test macro for the AltroChannelSelector component
 *
 * Usage:
 *   aliroot -b -q altro-channel-selection.C | tee altro-channel-selection.C
 *
 * The macro expects simulated TPC raw data in the form TPC_<ddlno>.dll in the
 * current directory, you might need to start the macro in one of the raw<x>
 * folders. You can easily change the sectors and readout partitions below.
 *
 * The simple test writes a fake file with the list of the selected channels.
 *
 * The function has two switch arguments:
 * - directDump=true/false <br>
 *   determines whether the input data should be forwarded directly. Selection
 *   of channels is disabled if \em true
 * - textDump=true/false   <br>
 *   write output either in ascii text dump using the TPCDigitDump or in
 *   binary data
 *
 * Please note that this macro uses also the TPC module, but this does not
 * imply dependencies to the libAliHLTTPC.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_rcu
 */
void altro_channel_selection(bool directDump=false, bool textDump=false)
{
  // this is just a tool to switch the logging systems
  AliHLTLogging log;
  //log.SwitchAliLog(0);

  AliHLTSystem gHLT;
  //gHLT.SetGlobalLoggingLevel(0x7c);

  // load the component library
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTTPC.so");
  gHLT.LoadComponentLibraries("libAliHLTRCU.so");

  // create a dummy pad selection list
  const char* dummySelectionList="/tmp/active-channels.dat";
  FILE* fp = fopen(dummySelectionList, "w");
  if (fp) {
    UShort_t channel=5;
    fwrite(&channel, sizeof(UShort_t), 1, fp);

    UShort_t channel=25;
    fwrite(&channel, sizeof(UShort_t), 1, fp);

    UShort_t channel=56;
    fwrite(&channel, sizeof(UShort_t), 1, fp);

    UShort_t channel=78;
    fwrite(&channel, sizeof(UShort_t), 1, fp);

    UShort_t channel=100;
    fwrite(&channel, sizeof(UShort_t), 1, fp);

    fclose(fp);
  } else {
    cout << "can not open file " << dummySelectionList << " for writing" << endl;
    return;
  }

  // the configuration
  int iMinSlice=0; 
  int iMaxSlice=0;
  int iMinPart=0;
  int iMaxPart=5;
  TString writerInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, selector, activepads;
      TString selectorInput;

      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x%02x%02x%02x%02x -datafile TPC_%d.ddl", slice, slice, part, part, ddlno);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "FilePublisher", NULL , arg.Data());

      // publisher for a dummy active pad list
      activepads.Form("APP_%02d_%d", slice, part);
      arg.Form("-datatype 'HWADDR16' 'TPC ' -dataspec 0x%02x%02x%02x%02x -datafile %s", slice, slice, part, part, dummySelectionList);
      AliHLTConfiguration appconf(activepads.Data(), "FilePublisher", NULL , arg.Data());


      if (selectorInput.Length()>0) selectorInput+=" ";
      selectorInput+=publisher; selectorInput+=" ";
      selectorInput+=activepads;

      // the selector configuration
      selector.Form("CHANNELSELECT_%02d_%d", slice, part);
      AliHLTConfiguration channelselect(selector.Data(), "AltroChannelSelector", selectorInput.Data(), "");

      // add either the raw file directly to output or the filtered one
      if (writerInput.Length()>0) writerInput+=" ";
      if (directDump) {
	writerInput+=publisher;
      } else {
	writerInput+=selector;
      }
    }
  }

  // the writer configuration
  if (textDump)
    AliHLTConfiguration digitdump("digitdump", "TPCDigitDump"   , writerInput.Data(), "-datafile digit.dump -specfmt=_0x%08x -subdir=out_%d -blcknofmt=_0x%x -idfmt=_0x%08x");
  else
    AliHLTConfiguration digitdump("digitdump", "FileWriter"   , writerInput.Data(), "-datafile RAW.ddl -specfmt=_0x%08x -subdir=out_%d -blcknofmt= -idfmt= -skip-datatype");

  // build the ask list and execute
  gHLT.BuildTaskList("digitdump");
  gHLT.Run();

  // delete temporary file
  TString shellcmd;
  shellcmd.Form("rm %s", dummySelectionList);
  gSystem->Exec(shellcmd.Data());
}
