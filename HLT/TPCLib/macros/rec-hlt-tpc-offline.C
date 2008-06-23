// $Id$
/**
 * @file rec-hlt-tpc-offline.C
 * @brief Test macro for the HLT TPC offline reco wrappers.
 *
 * The macro runs an HLT chain of TPC analysis, using the offline
 * algorithms ans appropriate wrappers.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-hlt-tpc-offline.C | tee rec-hlt-tpc-offline.log
 * </pre>
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The makro asumes that raw data is available in the rawx folders, either
 * simulated or real data.
 *
 * In the first section, an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * Matthias.Richter@ift.uib.no
 */
void rec_hlt_tpc_offline()
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  gSystem->Load("libHLTrec.so");
  AliHLTSystem* gHLT=AliHLTReconstructorBase::GetInstance();

  //gHLT.SwitchAliLog(0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //

  bool sectorClusterer=true; // run clusterer on sector or DDL level
  // check if the AliRawReaderMemory supports multiple buffers
  TClass* info=TClass::GetClass("AliRawReaderMemory");
  TList* methods=info->GetListOfAllPublicMethods();
  if (sectorClusterer && !methods->FindObject("AddBuffer")) {
    cerr << "warning: AliRawReaderMemory does not support multiple buffers, falling back to run clusterer on DDL level" << endl;
    sectorClusterer=false;
  }
 
  int iMinSlice=0;
  int iMaxSlice=17;
  int iMinPart=0;
  int iMaxPart=5;

  int DDLNoFromSlicePatch(int, int);

  TString writerInput;
  TString trackerInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    TString arg, clustererInput;
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString publisher, cf;

      // raw data publisher components
      int ddlno=DDLNoFromSlicePatch(slice, part);
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      if (!sectorClusterer) {
      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCOfflineClusterer", publisher.Data(), "");

      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      //if (writerInput.Length()>0) writerInput+=" ";
      //writerInput+=cf;
      } else {
	if (clustererInput.Length()>0) clustererInput+=" ";
	clustererInput+=publisher;
      }
    }
    if (sectorClusterer) {
      // cluster finder components
      cf.Form("CF_%02d", slice);
      AliHLTConfiguration cfconf(cf.Data(), "TPCOfflineClusterer", clustererInput.Data(), "");

      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
    }
  }

  // one global tracker component
  TString tracker;
  tracker.Form("Global_TR");
  AliHLTConfiguration trackerconf(tracker.Data(), "TPCOfflineTracker", trackerInput.Data(), "");
  if (writerInput.Length()>0) writerInput+=" ";
  writerInput+=tracker;

  // the writer configuration
  AliHLTConfiguration esdwconf("sink1", "EsdCollector"   , writerInput.Data(), "-directory hlt-tpc-offline");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput("./");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetFillESD("");
  rec.SetRunQA(kFALSE);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field,kTRUE);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=sink1");
  rec.Run();
}

int DDLNoFromSlicePatch(int slice, int part)
{
  int ddlno=768;
  if (part>1) ddlno+=72+4*slice+(part-2);
  else ddlno+=2*slice+part;

  return ddlno;
}
