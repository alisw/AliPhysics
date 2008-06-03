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
  int iMinSlice=0;
  int iMaxSlice=0;
  int iMinPart=0;
  int iMaxPart=5;

  TString writerInput;
  TString trackerInput;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, cf;

      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // cluster finder components
      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCOfflineClusterer", publisher.Data(), "");

      if (trackerInput.Length()>0) trackerInput+=" ";
      trackerInput+=cf;
      //if (writerInput.Length()>0) writerInput+=" ";
      //writerInput+=cf;
    }
  }

  // one global tracker component
  TString tracker;
  tracker.Form("Global_TR");
  AliHLTConfiguration trackerconf(tracker.Data(), "TPCOfflineTracker", trackerInput.Data(), "");
  if (writerInput.Length()>0) writerInput+=" ";
  writerInput+=tracker;

  // the writer configuration
  AliHLTConfiguration fwconf("sink1", "ROOTFileWriter"   , writerInput.Data(), "-datafile ESD.root -specfmt=_%d -subdir=out_%d -blcknofmt=_0x%x -idfmt=_0x%08x");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput("./");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetFillESD("HLT");
  rec.SetRunQA(kFALSE);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field,kTRUE);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=sink1");
  rec.Run();
}
