// $Id$
/**
 * @file rec-hlt-globalhisto.C
 * @brief Run reconstruction and fill the global histograms
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     rec-hlt-globalhisto.C'("file", "cdb", minEvent, maxEvent)'
 * Example:
 *   aliroot -b -q -l \
 *     rec-hlt-globalhisto.C'("raw.root", "local://$ALICE_ROOT/OCDB")'
 * </pre>
 * 
 * Macro runs HLT reconstruction on raw data, all other AliRoot
 * modules switched off. The output of the GlobalEsdConverter is filled
 * into the GlobalHisto component, it's output is attached to a Root file
 * writer producing the file HLT.Histo.root.
 *
 * A raw file, either simulated or real, is needed and the corresponding
 * OCDB has to be specified.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_qa
 */
void rec_hlt_globalhisto(const char *filename,
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

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  AliHLTConfiguration globalhisto("globalhisto", "GlobalHisto", "GLOBAL-esd-converter",
				  "-histogram TrackPt -size 1000 -expression Track_pt "
				  "-histogram TrackPhi -size 1000 -expression Track_phi "
				  "-histogram TrackEta -size 1000 -expression Track_eta "
				  "-histogram TrackCount -size 1000 -expression trackcount "
				  );
  AliHLTConfiguration writer("writer", "ROOTFileWriter", "globalhisto", "-datafile HLT.Histo.root -overwrite -concatenate-events");

  TString hltOptions="loglevel=0x7c ignore-hltout chains=writer";

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  rec.SetRunReconstruction("HLT");

  // QA options
  TString qaOptions="HLT";
  qaOptions+=":ALL";
  rec.SetRunQA(qaOptions) ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunMultFinder(kFALSE);
  rec.SetInput(filename);
  rec.SetOption("HLT", hltOptions);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}
