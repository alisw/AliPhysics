// $Id$
/**
 * @file raw-hlt-its.C
 * @brief HLT data replay in AliReconstruction
 *
 * <pre>
 * Usage: aliroot -b -q ra-hlt-it.C'("raw.root","cdburi", minevent, maxevent)'
 * Parameter:
 *     input     default 'raw.root'
 *     cdburi    default 'raw://'
 *     minevent  default 0  (no restriction)
 *     maxevent  default -1 (no restriction)
 *     usepredefined  default false
 * </pre>
 *
 * Replay of raw data through the ITS chain. The reduced configuration
 * consists of SPD clusterfinder, SPD vetexer and histogram component.
 * Histograms are written to file 'histos.root'.
 *
 * The predefined chain can be run by specifying the 'GLOBAL-vertexhisto'
 * configuration. In that case also TPC tracking, ITS tracking, generation
 * of ESD and GlobalVertexer are executed as input to the histogramming
 * component.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_its
 */
void rec_its_raw(const char* input="raw.root", 
		 const char* cdburi="raw://",
		 int minevent=0, int maxevent=-1,
		 bool usepredefined=false)
{
  AliReconstruction rec;

  TString uri=cdburi;
  if (!uri.BeginsWith("local://")) {
    TGrid::Connect("alien");
  }
  rec.SetDefaultStorage(cdburi);
  rec.SetRunQA(":") ;

  // AliReconstruction settings
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetLoadAlignFromCDB(0);

  rec.SetRunReconstruction("HLT");

  if (minevent>=0 && maxevent>=minevent) {
    rec.SetEventRange(minevent, maxevent);
  }

  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  if (!usepredefined) {
    // build up the chain only for the SPD CF and vertexer, including the histogramming
    AliHLTConfiguration vertexerSPD("vertexerSPD","ITSVertexerSPD", "ITS-SPD-CF","");
    AliHLTConfiguration vertexhisto("vertexhisto","GlobalVertexerHisto", "vertexerSPD","");

    AliHLTConfiguration filedump("file-dump","ROOTFileWriter", "vertexhisto vertexerSPD", "-datafile histos.root -concatenate-events -overwrite");
  } else {
    // use the default chain from the agents
    AliHLTConfiguration filedump("file-dump","ROOTFileWriter", "GLOBAL-vertexhisto", "-datafile histos.root -concatenate-events -overwrite");
  }

  rec.SetOption("HLT", "loglevel=0x7c ignore-hltout chains=file-dump");

  AliLog::Flush();
  rec.Run();
}
