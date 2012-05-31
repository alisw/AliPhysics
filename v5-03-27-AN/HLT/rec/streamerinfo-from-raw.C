// $Id$
/**
 * @file streamerinfo-from-raw.C
 * @brief Extract the streamer info for a raw file and AliRoot version
 *
 * This macro generates the streamer info for all objects in the HLTOUT
 * payload of the specified raw file. The AliRoot version needs to be the
 * same as the one the HLTOUT payload was generated with.
 *
 * The macro has two optional parameters:
 * - filename:  default 'raw://'
 * - OCDB URI:  default 'raw.root'
 *
 * If either file or URI refers to the GRID, the macro connects to alien.
 * Make sure that ROOT is compiled with ALIEN support and that the required
 * alien API libraries can be loaded.
 *
 * You can run this macro with defaults using the following shell command:
 * @code
 *   > aliroot -b -q $ALICE_ROOT/HLT/rec/streamerinfo-from-raw.C
 * @endcode
 * With a raw file from the GRID
 * @code
 *   > aliroot -b -q streamerinfo-from-raw.C'("alien:///alice/data/2009/LHC09d/000102925/raw/09000102925035.10.root")'
 * @endcode
 *
 * Streamer info is stored in a TObjArray wrapped into an AliCDBEntry
 * object, and stored to HLT.StreamerInfo.root. This file can be used
 * directly as OCDB object, it just needs to be renamed.
 *
 * @author Matthias.Richter@ift.uib.no
 */
void streamerinfo_from_raw(const char *filename="raw.root", const char* cdbUri="raw://")
{
  TString tmpStr1=filename;
  TString tmpStr2=cdbUri;
  if ((tmpStr1.Contains("://") && !tmpStr1.BeginsWith("local://")) ||
      !tmpStr2.BeginsWith("local://")) {
    // setup the GRID connection if raw file or OCDB are from GRID
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbUri);

  // Reconstruction settings
  AliReconstruction rec;
  rec.SetRunPlaneEff(kFALSE);
  rec.SetCleanESD(kFALSE);

  rec.SetInput(filename);
  rec.SetRunReconstruction("HLT");
  rec.SetOption("HLT","skip-hltout chains=schemaevo");

  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, "");
  AliHLTConfiguration collector("schemaevo", "ROOTSchemaEvolutionComponent"   , "hltout-publisher", "-file=HLT.StreamerInfo.root");

  // QA options
  rec.SetRunQA(":") ;

  AliLog::Flush();
  rec.Run();
}
