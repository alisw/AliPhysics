// $Id$
/**
 * @file publish-rawreader-data.C
 * @brief Publish data from the the RawReader provided by the AliReconstruction.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q 'publish-rawreader-data.C("raw.root", 768, 769)' | tee publish-rawreader-data.log
 * </pre>
 *
 * This macro runs shows a use case of the AliHLTRawReaderPublisherComponent
 * An HLT chain is run inside the AliReconstruction and the DDL payload is
 * extracted from the AliRawReader in a given  equipment range.
 * Input can be either a root file or the path to the directory containing
 * the 'rawx' sub folders.
 *
 * @author Matthias.Richter@ift.uib.no
 */

const char* defaultInput="./";
const int defaultMinId=768;
const int defaultMaxId=769;

void publish_rawreader_data(const char* input, int iMinDDLno, int iMaxDDLno)
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // some defaults
  const char* baseName="RAW.ddl";

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstruction
  if (!input) {
    cerr << "invalid path" << endl;
    cerr << "usage: aliroot -b -q 'publish-rawreader-data.C(\"raw.root\", 768, 769)'" << endl;
    return;
  }

  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetOption("HLT", "libAliHLTUtil.so loglevel=0x7c chains=sink1");

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
  TString writerInput;
  TString arg;

  arg.Form("-minid %d -maxid %d -skipempty -verbose", iMinDDLno, iMaxDDLno);
  AliHLTConfiguration pubconf("publisher", "AliRawReaderPublisher", NULL , arg.Data());
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="publisher";

  // the writer configuration
  arg.Form("-specfmt=_%%d -subdir=out%%d -blocknofmt= -idfmt= -datafile %s", baseName);
  AliHLTConfiguration fwconf("sink1", "FileWriter"   , writerInput.Data(), arg.Data());


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the reconstruction loop
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetFillESD("HLT");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunQA(":");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetLoadAlignFromCDB(0);
  rec.SetFillTriggerESD(kFALSE);
  rec.Run();
}

void publish_rawreader_data()
{
  cout << "runnig from defaults: input " << defaultInput << "  " << defaultMinId << "-" << defaultMaxId << endl;
  publish_rawreader_data(defaultInput, defaultMinId, defaultMaxId);
}
