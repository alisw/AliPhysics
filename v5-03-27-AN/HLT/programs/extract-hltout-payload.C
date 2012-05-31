// $Id$
/**
 * @file extract-hltout-payload.C
 * @brief Extraction of data blocks from HLTOUT
 *
 * <pre>
 * Usage: aliroot -b -q extract-hltout-payload.C'("raw.root", "selection", nofEvents)' \
 *            | tee extract-hltout-payload.log
 *
 * Defaults
 *     selection=""      -> no data type selection
 *     nofEvents=-1      -> all events
 * </pre>
 *
 * The raw input file can be accessed directly from the GRID, e.g
 * "alien:///alice/data/2010/LHC10b/000115887/raw/10000115887021.20.root"
 *
 * The macro stores all data blocks from the HLTOUT payload into separated
 * folders for each event. The file names are derived from data type and
 * specification of the block. Data is read from a raw reader. An input
 * file can be specified as the first argument, default is "./" and reads
 * ddl files through AliRawReaderFile.
 *
 * A selection criterion can be specified as second argument, the format
 * of the selection string follows the arguments of the 
 * AliHLTOUTPublisherComponent except from the quote which need to be
 * replaced by brackets due to CINT, e.g.
 *
 * <pre>
 * aliroot -b -q extract-hltout-payload.C'("raw.root", "-datatype {DDL_RAW } ISDD")'
 * aliroot -b -q extract-hltout-payload.C'("raw.root", "-origin {TPC }")'
 * </pre>
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_programs
 */
void extract_hltout_payload(const char* input, const char* selection="", int maxEvent=-1)
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup GRID if input is not a local file
  TString strfile=input;
  if (strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

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
  TString arg;

  // the publisher configuration
  arg.Form("%s", selection);
  arg.ReplaceAll("{", "'");
  arg.ReplaceAll("}", "'");
  AliHLTConfiguration publisher("hltout-publisher", "AliHLTOUTPublisher" , NULL, arg.Data());

  // the writer configuration
  arg.Form("-subdir=event_%%d -blocknofmt= -datafile hltout.dat -specfmt");
  AliHLTConfiguration collector("sink1", "FileWriter"   , "hltout-publisher", arg.Data());

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstruction

  AliHLTReconstructor hltRec;
  hltRec.Init("chains=sink1 ignore-ctp");

  AliRawReader* rawreader=AliRawReader::Create(input);
  if (!rawreader) {
    return;
  }
  rawreader->RewindEvents();
  int count=0;
  if (!rawreader->NextEvent()) {
    cout << "no events found in " << input << endl;
    return;
  }

  do {
    cout << "processing event " << count++ << endl;
    hltRec.Reconstruct(rawreader, NULL);
  } while (rawreader->NextEvent() && (maxEvent<0 || count<maxEvent));
}

void extract_hltout_payload()
{
  cerr << "==============================================================================" << endl;
  cerr << "usage: aliroot -b -q -l extract-hltout-payload.C'(input, selection, maxEvent)'" << endl << endl;
  cerr << "please provide input, e.g. \"raw.root\", or \"./\"" << endl;
  cerr << "optional data type selection, e.g \"-datatype {DDL_RAW } ISDD\", " << endl;
  cerr << "                                  \"-origin {TPC }\", \"-typeid {DDL_RAW }\"" << endl;
  cerr << "optional max event" << endl;
  cerr << "==============================================================================" << endl;
}
