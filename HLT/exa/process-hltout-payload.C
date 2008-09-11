// $Id$
/**
 * @file process-hltout-payload.C
 * @brief Standard processing of HLTOUT payload
 *
 * <pre>
 * Usage: aliroot -b -q process-hltout-payload.C'("raw.root")' | tee process-hltout-payload.log
 * </pre>
 *
 * Open a raw reader for the specified input and process HLTOUT payload
 * using the AliHLTReconstructor.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void process_hltout_payload(const char* input)
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstructor

  gSystem->Load("libHLTrec");
  AliHLTReconstructor hltRec;
  hltRec.Init("chains=");

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

  hltRec.ProcessHLTOUT(rawreader, NULL);
}

void process_hltout_payload()
{
  cerr << "===============================================================" << endl;
  cerr << "usage: aliroot -b -q -l process-hltout-payload.C'(\"raw.root\")'" << endl << endl;
  cerr << "please provide input, e.g. \"raw.root\", or \"./\"" << endl;
  cerr << "===============================================================" << endl;
}
