// $Id$
/**
 * @file extract-ddlraw.C
 * @brief Tool to extract DDL raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q 'extract-ddlraw.C("raw.root", 768, 769)' | tee extract-ddlraw.log
 * </pre>
 *
 * This macro is an example for the AliHLTRawReaderPublisherComponent.
 * It extracts the DDL payload from an AliRawReader in a given equipment
 * range. Input can be either a root file or the path to the directory
 * containing the 'rawx' sub folders.
 *
 * A light-weight AliReconstruction-like setup creates the appropriate
 * RawReader for the specified input and uses the standard AliHLTReconstructor
 * to run a small HLT chain. The chain utilizes the AliRawReaderPublisher
 * to extract the payload of the equipments and writes this to files.
 *
 * @note In this example the AliHLTRawReaderPublisherComponent does not set any data
 * type nor specification for the published data blocks. Please remember
 * to provide appropriate arguments via '-datatype' and '-dataspec'
 * arguments (see AliHLTRawReaderPublisherComponent).
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void extract_ddlraw(const char* input, int iMinDDLno, int iMaxDDLno)
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // some defaults
  const char* baseName="RAW.ddl";

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the RawReader
  if (!input) {
    cerr << "invalid path" << endl;
    cerr << "usage: aliroot -b -q 'extract-ddlraw.C(\"raw.root\", 768, 769)'" << endl;
    return;
  }

  AliRawReader* pRawReader=AliRawReader::Create(input);
  if (!pRawReader) {
    cout << "can not open RawReader for file " << input << endl;
    return;
  }
  if (!pRawReader->NextEvent()) {
    cerr << "no events available" << endl;
    return;
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
  // we show two possible configurations:
  //  1. having one publisher for each ddl, the configurations are created
  //     in a loop
  //  2. all in one publisher
  // can be easily switched with the following
  bool bAllInOne=false;

  TString writerInput;
  TString arg;

  if (!bAllInOne) {
    // create one publisher for each ddl
    for (int ddlno=iMinDDLno; ddlno<=iMaxDDLno; ddlno++) {
      TString arg, publisher;

      // raw data publisher components
      arg.Form("-minid %d -skipempty -verbose", ddlno);
      publisher.Form("DP_%d", ddlno);
      // see AliHLTRawReaderPublisherComponent
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      if (!writerInput.IsNull()) writerInput+=" ";
      writerInput+=publisher;
    }
  } else {
    // publish all ddls by the same component, this is much more
    // effective as it avoids repeated parsing through the data
    arg.Form("-minid %d -maxid %d -skipempty -verbose", iMinDDLno, iMaxDDLno);
    // see AliHLTRawReaderPublisherComponent
    AliHLTConfiguration pubconf("publisher", "AliRawReaderPublisher", NULL , arg.Data());
    if (!writerInput.IsNull()) writerInput+=" ";
    writerInput+="publisher";
  }

  // the writer configuration is the same for both
  arg.Form("-specfmt=_%%d -subdir=raw%%d -blcknofmt= -idfmt= -datafile %s", baseName);
  // see AliHLTFileWriter
  AliHLTConfiguration fwconf("sink1", "FileWriter"   , writerInput.Data(), arg.Data());


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the reconstructor setup
  AliHLTReconstructor hltRec;
  hltRec.SetOption("libAliHLTUtil.so loglevel=0x7c chains=sink1");
  if (hltRec.Init()<0) {
    cerr << "initialization of reconstructor failed" << endl;
    return;
  }

  // this is just a dummy ESD to provide valid parameters to the
  // reconstructor
  AliESDEvent* pESD = new AliESDEvent;
  pESD->CreateStdContent();

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the reconstruction loop
  Int_t event=0;
  UChar_t* pData=NULL;
  pRawReader->RewindEvents();
  while (pRawReader->NextEvent()) {
    cout << "=======================================================" << endl;
    cout << "event " << event << endl;
    cout << "-------------------------------------------------------" << endl;
    pRawReader->Reset();
    hltRec.Reconstruct(pRawReader, NULL);
    event++;
  }

  delete pESD;
}

void extract_ddlraw()
{
  cerr << "===============================================================" << endl;
  cerr << "usage: aliroot -b -q 'extract-ddlraw.C(\"raw.root\", 768, 769)'" << endl << endl;
  cerr << "please provide input, min and max equipment id" << endl;
  cerr << "===============================================================" << endl;
}
