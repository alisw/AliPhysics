// $Id$
/**
 * @file compareDigitReaders.C
 * @brief Compare the TPCDigitReaderPacked and TPCDigitReaderDecoder
 *
 * Usage:
 * <pre>
 *   aliroot -b -q 'compareDigitReaders.C("raw.root")' | tee compareDigitReaders.log
 * </pre>
 *
 * Test macro for comparison of the AliHLTTPCDigitReaderPacked and *Decoder.
 * The AliHLTTPCDigitReaderPacked is based on the AliTPCRawStream
 * (AliAltroRawStream). The AliHLTTPCDigitReaderDecoder is based on the
 * fast AliDecoder. The macro loops over all events of the specified data
 * file, instantiates the two readers and compares signal by signal the
 * output of both.
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 */

#ifndef __CINT__
const int sizeofAliRawDataHeader=sizeof(AliRawDataHeader);
#else
// cint does not handle sizeof correctly
const int sizeofAliRawDataHeader=32;
#endif

bool bVerbose=false;

int CompareReaders(AliRawReader* pRawReader);

int compareDigitReaders(const char* input, int iMinEvent=-1, int iMaxEvent=-1)
{
  int iResult=0;
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // some defaults

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the RawReader
  if (!input) {
    cerr << "invalid path" << endl;
    cerr << "usage: aliroot -b -q 'compareDigitReaders.C(\"raw.root\")'" << endl;
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

  pHLT->LoadComponentLibraries("libAliHLTTPC.so");

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the reconstruction loop
  Int_t event=0;
  UChar_t* pData=NULL;
  pRawReader->RewindEvents();
  while (pRawReader->NextEvent() && iResult>=0 && (iMaxEvent<0 || event<=iMaxEvent)) {
    if (event>=iMinEvent) {
      cout << "=======================================================" << endl;
      cout << "event " << event << endl;
      pRawReader->Reset();
      int verifiedDDLs=0;
      while (pRawReader->ReadHeader() && iResult>=0) {
	if ((iResult=CompareReaders(pRawReader))==0) {
	  verifiedDDLs++;
	}
      }
      if (iResult>=0) {
	cout << "event " << event << ": " << verifiedDDLs << " DDL(s) verified" <<endl;
	cout << "-------------------------------------------------------" << endl;
      }
    }
    event++;
  }

  return iResult;
}

int CompareReaders(AliRawReader* pRawReader)
{
  int iResult=0;

  int ddlid=pRawReader->GetEquipmentId();
  int dataSize=pRawReader->GetDataSize();
  if (dataSize<=0 || ddlid<768 || ddlid>983) return 1;

  const AliRawDataHeader* pHeader=pRawReader->GetDataHeader();
  if (pHeader==NULL) {
    cerr << "warning: can not get data header from RawReader, skipping data block ..." << endl;
    return 1;
  }

  TArrayC buffer(dataSize+sizeofAliRawDataHeader);
  UChar_t* pTgt=buffer.GetArray();
  memcpy(pTgt, pHeader, sizeofAliRawDataHeader);
  pTgt+=sizeofAliRawDataHeader;


  if (!pRawReader->ReadNext(pTgt, dataSize)) {
    cerr << "error: reading " << dataSize << " byte(s) from ReadNextData (ddl " << ddlid << ")" << endl;
    return -1;
  }

  int part=0;
  if (ddlid<840) {
    part=ddlid%2; 
  } else {
    part=(ddlid%4)+2;
  }

  if (bVerbose) cout << "data: " << dataSize << " byte(s) on ddl " << ddlid << endl;

  AliHLTTPCDigitReader* pPacked=new AliHLTTPCDigitReaderPacked; 
  AliHLTTPCDigitReader* pDecoder=new AliHLTTPCDigitReaderDecoder;

  pPacked->SetUnsorted(true);
  pTgt=buffer.GetArray();
  if ((iResult=pPacked->InitBlock(pTgt, buffer.GetSize(), part, 0))<0) {
    cerr << "error setting up DigitReaderPacked" << endl;
    delete pPacked;
    return iResult;
  }

  pDecoder->SetUnsorted(true);
  if ((iResult=pDecoder->InitBlock(pTgt, buffer.GetSize(), part, 0))<0) {
    cerr << "error setting up DigitReaderDecoder" << endl;
    delete pPacked;
    delete pDecoder;
    return iResult;
  }

  int iPrintedPart=-1;
  int iPrintedRow=-1;
  int iPrintedPad=-1;
  int iLastTime=-1;

  while (pPacked->Next()) {
    if (!pDecoder->Next()) {
      cerr << endl << "error: no more data on DigitReaderDecoder" << endl;
      return -1;
    }

    if (pPacked->GetAltroBlockHWaddr()!=pDecoder->GetAltroBlockHWaddr()) {
      cerr << endl << "error: hw address mismatch Packed " << pPacked->GetAltroBlockHWaddr() << " Decoder " << pDecoder->GetAltroBlockHWaddr() << endl;
      return -1;
    }
    if (pPacked->GetRow()!=pDecoder->GetRow()) {
      //cerr << endl << "warning: row mismatch Packed " << pPacked->GetRow() << " Decoder " << pDecoder->GetRow() << endl;
      //return -1;
    }
    if (pPacked->GetPad()!=pDecoder->GetPad()) {
      cerr << endl << "error: row mismatch Packed " << pPacked->GetPad() << " Decoder " << pDecoder->GetPad() << endl;
      return -1;
    }
    if (pPacked->GetTime()!=pDecoder->GetTime()) {
      cerr << endl << "error: row mismatch Packed " << pPacked->GetTime() << " Decoder " << pDecoder->GetTime() << endl;
      return -1;
    }

    if (bVerbose) {
      if ((iLastTime!=-1 && iLastTime!=pPacked->GetTime()+1 && iLastTime!=pPacked->GetTime()-1)) {
	cout << "    -> Time: " << iLastTime << endl;
      } else if ((iPrintedPad!=-1 && iPrintedPad!=pPacked->GetPad()) ||
		 (iPrintedRow!=-1 && iPrintedRow!=pPacked->GetRow())) {
	cout << endl;
      }

      if (iPrintedPart!=part) {
	iPrintedPart=part;
	cout << "====================================================================" << endl;
	cout << "    Partition: " << iPrintedPart << endl;
	iPrintedRow=-1;
      }
      if (iPrintedRow!=pPacked->GetRow()) {
	iPrintedRow=pPacked->GetRow();
	cout << "--------------------------------------------------------------------" << endl;
	cout << "Row: " << iPrintedRow << endl;
	iPrintedPad=-1;
      }
      if (iPrintedPad!=pPacked->GetPad()) {
	iPrintedPad=pPacked->GetPad();
	cout << "Row: " << iPrintedRow << "  Pad: " << iPrintedPad << "  HW address: " << pPacked->GetAltroBlockHWaddr() << endl;
	iLastTime=-1;
      }
      if (iLastTime!=pPacked->GetTime()+1 && iLastTime!=pPacked->GetTime()-1 ) {
	cout << "                     Time " << pPacked->GetTime() << ":  ";
      }
      iLastTime=pPacked->GetTime();
      cout << "  " << pPacked->GetSignal();
    }
  }
  if (bVerbose) cout << endl << endl;

  if (pDecoder->Next()) {
    cerr << endl << "error: additional data on DigitReaderDecoder" << endl;
    return -1;
  }

  delete pDecoder;
  delete pPacked;

  return iResult;
}

int compareDigitReaders()
{
  cerr << "===============================================================" << endl;
  cerr << "usage: aliroot -b -q 'compareDigitReaders.C(\"raw.root\")'" << endl << endl;
  cerr << "please provide input, and optional min and max event" << endl;
  cerr << "===============================================================" << endl;
}
