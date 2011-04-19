// $Id$

#if !defined(__CINT__)
#include "AliRawReader.h"
#include "AliAltroRawStreamV3.h"
#include <iostream>
using namespace std;
#endif

/**
 * @file simhlt.C
 * @brief Preliminary macro to run the AliHLTCOMPHuffmanAltro encoder
 */
void run_AliHLTCOMPHuffmanAltro(const char* input, int iMinEquipmentId=768, int iMaxEquipmentId=983, int nofEvents=-1)
{
  AliRawReader* pRawReader=AliRawReader::Create(input);
  if (!pRawReader) {
    cout << "can not open RawReader for file " << input << endl;
    return;
  }

  if (!pRawReader->NextEvent()) {
    cerr << "no events available" << endl;
    return;
  }
  pRawReader->RewindEvents();
  pRawReader->SelectEquipment(0, iMinEquipmentId, iMaxEquipmentId);

  gSystem->Load("libAliHLTComp.so");
  AliHLTCOMPHuffmanAltro encoder(kTRUE, kTRUE, NULL, 0);

  int event=0;
  while (pRawReader->NextEvent() &&
	 (nofEvents<0 || event<nofEvents)) {
    cout << "=======================================================" << endl;
    cout << "event " << event << endl;
    cout << "-------------------------------------------------------" << endl;
    UChar_t* data=NULL;
    while (pRawReader->ReadNextData(data)) {
      encoder.Reset();
      int size=pRawReader->GetEquipmentSize();
      int ddlid=pRawReader->GetEquipmentId();
      TArrayC buffer(size);
      memcpy(buffer.GetArray(), pRawReader->GetDataHeader(), 32);
      memcpy(buffer.GetArray()+32, data, size-32);
      cout << " DDL: " << ddlid << "  size " << size << endl;
      encoder.AddInputData((UChar_t*)buffer.GetArray(), size, ddlid);
      encoder.ProcessData();
      //encoder.CreateCodeTable();
    }
    event++;
  }
}
