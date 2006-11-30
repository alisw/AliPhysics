// @(#) $Id$

#ifndef ALIL3DDLDATAFILEHANDLER_H
#define ALIL3DDLDATAFILEHANDLER_H

#ifdef use_newio
#include "../RAW/AliRawEvent.h"
#include "../RAW/AliRawReader.h"
#include "../TPC/AliTPCRawStream.h"
#include <TString.h>
#else
#include "AliHLTDDLRawReaderFile.h"
#include "AliHLTDDLTPCRawStream.h"
#endif
#include "AliHLTMemHandler.h"

class AliHLTDDLDataFileHandler:public AliHLTMemHandler{
 
  public:

   AliHLTDDLDataFileHandler();
   virtual ~AliHLTDDLDataFileHandler();

#ifdef use_newio
   Bool_t SetReaderInput(AliRawEvent *rawevent);
   Bool_t SetReaderInput(Char_t *name,Int_t event=0);
   Bool_t IsDigit(Int_t i=0);
   AliHLTDigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t /*eventmerge*/=kFALSE){return DDLData2Memory(nrow,event);};
#else
   Bool_t SetReaderInput(Char_t *name,Bool_t add=kTRUE);
   Bool_t SetReaderInput(AliHLTDDLRawReaderFile *rf);
#endif

   void CloseReaderInput();
   void FreeAll(); //like AliHLTMemHandler::Free() or AliHLTFileHandler::FreeDigitsTree

   AliHLTDigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
   Bool_t DDLData2CompBinary(Int_t event=-1);

   AliTPCRawStream* GetTPCRawStream(){return fTPCStream;}

  private:
#ifdef use_newio
   TString          fFilename; // IO file name
   Int_t            fEvent;    // event number
   AliRawReader    *fReader;   // raw reader
   AliTPCRawStream *fTPCStream;// TPC raw stream
#else
   AliHLTDDLRawReaderFile *fReader; // raw reader
   AliHLTDDLTPCRawStream *fTPCStream; // TPC raw stream
#endif

   ClassDef(AliHLTDDLDataFileHandler,1)   //DDL Data Filehandler class
};

typedef AliHLTDDLDataFileHandler AliL3DDLDataFileHandler; // for backward compatibility

#endif
