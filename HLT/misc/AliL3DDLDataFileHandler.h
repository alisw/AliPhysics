// @(#) $Id$

#ifndef ALIL3DDLDATAFILEHANDLER_H
#define ALIL3DDLDATAFILEHANDLER_H

#ifdef use_newio
#include "../RAW/AliRawReader.h"
#include "../RAW/AliTPCRawStream.h"
#include <TString.h>
#else
#include "AliL3DDLRawReaderFile.h"
#include "AliL3DDLTPCRawStream.h"
#endif
#include "AliL3MemHandler.h"

class AliL3DDLDataFileHandler:public AliL3MemHandler{
 
  private:
#ifdef use_newio
   TString          fFilename;
   Int_t            fEvent;
   AliRawReader    *fReader;
   AliTPCRawStream *fTPCStream;
#else
   AliL3DDLRawReaderFile *fReader;
   AliL3DDLTPCRawStream *fTPCStream;
#endif

  public:

   AliL3DDLDataFileHandler();
   virtual ~AliL3DDLDataFileHandler();

#ifdef use_newio
   Bool_t SetReaderInput(Char_t *name,Int_t event=0);
   Bool_t IsDigit(Int_t i=0);
   AliL3DigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t /*eventmerge*/=kFALSE){return DDLData2Memory(nrow,event);};
#else
   Bool_t SetReaderInput(Char_t *name,Bool_t add=kTRUE);
   Bool_t SetReaderInput(AliL3DDLRawReaderFile *rf);
#endif

   void CloseReaderInput();
   void FreeAll(); //like AliL3MemHandler::Free() or AliL3FileHandler::FreeDigitsTree

   AliL3DigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
   Bool_t DDLData2CompBinary(Int_t event=-1);

   ClassDef(AliL3DDLDataFileHandler,1)   //DDL Data Filehandler class
};
#endif
