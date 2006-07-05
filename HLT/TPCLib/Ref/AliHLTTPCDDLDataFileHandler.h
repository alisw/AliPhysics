// @(#) $Id$

#ifndef ALIHLTTPCDDLDATAFILEHANDLER_H
#define ALIHLTTPCDDLDATAFILEHANDLER_H

#ifdef use_newio
#include "AliRawReader.h"
#include "AliTPCRawStream.h"
#include <TString.h>
#else
#include "AliHLTTPCDDLRawReaderFile.h"
#include "AliHLTTPCDDLTPCRawStream.h"
#endif
#include "AliHLTTPCMemHandler.h"

class AliHLTTPCDDLDataFileHandler:public AliHLTTPCMemHandler{
 
  private:
#ifdef use_newio
   TString          fFilename;
   Int_t            fEvent;
   AliRawReader    *fReader;
   AliTPCRawStream *fTPCStream;
#else
   AliHLTTPCDDLRawReaderFile *fReader;
   AliHLTTPCDDLTPCRawStream *fTPCStream;
#endif

  public:

   AliHLTTPCDDLDataFileHandler();
   virtual ~AliHLTTPCDDLDataFileHandler();

#ifdef use_newio
   Bool_t SetReaderInput(Char_t *name,Int_t event=0);
   Bool_t IsDigit(Int_t i=0);
   AliHLTTPCDigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t eventmerge=kFALSE){return DDLData2Memory(nrow,event);};
#else
   Bool_t SetReaderInput(Char_t *name,Bool_t add=kTRUE);
   Bool_t SetReaderInput(AliHLTTPCDDLRawReaderFile *rf);
#endif

   void CloseReaderInput();
   void FreeAll(); //like AliHLTTPCMemHandler::Free() or AliHLTTPCFileHandler::FreeDigitsTree

   AliHLTTPCDigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
   Bool_t DDLData2CompBinary(Int_t event=-1);

   ClassDef(AliHLTTPCDDLDataFileHandler,1)   //DDL Data Filehandler class
};
#endif
