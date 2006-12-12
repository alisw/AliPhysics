// @(#) $Id$

#ifndef ALIHLTTPCDDLDATAFILEHANDLER_H
#define ALIHLTTPCDDLDATAFILEHANDLER_H

#include "AliRawEvent.h"
#include "AliRawReader.h"
#include "AliTPCRawStream.h"
#include <TString.h>
#include "AliHLTTPCMemHandler.h"

class AliHLTTPCDDLDataFileHandler:public AliHLTTPCMemHandler{
 
  public:

   AliHLTTPCDDLDataFileHandler();
   virtual ~AliHLTTPCDDLDataFileHandler();

   Bool_t SetReaderInput(AliRawEvent *rawevent);
   Bool_t SetReaderInput(Char_t *name,Int_t event=0);
   Bool_t IsDigit(Int_t i=0) const;
   AliHLTTPCDigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t /*eventmerge*/=kFALSE){return DDLData2Memory(nrow,event);};

   void CloseReaderInput();
   void FreeAll(); //like AliHLTTPCMemHandler::Free() or AliHLTTPCFileHandler::FreeDigitsTree

   AliHLTTPCDigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
   Bool_t DDLData2CompBinary(Int_t event=-1);

   AliTPCRawStream* GetTPCRawStream(){return fTPCStream;}

  private:
   TString          fFilename; // IO file name
   Int_t            fEvent;    // event number
   AliRawReader    *fReader;   // raw reader
   AliTPCRawStream *fTPCStream;// TPC raw stream

   ClassDef(AliHLTTPCDDLDataFileHandler,1)   //DDL Data Filehandler class
};
#endif
