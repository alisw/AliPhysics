// @(#) $Id$

#ifndef ALIL3DDLDATAFILEHANDLER_H
#define ALIL3DDLDATAFILEHANDLER_H

#include "../RAW/AliRawReader.h"
#include "../RAW/AliTPCRawStream.h"
#include "AliL3MemHandler.h"

class AliL3DDLDataFileHandler:public AliL3MemHandler{
 private:
  
  TString          fFilename;
  Int_t            fEvent;
  AliRawReader    *fReader;
  AliTPCRawStream *fTPCStream;

 public:
  Bool_t SetReaderInput(Char_t *name,Int_t event=0);
  void CloseReaderInput();

  AliL3DDLDataFileHandler();
  virtual ~AliL3DDLDataFileHandler();

  void FreeAll(); //like AliL3MemHandler::Free() or AliL3FileHandler::FreeDigitsTree
  Bool_t IsDigit(Int_t i=0);
  AliL3DigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
  Bool_t DDLData2CompBinary(Int_t event=-1);
  AliL3DigitRowData *AliAltroDigits2Memory(UInt_t & nrow,Int_t event=0,Bool_t eventmerge=kFALSE){return DDLData2Memory(nrow,event);};
  ClassDef(AliL3DDLDataFileHandler,1)   //DDL Data Filehandler class
};
#endif
