// @(#) $Id$

#ifndef ALIL3DDLDATAFILEHANDLER_H
#define ALIL3DDLDATAFILEHANDLER_H

#include "AliL3DDLRawReaderFile.h"
#include "AliL3DDLTPCRawStream.h"
#include "AliL3MemHandler.h"

class AliL3DDLDataFileHandler:public AliL3MemHandler{
 private:
  
  AliL3DDLRawReaderFile *fReader;
  AliL3DDLTPCRawStream *fTPCStream;

 public:
  Bool_t SetReaderInput(Char_t *name,Bool_t add=kTRUE);
  Bool_t SetReaderInput(AliL3DDLRawReaderFile *rf);
  void CloseReaderInput();

  AliL3DDLDataFileHandler();
  virtual ~AliL3DDLDataFileHandler();

  void FreeAll(); //like AliL3MemHandler::Free() or AliL3FileHandler::FreeDigitsTree
  AliL3DigitRowData* DDLData2Memory(UInt_t &nrow,Int_t event=-1);
  Bool_t DDLData2CompBinary(Int_t event=-1);

  ClassDef(AliL3DDLDataFileHandler,1)   //DDL Data Filehandler class
};
#endif
