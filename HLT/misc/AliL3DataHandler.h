// @(#) $Id$

#ifndef ALIL3DATAHANDLER_H
#define ALIL3DATAHANDLER_H

#include "AliL3MemHandler.h"

class AliL3TransBit;

class AliL3DataHandler : public AliL3MemHandler {

 public:
  AliL3DataHandler();
  ~AliL3DataHandler();
    
  void Convert10to8Bit();
  Bool_t Memory2CompBinary(UInt_t nrow,AliL3DigitRowData *data);
  AliL3DigitRowData *CompBinary2Memory(UInt_t &nrows);

 private:
  
  AliL3TransBit *fBitTransformer; //! bit transsformer
  
  void Write(Byte_t *comp,UInt_t &index,UShort_t value);
  Short_t Read(Byte_t *comp,UInt_t &index);
  Short_t Test(Byte_t *comp,UInt_t index);
  Bool_t Memory2CompMemory(UInt_t nrow,AliL3DigitRowData *data,Byte_t *comp);
  UInt_t GetCompMemorySize(UInt_t row,AliL3DigitRowData *data);
  UInt_t GetMemorySize(UInt_t nrow,Byte_t *comp);
  Bool_t CompMemory2CompBinary(UInt_t nrow,Byte_t *comp,UInt_t size);
  Bool_t CompBinary2CompMemory(UInt_t &nrow,Byte_t *comp);
  UInt_t CompMemory2Memory(UInt_t nrow,AliL3DigitRowData *data,Byte_t *comp);

  ClassDef(AliL3DataHandler,1) //Data handler class
};

#endif
