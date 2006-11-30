// @(#) $Id$

#ifndef ALIL3DATAHANDLER_H
#define ALIL3DATAHANDLER_H

#include "AliHLTMemHandler.h"

class AliHLTTransBit;

class AliHLTDataHandler : public AliHLTMemHandler {

 public:
  AliHLTDataHandler();
  ~AliHLTDataHandler();
    
  void Convert10to8Bit();
  Bool_t Memory2CompBinary(UInt_t nrow,AliHLTDigitRowData *data);
  AliHLTDigitRowData *CompBinary2Memory(UInt_t &nrows);

 private:
  
  AliHLTTransBit *fBitTransformer; //! bit transsformer
  
  void Write(Byte_t *comp,UInt_t &index,UShort_t value);
  Short_t Read(Byte_t *comp,UInt_t &index);
  Short_t Test(Byte_t *comp,UInt_t index);
  Bool_t Memory2CompMemory(UInt_t nrow,AliHLTDigitRowData *data,Byte_t *comp);
  UInt_t GetCompMemorySize(UInt_t row,AliHLTDigitRowData *data);
  UInt_t GetMemorySize(UInt_t nrow,Byte_t *comp);
  Bool_t CompMemory2CompBinary(UInt_t nrow,Byte_t *comp,UInt_t size);
  Bool_t CompBinary2CompMemory(UInt_t &nrow,Byte_t *comp);
  UInt_t CompMemory2Memory(UInt_t nrow,AliHLTDigitRowData *data,Byte_t *comp);

  ClassDef(AliHLTDataHandler,1) //Data handler class
};

typedef AliHLTDataHandler AliL3DataHandler; // for backward compatibility

#endif
