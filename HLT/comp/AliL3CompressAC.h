// @(#) $Id$

#ifndef AliL3_CompressAC
#define AliL3_CompressAC

#include "AliL3Compress.h"
#include "bitio.h"

class AliL3CompressAC : public AliL3Compress {
  
 private:
  UChar_t *fCount;  //!
  UInt_t *fTotals;  //!
  UShort_t fMax;
  UInt_t fScale;

  UInt_t fRange;
  UShort_t fLow;
  UShort_t fHigh;
  UShort_t fUnderflowBits;
  UShort_t fCode;

  void ClearArrays();
  void BuildModel(BIT_FILE *output);
  void RebuildModel(BIT_FILE *input);
  void FillTotals();

  void InitEncoder();
  void InitDecoder(BIT_FILE *input);
  void ConvertIntToSymbol(Int_t value);
  UInt_t ConvertSymbolToInt();
  void EncodeSymbol(BIT_FILE *output);
  void RemoveSymbolFromStream(BIT_FILE *input,Int_t j);  
  void FlushEncoder(BIT_FILE *output);
  
 public:
  AliL3CompressAC();
  AliL3CompressAC(Int_t slice,Int_t patch,Char_t *path="./",Bool_t writeshape=kFALSE,Int_t event=-1);
  virtual ~AliL3CompressAC();
  
  Bool_t CompressFile();
  Bool_t ExpandFile();
  void PrintCompRatio(STDOF *outfile=0);
  void PrintTotals();
  
  ClassDef(AliL3CompressAC,1) 

};

#endif
