// @(#) $Id$

#ifndef AliL3_CompressAC
#define AliL3_CompressAC

#include "AliL3Compress.h"
#include "bitio.h"

class AliL3CompressAC : public AliL3Compress {
  
 public:
  AliL3CompressAC();
  AliL3CompressAC(Int_t slice,Int_t patch,Char_t *path="./",Bool_t writeshape=kFALSE,Int_t event=-1);
  virtual ~AliL3CompressAC();
  
  Bool_t CompressFile();
  Bool_t ExpandFile();
  void PrintCompRatio(STDOF *outfile=0);
  void PrintTotals() const;
  
 private:
  UChar_t *fCount;  //! Array of counts
  UInt_t *fTotals;  //! Array of totals, which is actually the model being used during encoding/decoding
  UShort_t fMax; // Max number of counts
  UInt_t fScale; // Scaling factor (not used?)

  UInt_t fRange; // Range (not used?)
  UShort_t fLow; // Low bytes?
  UShort_t fHigh; // High bytes?
  UShort_t fUnderflowBits; // Underflow (not used?)
  UShort_t fCode; // Coded informatio?

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
  
  ClassDef(AliL3CompressAC,1) 

};

#endif
