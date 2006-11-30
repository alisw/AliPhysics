// @(#) $Id$

#ifndef ALIL3RAWDATAFILEHANDLER_H
#define ALIL3RAWDATAFILEHANDLER_H

#include "AliHLTMemHandler.h"

class AliHLTRawDataFileHandler:public AliHLTMemHandler {
 private:
  
  FILE *fMapping;//!
  ifstream *fInRaw;//!
  ifstream *fInRawPed;//!
  ofstream *fOutRaw;//!

  //found in mapping file
  UInt_t fNChannels;   //number of channels 
  Byte_t *fRow;//!     //store the channeln to row mapping
  Byte_t *fPad;//!     //store the channel to pad mapping
  Short_t **fRowPad;//! //store the row-and-pad to channel mapping

  Int_t fRowMinUsed; //min row val (found from mapping)
  Int_t fRowMaxUsed; //max row val (found from mapping)
  Int_t fPadMinUsed; //min pad val (found from mapping)
  Int_t fPadMaxUsed; //max pad val (found from mapping)
  Int_t fNPads[159]; //pads used (found from mapping)
  Int_t fNTimeBins;  //stored in data header

  Short_t **fPedestals;//! pedestal values, if not used, fPedVal is used
  Short_t fPedVal; //ped val if not used per channel

  Short_t **fCharges;//! charge values read from pointer or from file

  Bool_t fConvert; //convert big/little
  Int_t Convert4(Int_t i) const;     //big2little and vice versa
  Short_t Convert2(Short_t i) const; //big2little and vice versa

 public:
  AliHLTRawDataFileHandler();
  virtual ~AliHLTRawDataFileHandler();

  void FreeAll(); //like AliHLTMemHandler::Free() or AliHLTFileHandler::FreeDigitsTree

  Bool_t SetRawInput(Char_t *name);
  Bool_t SetRawInput(STDIF *file);
  void CloseRawInput(); 
  Int_t ReadRawInput();

  Int_t ReadRawInputPointer(const Char_t *ptr);

  Short_t** GetRawData(Int_t &channels, Int_t & timebins);

  Bool_t SetRawOutput(Char_t *name);
  Bool_t SetRawOutput(STDOF *file);
  void CloseRawOutput(); 
  Int_t StoreRawData(Short_t **charges);

  Bool_t SetMappingFile(Char_t *name);
  Bool_t SetMappingFile(FILE *file);
  void CloseMappingFile(); 
  Int_t ReadMappingFile();
  
  Bool_t SetRawPedestalsInput(Char_t *name);
  Bool_t SetRawPedestalsInput(STDIF *file);
  void CloseRawPedestalsInput(); 
  Int_t ReadRawPedestalsInput();

  void SetBig2Little(Bool_t b){fConvert=b;}
  void SetPedVal(Short_t val){fPedVal=val;}

  Int_t GetRowMinUsed() const {return fRowMinUsed;} //smallest row number used in the test
  Int_t GetRowMaxUsed() const {return fRowMaxUsed;} //hightest row number used in the test
  Int_t GetPadMinUsed() const {return fPadMinUsed;} 
  Int_t GetPadMaxUsed() const {return fPadMaxUsed;} 
  Short_t GetPedVal()   const {return fPedVal;}
  Int_t GetNChannels()  const {return fNChannels;}

  AliHLTDigitRowData* RawData2Memory(UInt_t &nrow,Int_t event=-1);
  Bool_t RawData2CompBinary(Int_t event=-1);

  ClassDef(AliHLTRawDataFileHandler,1)   //RawData Filehandler class
};

typedef AliHLTRawDataFileHandler AliL3RawDataFileHandler; // for backward compatibility

#endif

