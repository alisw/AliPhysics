// @(#) $Id$

#ifndef ALIL3RAWDATAFILEHANDLER_H
#define ALIL3RAWDATAFILEHANDLER_H

#include "AliL3MemHandler.h"

class AliL3RawDataFileHandler:public AliL3MemHandler{
 private:
  
  FILE *fMapping;//!
  ifstream *fInRaw;//!
  ifstream *fInRawPed;//!
  ofstream *fOutRaw;//!

  //found in mapping file
  UInt_t fNChannels; 
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
  Short_t fPedVal;

  Bool_t fConvert;
  Int_t Convert4(Int_t i);     //big2little and vice versa
  Short_t Convert2(Short_t i); //big2little and vice versa

 public:
  AliL3RawDataFileHandler();
  virtual ~AliL3RawDataFileHandler();

  void FreeAll(); //like AliL3MemHandler::Free() or AliL3FileHandler::FreeDigitsTree

  Bool_t SetRawInput(Char_t *name);
  Bool_t SetRawInput(ifstream *file);
  void CloseRawInput(); 
  Int_t ReadRawInput();
  Short_t** GetRawData(Int_t &channels, Int_t & timebins);

  Bool_t SetRawOutput(Char_t *name);
  Bool_t SetRawOutput(ofstream *file);
  void CloseRawOutput(); 
  Int_t StoreRawData(Short_t **charges);

  Bool_t SetMappingFile(Char_t *name);
  Bool_t SetMappingFile(FILE *file);
  void CloseMappingFile(); 
  Int_t ReadMappingFile();
  
  Bool_t SetRawPedestalsInput(Char_t *name);
  Bool_t SetRawPedestalsInput(ifstream *file);
  void CloseRawPedestalsInput(); 
  Int_t ReadRawPedestalsInput();

  void SetBig2Little(Bool_t b){fConvert=b;}
  void SetPedVal(Short_t val){fPedVal=val;}

  Int_t GetRowMinUsed(){return fRowMinUsed;} //smallest row number used in the test
  Int_t GetRowMaxUsed(){return fRowMaxUsed;} //hightest row number used in the test
  Int_t GetPadMinUsed(){return fPadMinUsed;} 
  Int_t GetPadMaxUsed(){return fPadMaxUsed;} 
  Short_t GetPedVal(){return fPedVal;}
  Int_t GetNChannels(){return fNChannels;}
  //Int_t GetNTimeBins(){return fNTimeBins;}

  AliL3DigitRowData* RawData2Memory(UInt_t &nrow,Int_t event=-1);

  ClassDef(AliL3RawDataFileHandler,1)   //RawData Filehandler class
};
#endif

