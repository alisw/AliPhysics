/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////////////////////
// Class used for generating the files containung data, required by the Data Challenge //
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef AliITSDDLRAWDATA_H
#define AliITSDDLRAWDATA_H

class AliITS;
class TTree;

class AliITSDDLRawData:public TObject{
 public:
  AliITSDDLRawData();//default constructor
  virtual ~AliITSDDLRawData(){;}//destructor
  AliITSDDLRawData(const AliITSDDLRawData &source); // copy constructor
  AliITSDDLRawData& operator=(const AliITSDDLRawData &source); // ass. op.
  // This method generates the files with the Silicon pixel detector data
  Int_t RawDataSPD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber=2);
  // This method generates the files with the Silicon drift detector data
  Int_t RawDataSDD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber=4);
  // This method generates the files with the Silicon pixel detector data
  Int_t RawDataSSD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber=2);
  //A debugging method used to test the files generated for the SPD.
  void  TestFormat();
 private: 
  //This method formats and stores in buf all the digits of a SPD module
  void  GetDigitsSPD(TClonesArray *ITSdigits, Int_t mod, ULong_t *buf);
  //This method formats and stores in buf all the digits of a SDD module
  void  GetDigitsSDD(TClonesArray *ITSdigits, Int_t mod, ULong_t *buf);
  //This method formats and stores in buf all the digits of a SSD module
  void  GetDigitsSSD(TClonesArray *ITSdigits, Int_t mod, ULong_t *buf);
  //This method stores the value of the variable Word of StopBit-StartBit+1 bits 
  //in BaseWord, starting from the bit StartBit
  void  PackWord(ULong_t &BaseWord, ULong_t Word, Int_t StartBit, Int_t StopBit);
  //This method extracts a group of adiacents bits, specified by StartBit and StopBit, 
  //from the word PackedWord. The resulting word is saved in the Word variable
  void  UnpackWord(ULong_t PackedWord, Int_t StartBit, Int_t StopBit, ULong_t &Word);
  //The following three methods are used to store the data according to the 
  //Silicon pixel detector data format
  void  WriteChipHeader(Int_t ChipAddr,Int_t EventCnt,ULong_t &BaseWord);
  void  WriteChipTrailer(ULong_t *buf,Int_t ChipHitCount,ULong_t &BaseWord);
  void  WriteHit(ULong_t *buf,Int_t RowAddr,Int_t HitAddr,ULong_t &BaseWord);
  //Methods used for reading and dubugging SPD data files
  void  ReadChipHeader(Int_t &ChipAddr,Int_t &EventCnt,ULong_t BaseWord);
  void  ReadChipTrailer(Int_t &ChipHitCount,ULong_t BaseWord);
  void  DecodeWord(ULong_t Code,ULong_t BaseWord,Int_t FirstHalf,ULong_t &Decoded1,ULong_t &Decoded2);
  Int_t fIndex;
  Int_t fHalfStaveModule;
  ClassDef(AliITSDDLRawData,1)
};
    
#endif
