#ifndef ALIFITRAWDATA_H
#define ALIFITRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts FIT digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliFIT;
class AliFITDigit;
class AliFstream;
class TFile;
class TMap;
class AliRawDataHeaderSim;
class AliFITRawData : public TObject {

 public:

  AliFITRawData();                                   // default constructor
  AliFITRawData(const AliFITRawData &r);                 // copy constructor
  virtual ~AliFITRawData();                                // destructor
  AliFITRawData &operator=(const AliFITRawData &r);      // ass. op.

  Int_t RawDataFIT (TBranch* branch); 
  // This method generates the files with the FIT detector data
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
  // To set the verbose level
  void  GetDigits();
  //This method formats and stores in buf all the digits of  FIT module

   void  WriteDataHeader(Bool_t dummy, Bool_t compressed);
   void  WriteDRMDataHeader();
   void  WriteTRMDataHeader(UInt_t slotID, Int_t nWords, Int_t positionOfTRMHeader);
   //this method is used to write the data header
   void  WriteTrailer(UInt_t slot, Int_t word1, UInt_t word2, UInt_t word3); 
   void  WriteChainDataHeader(UInt_t chainNumber,UInt_t slotID);
   void  WriteChainDataTrailer(UInt_t chainNumber);
   void  FillTime(Int_t ch, Int_t iTDC, Int_t time);

 
  
 protected:

  TClonesArray *fFITdigitArray;   //Pointer to the FIT digits
  Int_t fVerbose;           //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  UInt_t fBuffer[512];      // buffer for writing rawdata
  Int_t fIndex;             //number of 32 words to be stored into the output file
  Int_t fEventNumber;       // current event number
  Int_t fDataHeaderPos;
   Int_t fAllData[500];
  AliFstream* fFile;        //logical name of the I/O file
 
  ClassDef(AliFITRawData,1)             //  FIT raw data class

};


#endif
