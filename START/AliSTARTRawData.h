#ifndef ALISTARTRAWDATA_H
#define ALISTARTRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts START digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TArrayI.h"
class AliSTART;
class AliSTARTdigit;
//class AliRawDataHeader;
class AliSTARTRawData : public TObject {

 public:

  AliSTARTRawData();                                         // default constructor
  AliSTARTRawData(const AliSTARTRawData &r);                 // copy constructor
  virtual ~AliSTARTRawData();                                // destructor
  AliSTARTRawData &operator=(const AliSTARTRawData &r);      // ass. op.

   Int_t RawDataSTART (AliSTARTdigit *fDigits); 
  // This method generates the files with the TOF detector data
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
  // To set the verbose level
  void  GetDigits(AliSTARTdigit *fDigits, UInt_t *buf);
  //This method formats and stores in buf all the digits of a TOF module

  void  PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit);
  //This method stores the value of the variable Word of StopBit-StartBit+1 bits 
  //in BaseWord, starting from the bit StartBit


 //START digits arrays

  AliSTARTdigit *fDigits;  //! The START digits manager

  TArrayI *timeTDC() {return ftimeTDC;} 
  TArrayI *ADC() {return fADC;} 
  
  
 protected:

  Int_t fVerbose;            //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  Int_t fIndex;              //number of 32 words to be stored into the output file

  Int_t fBestTimeRight      ; //smallest time on the right side
  Int_t fBestTimeLeft      ; //smallest time on the left side
  Int_t fTimeDiff     ; //time difference 
  Int_t fMeanTime      ; // average time - ALICE start signal 
  TArrayI *ftimeTDC    ; //array of TDC signal from right side
  TArrayI *fADC    ;   //array of ADC signal from right sida 

  ClassDef(AliSTARTRawData,1)             //  START raw data class

};
#endif
