#ifndef ALIT0RAWDATA_H
#define ALIT0RAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts T0 digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliT0;
class AliT0digit;
class AliFstream;
class TFile;
class TMap;
class AliRawDataHeaderSim;
class AliT0RawData : public TObject {

 public:

  AliT0RawData();                                         // default constructor
  AliT0RawData(const AliT0RawData &r);                 // copy constructor
  virtual ~AliT0RawData();                                // destructor
  AliT0RawData &operator=(const AliT0RawData &r);      // ass. op.

   Int_t RawDataT0 (AliT0digit *fDigits); 
  // This method generates the files with the TOF detector data
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;}
  // To set the verbose level
  void  GetDigits(AliT0digit *fDigits);
  //This method formats and stores in buf all the digits of a TOF module

   void  WriteDataHeader(Bool_t dummy, Bool_t compressed);
   void  WriteDRMDataHeader();
   void  WriteTRMDataHeader(UInt_t slotID, Int_t nWords, Int_t positionOfTRMHeader);
   //this method is used to write the data header
   void  WriteTrailer(UInt_t slot, Int_t word1, UInt_t word2, UInt_t word3); 
   void  WriteChainDataHeader(UInt_t chainNumber,UInt_t slotID);
   void  WriteChainDataTrailer(UInt_t chainNumber);
   void  FillTime(Int_t ch, Int_t iTDC, Int_t time);

 //T0 digits arrays


  TArrayI *TimeLED() {return fTimeLED;}
  TArrayI *ADC1() {return fADC1;}
  TArrayI *TimeCFD() {return fTimeCFD;}
  TArrayI *ADC0() {return fADC0;}

 
  
 protected:

  Int_t fVerbose;           //Verbose level (0:no msg, 1:msg, 2:digits in txt files)
  UInt_t fBuffer[512];      // buffer for writing rawdata
  Int_t fIndex;             //number of 32 words to be stored into the output file
  Int_t fEventNumber;       // current event number
  TArrayI * fTimeCFD;       //TDC on the each PMT
  TArrayI *  fADC1;         //QTC (ADC) on the each PMT
  TArrayI * fTimeLED;       // TDC with amplified signal
  TArrayI *  fADC0;         //QTC amplified
  AliFstream* fFile;        //logical name of the I/O file
  UInt_t fDataHeaderPos;    //Data header position
  UInt_t fDRMDataHeaderPos; //Data DRM header position
  UInt_t fTRMDataHeaderPos; //Data TRM header position
  AliT0Parameters *fParam;  // pointer to T0parameters
  TMap fLookUp;             // LookUp Table
  
  ClassDef(AliT0RawData,2)             //  T0 raw data class

};

typedef AliT0RawData AliSTARTRawData; // for backward compatibility

#endif
