#ifndef ALISTARTRAWREADER_H
#define ALISTARTRAWREADER_H
 
#include <TTask.h>
#include <Riostream.h>
#include "TArrayI.h"

class AliRawReader;
 
class AliSTARTRawReader : public TTask {
  public :

  AliSTARTRawReader() ;

  virtual  ~AliSTARTRawReader();


  Int_t GetPMTId () {return fPMTId;}
  UInt_t UnpackWord(UInt_t PackedWord, Int_t StartBit, Int_t StopBit); // unpack packed words
  Bool_t NextThing(AliRawReader *rawReader); //read next raw digit
  
  TArrayI *TimeTDC1() {return fTimeTDC1;} 
  TArrayI *TimeTDC2() {return fTimeTDC2;} 
  TArrayI *ChargeADC1() {return fChargeADC1;} 
  TArrayI *ChargeADC2() {return fChargeADC2;} 
  virtual void GetTime (TArrayI &o);
  virtual void GetADC (TArrayI &o);
  

protected :

  UInt_t           fData;         // data read for file
 
  AliRawReader*    fRawReader;    // object for reading the raw data
 Int_t fPMTId ;          // PMT number
 TArrayI *fTimeTDC1 ;     //TDC signal
 TArrayI *fChargeADC1 ;   //ADC signal
 TArrayI *fTimeTDC2  ;    //amplified TDC signal
 TArrayI *fChargeADC2 ;   //amplified ADC signal
  
 ClassDef(AliSTARTRawReader, 0) //class for reading START Raw data
};
 
#endif
