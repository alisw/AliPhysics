#ifndef ALIMUONDDLTRACKER_H
#define ALIMUONDDLTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>
static const Int_t  bufsize = 1024;

class AliMUONDDLTracker : public TObject {

public:
   AliMUONDDLTracker();
   virtual ~AliMUONDDLTracker(){;}

   UInt_t  GetRawData(Int_t n) const {return fData[n];}
   Int_t   GetLength()  const {return fLength;}
   Int_t   GetBusPatchId()  const {return fBusPatchId;}
   Int_t   GetEoD()     const {return fEndOfDDL;}  

   Char_t   GetParity(Int_t n)    {return (Char_t)(fData[n] >> 29) &  0x7;}
   UShort_t GetManuId(Int_t n)    {return (UShort_t)(fData[n] >> 18) &  0x7FF;}
   Char_t   GetChannelId(Int_t n) {return (Char_t)(fData[n] >> 12) & 0x3F;}
   UShort_t GetCharge(Int_t n)    {return (UShort_t)(fData[n] & 0xFFF);}

   void    SetRawData(UInt_t w) {fData[fLength++] = w;}
   void    SetLength(Int_t l) {fLength = l;}
   void    SetBusPatchId(Int_t b) {fBusPatchId = b;}  
   void    SetEoD(Int_t e) {fEndOfDDL = e;}  

   Int_t* GetAddress() {return &fLength;}

 private:
   
   Int_t     fLength;           // length of data
   Int_t     fBusPatchId;       // bus patch id
   UInt_t    fData[bufsize];    // data 
   Int_t     fEndOfDDL  ;       // end of DDL

   ClassDef(AliMUONDDLTracker,1)  // MUON DDL Tracker
};
#endif
