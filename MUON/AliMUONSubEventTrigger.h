#ifndef ALIMUONSUBEVENTTRIGGER_H
#define ALIMUONSUBEVENTTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>

class AliMUONSubEventTrigger : public TObject{
 
public:
   AliMUONSubEventTrigger();
   virtual ~AliMUONSubEventTrigger(){;}


   UInt_t  GetRegWord() const {return fRegWord;}
   Int_t   GetRegInput(Int_t n)  const {return fRegInput[n];}
   Int_t   GetRegOutput() const {return fRegOutput;}
   Int_t   GetLocalData(Int_t n)  const {return fLocalData[n];}

   //serialNB:5,DarcId:5,Version:8,MBZ:14
   Char_t   GetSerialNb()  {return (Char_t)(fRegWord >> 27) &  0x1F;}
   Char_t   GetRegId()     {return (Char_t)(fRegWord >> 22) &  0x1F;}
   Char_t   GetVersion()   {return (Char_t)(fRegWord >> 14) &  0xFF;}

   UShort_t GetX2(Int_t n) {return (fLocalData[16*n]     >> 16) &  0xFFFF;}
   UShort_t GetX1(Int_t n) {return (fLocalData[16*n])           &  0xFFFF;}
   UShort_t GetX4(Int_t n) {return (fLocalData[16*n + 1] >> 16) &  0xFFFF;}
   UShort_t GetX3(Int_t n) {return (fLocalData[16*n + 1])       &  0xFFFF;}

   UShort_t GetY2(Int_t n) {return (fLocalData[16*n + 2] >> 16) &  0xFFFF;}
   UShort_t GetY1(Int_t n) {return (fLocalData[16*n + 2])       &  0xFFFF;}
   UShort_t GetY4(Int_t n) {return (fLocalData[16*n + 3] >> 16) &  0xFFFF;}
   UShort_t GetY3(Int_t n) {return (fLocalData[16*n + 3])       &  0xFFFF;}

   Char_t   GetLocalId(Int_t n)  {return fLocalData[16*n + 4] >> 19 &  0xF;}
   Char_t   GetLocalDec(Int_t n) {return fLocalData[16*n + 4] >> 15 &  0xF;}
   Char_t   GetTriggerY(Int_t n) {return fLocalData[16*n + 4] >> 14 &  0x1;}
   Char_t   GetYPos(Int_t n)     {return fLocalData[16*n + 4] >> 10 &  0xF;}
   Char_t   GetXDev(Int_t n)     {return fLocalData[16*n + 4] >> 5  &  0x1F;}
   Char_t   GetXPos(Int_t n)     {return fLocalData[16*n + 4]       &  0x1F;}

   void    SetRegWord(UInt_t w) {fRegWord = w;}
   void    SetRegInput(Int_t in, Int_t n) {fRegInput[n] = in;}
   void    SetRegOutput(Int_t out) { fRegOutput = out;}  
   void    SetLocalData(Int_t d, Int_t n) {fLocalData[n] = d;}
   
   UInt_t* GetAddress() {return &fRegWord;}

 private:
   
   UInt_t    fRegWord;          // first word
   Int_t     fLocalData[16*5];  // local data
   Int_t     fRegInput[2];      // regional input
   Int_t     fRegOutput;        // regional output
 
   ClassDef(AliMUONSubEventTrigger,1)  // MUON Pad Hit
};
#endif
