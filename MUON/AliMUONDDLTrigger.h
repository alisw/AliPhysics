#ifndef ALIMUONDDLTRIGGER_H
#define ALIMUONDDLTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>
#include "AliRawDataHeader.h"

class AliMUONDDLTrigger : public TObject {
 
public:
   AliMUONDDLTrigger();
   virtual ~AliMUONDDLTrigger();


   UInt_t  GetDDLWord()            const {return fddlWord;}
   Int_t   GetGlobalInput(Int_t n) const {return fGlobalInput[n];}
   Int_t   GetGlobalOuput()        const {return fGlobalOutput;}
   Int_t   GetEoD()                const {return fEndOfDDL;}  

   //DarcId:2,Version:8,SerialNb:4,EventType:4,MBZ:16;
   Char_t   GetDarcId()    {return (Char_t)(fddlWord >> 30) &  0x10;}
   Char_t   GetVersion()   {return (Char_t)(fddlWord >> 24) &  0xFF;}
   Char_t   GetSerialNb()  {return (Char_t)(fddlWord >> 20) &  0xF;}
   Char_t   GetEventType() {return (Char_t)(fddlWord >> 16) &  0xF;}

   void    SetDDLWord(UInt_t w) {fddlWord = w;}
   void    SetGlobalInput(Int_t in, Int_t n) {fGlobalInput[n] = in;}
   void    SetGlobalOutput(Int_t out) {fGlobalOutput = out;}
   void    SetEoD(Int_t e) {fEndOfDDL = e;}  

   UInt_t* GetEnhancedHeader() {return &fddlWord;}

   AliRawDataHeader GetHeader(){return fHeader;}
   Int_t GetHeaderSize() {return sizeof(AliRawDataHeader)/4;} // in words

 private:

   UInt_t    fddlWord;           // first word
   Int_t     fGlobalInput[4];    // global input
   Int_t     fGlobalOutput;      // global ouput
   Int_t     fEndOfDDL;          // end of DDL

   AliRawDataHeader fHeader;   // header of DDL

   ClassDef(AliMUONDDLTrigger,1)  // MUON DDL Trigger
};
#endif
