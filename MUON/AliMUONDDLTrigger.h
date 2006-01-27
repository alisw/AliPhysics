#ifndef ALIMUONDDLTRIGGER_H
#define ALIMUONDDLTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONDDLTrigger
/// \brief MUON DDL Trigger

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

   //DarcId:4,SerialNb:4,Version:8,EventType:4,GlobalFlag:4,MBZ:8;
   Char_t   GetDarcId()     {return (Char_t)(fddlWord >> 28) &  0xF;}
   Char_t   GetSerialNb()   {return (Char_t)(fddlWord >> 24) &  0xF;}
   Char_t   GetVersion()    {return (Char_t)(fddlWord >> 16) &  0xFF;}
   Char_t   GetEventType()  {return (Char_t)(fddlWord >> 12) &  0xF;}
   Char_t   GetGlobalFlag() {return (Char_t)(fddlWord >>  8) &  0xF;}

   void    SetDDLWord(UInt_t w) {fddlWord = w;}
   void    SetGlobalInput(Int_t in, Int_t n) {fGlobalInput[n] = in;}
   void    SetGlobalOutput(Int_t out) {fGlobalOutput = out;}
   void    SetEoD(Int_t e) {fEndOfDDL = e;}  

   Int_t GetHeaderLength() const {return fgkHeaderLength;}


   UInt_t* GetEnhancedHeader() {return &fddlWord;}
   Int_t*  GetGlobalInput()    {return &fGlobalInput[0];}


   AliRawDataHeader GetHeader(){return fHeader;}
   Int_t GetHeaderSize() {return sizeof(AliRawDataHeader)/4;} // in words

 private:

   UInt_t    fddlWord;           // first word
   Int_t     fGlobalInput[4];    // global input
   Int_t     fGlobalOutput;      // global ouput

   static const Int_t fgkHeaderLength; // header length

   Int_t     fEndOfDDL;          // end of DDL

   AliRawDataHeader fHeader;   // header of DDL

   ClassDef(AliMUONDDLTrigger,1)  // MUON DDL Trigger
};
#endif
