#ifndef ALIMUONDDLTRACKER_H
#define ALIMUONDDLTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>
#include "AliRawDataHeader.h"


class AliMUONDDLTracker : public TObject {

public:
   AliMUONDDLTracker();
   virtual ~AliMUONDDLTracker(){;}

   // Block header
   Int_t   GetTotalBlkLength() const {return fTotalBlkLength;}
   Int_t   GetBlkLength()      const {return fBlkLength;}
   Int_t   GetDspId()          const {return fDSPId;}
   Int_t   GetBlkTriggerWord(Int_t n) const {return fBlkTriggerWord[n];}
   const Int_t   GetPadding()  const {return fPadding;}

   void    SetTotalBlkLength(Int_t l) {fTotalBlkLength = l;}
   void    SetBlkLength(Int_t l)      {fBlkLength = l;}
   void    SetDSPId(Int_t d)          {fDSPId = d;}  
   void    SetBlkTriggerWord(Int_t w, Int_t n) {fBlkTriggerWord[n] = w;}

   // DSP header
   Int_t   GetTotalDspLength() const {return fTotalDspLength;}
   Int_t   GetDspLength()      const {return fDspLength;}
   Int_t   GetDspId1()         const {return fDSPId1;}
   Int_t   GetDspTriggerWord(Int_t n) const {return fDspTriggerWord[n];}
   Int_t   GetEventWord()      const {return fEventWord;}

   void    SetTotalDspLength(Int_t l) {fTotalDspLength = l;}
   void    SetDspLength(Int_t l)      {fDspLength = l;}
   void    SetDSPId1(Int_t d)         {fDSPId1 = d;}  
   void    SetDspTriggerWord(Int_t w, Int_t n) {fDspTriggerWord[n] = w;}
   void    SetEventWord(Int_t w)      {fEventWord = w;}

   Int_t* GetBlkHeader() {return &fTotalBlkLength;}
   Int_t* GetDspHeader() {return &fTotalDspLength;}

   AliRawDataHeader GetHeader(){return fHeader;}
   Int_t GetHeaderSize() {return sizeof(AliRawDataHeader)/4;} // in words

   const Int_t   GetEoD()      const {return fEndOfDDL;}  

 private:

   // block header
   Int_t     fTotalBlkLength;    // total length of block structure
   Int_t     fBlkLength;         // length of raw data
   Int_t     fDSPId;             // Dsp id
   Int_t     fBlkTriggerWord[4]; // counter trigger word
   Int_t     fPadding;           // padding dummy word for 64 bits transfer

   // Dsp header
   Int_t     fTotalDspLength;     // total length of block structure
   Int_t     fDspLength;          // length of raw data
   Int_t     fDSPId1;             // Dsp id ??
   Int_t     fDspTriggerWord[4];  // counter trigger word ?
   Int_t     fEventWord;          // nb word odd = 1, even = 0

   static const Int_t fEndOfDDL;  // end of DDL


   AliRawDataHeader fHeader;   // header of DDL
 
 
   ClassDef(AliMUONDDLTracker,1)  // MUON DDL Tracker
};
#endif
