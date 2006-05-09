#ifndef ALIMUONDSPHEADER_H
#define ALIMUONDSPHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONDspHeader
/// \brief MUON DSP header for tracker event
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONBusStruct;

class AliMUONDspHeader : public TObject {

public:
   AliMUONDspHeader();
   AliMUONDspHeader(const AliMUONDspHeader& event);
   AliMUONDspHeader& operator=(const AliMUONDspHeader& event);

   virtual ~AliMUONDspHeader();

   // DSP header
   Int_t   GetTotalLength()    const {return fTotalLength;}
   Int_t   GetLength()         const {return fLength;}
   Int_t   GetDspId()          const {return fDspId;}
   Int_t   GetTriggerWord(Int_t n) const {return fTriggerWord[n];}
   Int_t   GetEventWord()      const {return fEventWord;}
   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}

   void    SetTotalLength(Int_t l) {fTotalLength = l;}
   void    SetLength(Int_t l)      {fLength = l;}
   void    SetDspId(Int_t d)       {fDspId = d;}  
   void    SetTriggerWord(Int_t w, Int_t n) {fTriggerWord[n] = w;}
   void    SetEventWord(Int_t w)   {fEventWord = w;}

   Int_t*  GetHeader() {return &fTotalLength;}

   void    AddBusPatch(const AliMUONBusStruct& busPatch);

   // get TClonesArray
   TClonesArray*  GetBusPatchArray()  const {return fBusPatchArray;}

   // get entries
   Int_t GetBusPatchEntries()  const {return fBusPatchArray->GetEntriesFast();}

   // get entry
   AliMUONBusStruct* GetBusPatchEntry(Int_t i) const {
     return (AliMUONBusStruct*)fBusPatchArray->At(i);}

   // clear
   void Clear(Option_t* opt);

 private:

   // Dsp header
   Int_t     fTotalLength;      // total length of block structure
   Int_t     fLength;           // length of raw data
   Int_t     fDspId;            // Dsp id ??
   Int_t     fTriggerWord[4];   // counter trigger word ?
   Int_t     fEventWord;        // nb word odd = 1, even = 0
   static const Int_t fgkHeaderLength; // header length

    TClonesArray* fBusPatchArray;   // array of buspatch structure

   ClassDef(AliMUONDspHeader,1)  // MUON Dsp header for Tracker event
};
#endif
