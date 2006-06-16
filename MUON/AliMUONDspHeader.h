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
   Int_t   GetDataKey()        const {return fDataKey;}
   Int_t   GetTotalLength()    const {return fTotalLength;}
   Int_t   GetLength()         const {return fLength;}
   Int_t   GetDspId()          const {return fDspId;}
   Int_t   GetBlkL1ATrigger()  const {return fBlkL1ATrigger;}
   Int_t   GetMiniEventId()    const {return fMiniEventId;}
   Int_t   GetL1ATrigger()     const {return fL1ATrigger;}
   Int_t   GetL1RTrigger()     const {return fL1RTrigger;}
   UInt_t  GetPaddingWord()    const {return fPaddingWord;}
   Int_t   GetErrorWord()      const {return fErrorWord;}

   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}
   UInt_t  GetDefaultDataKey() const {return fgkDefaultDataKey;}
   UInt_t  GetDefaultPaddingWord() const {return fgkDefaultPaddingWord;}

   void    SetDataKey(Int_t d)        {fDataKey = d;}
   void    SetTotalLength(Int_t l)    {fTotalLength = l;}
   void    SetLength(Int_t l)         {fLength = l;}
   void    SetDspId(Int_t d)          {fDspId = d;}  
   void    SetBlkL1ATrigger(Int_t l1) {fBlkL1ATrigger = l1;}
   void    SetMiniEventId(Int_t id)   {fMiniEventId = id;}
   void    SetL1ATrigger(Int_t l1a)   {fL1ATrigger = l1a;}
   void    SetL1RTrigger(Int_t l1r)   {fL1RTrigger = l1r;}
   void    SetPaddingWord(UInt_t w)   {fPaddingWord = w;}
   void    SetErrorWord(Int_t w)      {fErrorWord = w;}

   Int_t*  GetHeader() {return &fDataKey;}

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
   Int_t     fDataKey;          ///< Data key word for FRT header 
   Int_t     fTotalLength;      ///< total length of block structure
   Int_t     fLength;           ///< length of raw data
   Int_t     fDspId;            ///< Dsp id
   Int_t     fBlkL1ATrigger;    ///< L1 accept in Block Structure (CRT)
   Int_t     fMiniEventId;      ///< Mini Event Id in bunch crossing 
   Int_t     fL1ATrigger;       ///< Number of L1 accept in DSP Structure (FRT)
   Int_t     fL1RTrigger;       ///< Number of L1 reject in DSP Structure (FRT)
   Int_t     fPaddingWord;      ///< padding dummy word for 64 bits transfer
   Int_t     fErrorWord;        ///< Error word

   static const Int_t  fgkHeaderLength; ///< header length
   static const UInt_t fgkDefaultDataKey; ///< default data key word for FRT header 
   static const UInt_t fgkDefaultPaddingWord; ///< default padding word value 

   TClonesArray* fBusPatchArray;   ///< array of buspatch structure

   ClassDef(AliMUONDspHeader,2)  // MUON Dsp header for Tracker event
};
#endif
