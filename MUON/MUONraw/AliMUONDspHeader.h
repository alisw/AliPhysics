#ifndef ALIMUONDSPHEADER_H
#define ALIMUONDSPHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONDspHeader
/// \brief MUON DSP header for tracker event
///
//  Author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONBusStruct;

class AliMUONDspHeader : public TObject {

public:
   AliMUONDspHeader();
  AliMUONDspHeader(TRootIOCtor* dummy);
   AliMUONDspHeader(const AliMUONDspHeader& event);
   AliMUONDspHeader& operator=(const AliMUONDspHeader& event);

   virtual ~AliMUONDspHeader();

   // DSP header
           /// Return Data key word for FRT header
   Int_t   GetDataKey()        const {return fDataKey;}
           /// Return total length of block structure
   Int_t   GetTotalLength()    const {return fTotalLength;}
           /// Return length of raw data
   Int_t   GetLength()         const {return fLength;}
           /// Return Dsp id
   Int_t   GetDspId()          const {return fDspId;}
           /// Return L1 accept in Block Structure (CRT)
   Int_t   GetBlkL1ATrigger()  const {return fBlkL1ATrigger;}
           /// Return Mini Event Id in bunch crossing
   Int_t   GetMiniEventId()    const {return fMiniEventId;}
           /// Return Number of L1 accept in DSP Structure (FRT)
   Int_t   GetL1ATrigger()     const {return fL1ATrigger;}
           /// Return Number of L1 reject in DSP Structure (FRT)
   Int_t   GetL1RTrigger()     const {return fL1RTrigger;}
           /// Return padding dummy word for 64 bits transfer
   UInt_t  GetPaddingWord()    const {return fPaddingWord;}
           /// Return Error word
   Int_t   GetErrorWord()      const {return fErrorWord;}

           /// Return header length
   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}
           /// Return default data key word for FRT header
   UInt_t  GetDefaultDataKey() const {return fgkDefaultDataKey;}
           /// Return default padding word value
   UInt_t  GetDefaultPaddingWord() const {return fgkDefaultPaddingWord;}

           /// Set Data key word for FRT header
   void    SetDataKey(Int_t d)        {fDataKey = d;}
           /// Set total length of block structure
   void    SetTotalLength(Int_t l)    {fTotalLength = l;}
           /// Set length of raw data
   void    SetLength(Int_t l)         {fLength = l;}
           /// Set Dsp id
   void    SetDspId(Int_t d)          {fDspId = d;}  
           /// Set L1 accept in Block Structure (CRT)
   void    SetBlkL1ATrigger(Int_t l1) {fBlkL1ATrigger = l1;}
           /// Set Mini Event Id in bunch crossing
   void    SetMiniEventId(Int_t id)   {fMiniEventId = id;}
           /// Set Number of L1 accept in DSP Structure (FRT)
   void    SetL1ATrigger(Int_t l1a)   {fL1ATrigger = l1a;}
           /// Set Number of L1 reject in DSP Structure (FRT)
   void    SetL1RTrigger(Int_t l1r)   {fL1RTrigger = l1r;}
           /// Set padding dummy word for 64 bits transfer
   void    SetPaddingWord(UInt_t w)   {fPaddingWord = w;}
           /// Set Error word
   void    SetErrorWord(Int_t w)      {fErrorWord = w;}

           /// Return header
   Int_t*  GetHeader() {return &fDataKey;}

   void    AddBusPatch(const AliMUONBusStruct& busPatch);

   /// get TClonesArray
   TClonesArray*  GetBusPatchArray()  const {return fBusPatchArray;}

   /// get entries
   Int_t GetBusPatchEntries()  const {return fBusPatchArray->GetEntriesFast();}

   /// get entry
   AliMUONBusStruct* GetBusPatchEntry(Int_t i) const {
     return (AliMUONBusStruct*)fBusPatchArray->At(i);}

   // clear
   void Clear(Option_t* opt);

   // print out
   void Print(Option_t* /*opt*/) const;

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
