#ifndef ALIMUONSUBEVENTTRACKER_H
#define ALIMUONSUBEVENTTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONSubEventTracker
/// \brief MUON subevent tracker

#include <TObject.h>

class AliMUONSubEventTracker : public TObject {

public:
   AliMUONSubEventTracker ();
   virtual ~AliMUONSubEventTracker ();
   AliMUONSubEventTracker(const AliMUONSubEventTracker& rhs);
   AliMUONSubEventTracker& operator=(const AliMUONSubEventTracker& rhs);

   Int_t   GetTotalLength() const {return fTotalLength;}
   Int_t   GetLength()      const {return fLength;}
   Int_t   GetBusPatchId()  const {return fBusPatchId;}
   Int_t   GetTriggerWord() const {return fTriggerWord;}
   UInt_t  GetData(Int_t n) const {return fData[n];}
   UInt_t* GetData()        const {return fData;}

   Char_t   GetParity(Int_t n)    {return (Char_t)(fData[n] >> 29) &  0x7;}
   UShort_t GetManuId(Int_t n)    {return (UShort_t)(fData[n] >> 18) &  0x7FF;}
   Char_t   GetChannelId(Int_t n) {return (Char_t)(fData[n] >> 12) & 0x3F;}
   UShort_t GetCharge(Int_t n)    {return (UShort_t)(fData[n] & 0xFFF);}

   void    SetTotalLength(Int_t l) {fTotalLength = l;}
   void    SetLength(Int_t l)      {fLength = l;}
   void    SetBusPatchId(Int_t b)  {fBusPatchId = b;}  
   void    SetTriggerWord(Int_t w) {fTriggerWord = w;}
   void    SetData(UInt_t d, Int_t n)  {fData[n] = d; fLength++;}

   void    AddData(UInt_t d);
   void    SetAlloc(Int_t size);

   Bool_t  IsSortable() const {return kTRUE;}
   Int_t   Compare(const TObject *obj) const;

   Int_t GetHeaderLength() const {return fgkHeaderLength;}

   Int_t* GetAddress() {return &fTotalLength;}

 private:
   Int_t     fTotalLength;   // total length of buspatch structure
   Int_t     fLength;        // length of raw data
   Int_t     fBusPatchId;    // bus patch id
   Int_t     fTriggerWord ;  // counter trigger word

   static const Int_t fgkHeaderLength;   // header length in word

   UInt_t*   fData;          // data 

   Int_t     fBufSize;      // initial size for data array

   void ResizeData(Int_t size = 0);

   ClassDef(AliMUONSubEventTracker,1)  // MUON DDL Tracker
};
#endif
