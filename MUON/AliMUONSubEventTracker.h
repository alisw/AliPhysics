#ifndef ALIMUONSUBEVENTTRACKER_H
#define ALIMUONSUBEVENTTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONSubEventTracker
/// \brief MUON subevent tracker
///
/// Bus patch structure for tracker raw data
/// each Dsp contains at most 5 bus patch structure
/// Beside the total length and length of the below data
/// the header of the block contains the bus patch id, trigger words 
/// and data structure itself (11bits for manu id, 6 bits for channel id and
/// 12 bits for charge)
///
/// \author Christian Finck

#include <TObject.h>

class AliMUONSubEventTracker : public TObject {

public:
   AliMUONSubEventTracker ();
   virtual ~AliMUONSubEventTracker ();
   AliMUONSubEventTracker(const AliMUONSubEventTracker& rhs);
   AliMUONSubEventTracker& operator=(const AliMUONSubEventTracker& rhs);

   Int_t   GetTotalLength() const {return fTotalLength;}
   Int_t   GetLength()      const {return fLength;}
   Int_t   GetBufSize()     const {return fBufSize;}
   Int_t   GetBusPatchId()  const {return fBusPatchId;}
   Int_t   GetTriggerWord() const {return fTriggerWord;}
   UInt_t  GetData(Int_t n) const;
   UInt_t* GetData()        const {return fData;}

   Char_t   GetParity(Int_t n) const;
   UShort_t GetManuId(Int_t n) const;
   Char_t   GetChannelId(Int_t n) const;
   UShort_t GetCharge(Int_t n) const;

   void    SetBusPatchId(Int_t b)  {fBusPatchId = b;}  
   void    SetTriggerWord(Int_t w) {fTriggerWord = w;}

   void    AddData(UInt_t d);
   void    SetAlloc(Int_t size);

   Bool_t  IsSortable() const {return kTRUE;}
   Int_t   Compare(const TObject *obj) const;

   Int_t   GetHeaderLength() const {return fgkHeaderLength;}

   Int_t*  GetBusPatchHeader() {return &fTotalLength;}

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
