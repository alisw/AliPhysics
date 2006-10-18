#ifndef ALIMUONBUSSTRUCT_H
#define ALIMUONBUSSTRUCT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONBusStruct
/// \brief MUON buspatch structure for tracker
///
/// \author Christian Finck

#include <TObject.h>

class AliMUONBusStruct : public TObject {

public:
   AliMUONBusStruct ();
   virtual ~AliMUONBusStruct ();
   AliMUONBusStruct(const AliMUONBusStruct& rhs);
   AliMUONBusStruct& operator=(const AliMUONBusStruct& rhs);

   // header
   Int_t   GetDataKey()     const {return fDataKey;}
   Int_t   GetTotalLength() const {return fTotalLength;}
   Int_t   GetLength()      const {return fLength;}
   Int_t   GetBusPatchId()  const {return fBusPatchId;}

   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}
   UInt_t  GetDefaultDataKey() const {return fgkDefaultDataKey;}

   Int_t*  GetHeader() {return &fDataKey;}

   // data
   Int_t   GetBufSize()     const {return fBufSize;}
   UInt_t* GetData()        const {return fData;}
   Int_t   GetBlockId()     const {return fBlkId;}
   Int_t   GetDspId()       const {return fDspId;}

   Char_t   GetParity(Int_t n) const;
   UShort_t GetManuId(Int_t n) const;
   Char_t   GetChannelId(Int_t n) const;
   UShort_t GetCharge(Int_t n) const; 
   UInt_t   GetData(Int_t n) const;

   // header setter
   void    SetDataKey(Int_t d)     {fDataKey = d;}
   void    SetTotalLength(Int_t l) {fTotalLength = l;}
   void    SetLength(Int_t l)      {fLength = l;}
   void    SetBusPatchId(Int_t b)  {fBusPatchId = b;} 
 
   // data 
   void    SetData(UInt_t d, Int_t n)  {fData[n] = d;}
   void    SetBlockId(Int_t b)     {fBlkId = b;}  
   void    SetDspId(Int_t d)       {fDspId = d;}  

   void    AddData(UInt_t d);
   void    SetAlloc(Int_t size);

   // TClonesArray
   Bool_t  IsSortable() const {return kTRUE;}
   Int_t   Compare(const TObject *obj) const;
   void    Clear(Option_t* opt);

 private:
   Int_t     fDataKey;       ///< Data key word for bus patch header 
   Int_t     fTotalLength;   ///< total length of buspatch structure
   Int_t     fLength;        ///< length of raw data
   Int_t     fBusPatchId;    ///< bus patch id

   static const Int_t  fgkHeaderLength;   ///< header length in word
   static const UInt_t fgkDefaultDataKey; ///< default data key word for Bus Patch Header 

   Int_t     fBufSize;       ///< initial size for data array

   UInt_t*   fData;          ///< data 

   Int_t     fDspId;         ///< Dsp number for monitoring
   Int_t     fBlkId;         ///< block numer for monitoring

   void ResizeData(Int_t size = 0);

   ClassDef(AliMUONBusStruct,3)  // MUON DDL Tracker
};
#endif
