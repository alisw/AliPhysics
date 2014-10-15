#ifndef ALIMUONBUSSTRUCT_H
#define ALIMUONBUSSTRUCT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONBusStruct
/// \brief MUON buspatch structure for tracker
///
//  Author Christian Finck

#include <TObject.h>

class AliMUONBusStruct : public TObject {

public:
   AliMUONBusStruct ();
   virtual ~AliMUONBusStruct ();
   AliMUONBusStruct(const AliMUONBusStruct& rhs);
   AliMUONBusStruct& operator=(const AliMUONBusStruct& rhs);

           /// Return header length in word
   static Int_t   GetHeaderLength()   {return fgkHeaderLength;}
           /// Return default data key word for Bus Patch Header
   static UInt_t  GetDefaultDataKey() {return fgkDefaultDataKey;}

   // header
           /// Return Data key word for bus patch header
   Int_t   GetDataKey()     const {return fDataKey;}
           /// Return total length of buspatch structure
   Int_t   GetTotalLength() const {return fTotalLength;}
           /// Return length of raw data
   Int_t   GetLength()      const {return fLength;}
           /// Return bus patch id
   Int_t   GetBusPatchId()  const {return fBusPatchId;}

           /// Return header
   Int_t*  GetHeader() {return &fDataKey;}

   // data
           /// Return initial size for data array
   Int_t   GetBufSize()     const {return fBufSize;}
           /// Return data
   UInt_t* GetData()        const {return fData;}
           /// Return block numer for monitoring
   Int_t   GetBlockId()     const {return fBlkId;}
           /// Return Dsp number for monitoring
   Int_t   GetDspId()       const {return fDspId;}

   Char_t   GetParity(Int_t n) const;
   UShort_t GetManuId(Int_t n) const;
   UChar_t   GetChannelId(Int_t n) const;
   UShort_t GetCharge(Int_t n) const; 
   UInt_t   GetData(Int_t n) const;

   // header setter
           /// Set Data key word for bus patch header 
   void    SetDataKey(Int_t d)     {fDataKey = d;}
           /// Set total length of buspatch structure 
   void    SetTotalLength(Int_t l) {fTotalLength = l;}
           /// Set length of raw data 
   void    SetLength(Int_t l)      {fLength = l;}
           /// Set bus patch id 
   void    SetBusPatchId(Int_t b)  {fBusPatchId = b;} 
 
   // data 
           /// Set data 
   void    SetData(UInt_t d, Int_t n)  {fData[n] = d;}
           /// Set block numer for monitoring 
   void    SetBlockId(Int_t b)     {fBlkId = b;}  
           /// Set Dsp number for monitoring 
   void    SetDspId(Int_t d)       {fDspId = d;}  

   void    AddData(UInt_t d);
   void    SetAlloc(Int_t size);

   // TClonesArray
           /// Return true as  Compare() is implemented
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
   static const Int_t  fgkManuNofChannels;///< max number of channels per manu;

   Int_t     fBufSize;       ///< initial size for data array

   UInt_t*   fData;          ///< data 

   Int_t     fDspId;         ///< Dsp number for monitoring
   Int_t     fBlkId;         ///< block numer for monitoring

   void ResizeData(Int_t size = 0);

   void Print(Option_t* opt) const;


   ClassDef(AliMUONBusStruct,3)  // MUON DDL Tracker
};

#endif
