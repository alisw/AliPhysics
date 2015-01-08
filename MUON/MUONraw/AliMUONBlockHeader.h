#ifndef ALIMUONBLOCKHEADER_H
#define ALIMUONBLOCKHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONBlockHeader
/// \brief MUON block (Crocus CRT)  header for tracker event
///
//  Author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONDspHeader;

class AliMUONBlockHeader : public TObject {

public:
   AliMUONBlockHeader();
   AliMUONBlockHeader(TRootIOCtor* dummy);
   AliMUONBlockHeader(const AliMUONBlockHeader &event);
   AliMUONBlockHeader& operator=(const AliMUONBlockHeader &event);

   virtual ~AliMUONBlockHeader();

   //
   // Block header
   //
           /// Return data key word for CRT header
   Int_t   GetDataKey()        const {return fDataKey;}
           /// Return total length of block structure (w/o padding word)
   Int_t   GetTotalLength()    const {return fTotalLength;}
           /// Return length of raw data
   Int_t   GetLength()         const {return fLength;}
           /// Return Dsp id
   Int_t   GetDspId()          const {return fDspId;}
           /// Return L0 trigger word
   Int_t   GetL0Trigger()      const {return fL0Trigger;}
           /// Return Bunch Crossing for mini-event id (see TDR chapter 8)
   Int_t   GetMiniEventId()    const {return fMiniEventId;}
           /// Return Event Id in bunch crossing
   Int_t   GetEventId1()       const {return fEventId1;}
           /// Return Event Id in orbit number
   Int_t   GetEventId2()       const {return fEventId2;}

           /// Return header length in word
   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}
           /// Return default data key word for CRT header
   UInt_t  GetDefaultDataKey() const {return fgkDefaultDataKey;}
         /// Return  data key end word for CRT header
   UInt_t  GetDdlDataKey() const {return fgkDdlDataKey;}

           /// Set data key word for CRT header
   void    SetDataKey(Int_t d)     {fDataKey = d;}
           /// Set total length of block structure (w/o padding word)
   void    SetTotalLength(Int_t l) {fTotalLength = l;}
           /// Set length of raw data
   void    SetLength(Int_t l)      {fLength = l;}
           /// Set Dsp id
   void    SetDspId(Int_t d)       {fDspId = d;}  
           /// Set L0 trigger word
   void    SetL0Trigger(Int_t l)   {fL0Trigger = l;}
           /// Set Bunch Crossing for mini-event id (see TDR chapter 8)
   void    SetMiniEventId(Int_t e) {fMiniEventId = e;}
           /// Set Event Id in bunch crossing
   void    SetEventId1(Int_t e)    {fEventId1 = e;}
           /// Set Event Id in orbit number
   void    SetEventId2(Int_t e)    {fEventId2 = e;}
   
   
   /// Return header
   Int_t* GetHeader() {return &fDataKey;}

   void   AddDspHeader(const AliMUONDspHeader& dspHeader);

   /// get TClonesArray
   TClonesArray*      GetDspHeaderArray()  const {return fDspHeaderArray;}

   /// get entries
   Int_t              GetDspHeaderEntries() const {return fDspHeaderArray->GetEntriesFast();}

   /// get entry
   AliMUONDspHeader*  GetDspHeaderEntry(Int_t i) const {
     return (AliMUONDspHeader*)fDspHeaderArray->At(i);}

   // clear
   void Clear(Option_t* opt);

   // print out
   void Print(Option_t* /*opt*/) const;

 private:

   // block header
   Int_t     fDataKey;        ///< Data key word for CRT header 
   Int_t     fTotalLength;    ///< total length of block structure (w/o padding word)
   Int_t     fLength;         ///< length of raw data
   Int_t     fDspId;          ///< Dsp id
   Int_t     fL0Trigger;      ///< L0 trigger word
   Int_t     fMiniEventId;    ///< Bunch Crossing for mini-event id (see TDR chapter 8)
   Int_t     fEventId1;       ///< Event Id in bunch crossing
   Int_t     fEventId2;       ///< Event Id in orbit number


   static const Int_t  fgkHeaderLength;    ///< header length in word
   static const UInt_t fgkDefaultDataKey;  ///< default data key word for CRT header 
   static const UInt_t fgkDdlDataKey;      ///< data key end word for CRT header 

   TClonesArray*  fDspHeaderArray;  ///< array of block header

   ClassDef(AliMUONBlockHeader,3)  // MUON block header for Tracker event
};
#endif
