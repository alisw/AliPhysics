#ifndef ALIMUONBLOCKHEADER_H
#define ALIMUONBLOCKHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONBlockHeader
/// \brief MUON block (Crocus CRT)  header for tracker event
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONDspHeader;

class AliMUONBlockHeader : public TObject {

public:
   AliMUONBlockHeader();
   AliMUONBlockHeader(const AliMUONBlockHeader &event);
   AliMUONBlockHeader& operator=(const AliMUONBlockHeader &event);

   virtual ~AliMUONBlockHeader();

   // Block header
   Int_t   GetDataKey()        const {return fDataKey;}
   Int_t   GetTotalLength()    const {return fTotalLength;}
   Int_t   GetLength()         const {return fLength;}
   Int_t   GetDspId()          const {return fDspId;}
   Int_t   GetL0Trigger()      const {return fL0Trigger;}
   Int_t   GetMiniEventId()    const {return fMiniEventId;}
   Int_t   GetEventId1()       const {return fEventId1;}
   Int_t   GetEventId2()       const {return fEventId2;}

   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}
   UInt_t  GetDefaultDataKey() const {return fgkDefaultDataKey;}

   void    SetDataKey(Int_t d)     {fDataKey = d;}
   void    SetTotalLength(Int_t l) {fTotalLength = l;}
   void    SetLength(Int_t l)      {fLength = l;}
   void    SetDspId(Int_t d)       {fDspId = d;}  
   void    SetL0Trigger(Int_t l)   {fL0Trigger = l;}
   void    SetMiniEventId(Int_t e) {fMiniEventId = e;}
   void    SetEventId1(Int_t e)    {fEventId1 = e;}
   void    SetEventId2(Int_t e)    {fEventId2 = e;}

   Int_t* GetHeader() {return &fDataKey;}

   void   AddDspHeader(const AliMUONDspHeader& dspHeader);

   // get TClonesArray
   TClonesArray*      GetDspHeaderArray()  const {return fDspHeaderArray;}

   // get entries
   Int_t              GetDspHeaderEntries() const {return fDspHeaderArray->GetEntriesFast();}

   // get entry
   AliMUONDspHeader*  GetDspHeaderEntry(Int_t i) const {
     return (AliMUONDspHeader*)fDspHeaderArray->At(i);}

   // clear
   void Clear(Option_t* opt);

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


   static const Int_t fgkHeaderLength;   ///< header length in word
   static const UInt_t fgkDefaultDataKey; ///< default data key word for CRT header 

   TClonesArray*  fDspHeaderArray;  ///< array of block header

   ClassDef(AliMUONBlockHeader,2)  // MUON block header for Tracker event
};
#endif
