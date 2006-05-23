#ifndef ALIMUONBLOCKHEADER_H
#define ALIMUONBLOCKHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONBlockHeader
/// \brief MUON block header for tracker event
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
   Int_t   GetTotalLength()    const {return fTotalLength;}
   Int_t   GetLength()         const {return fLength;}
   Int_t   GetDspId()          const {return fDspId;}
   Int_t   GetTriggerWord(Int_t n) const {return fTriggerWord[n];}
   Int_t   GetPadding()        const {return fPadding;}
   Int_t   GetHeaderLength()   const {return fgkHeaderLength;}

   void    SetTotalLength(Int_t l) {fTotalLength = l;}
   void    SetLength(Int_t l)      {fLength = l;}
   void    SetDspId(Int_t d)       {fDspId = d;}  
   void    SetTriggerWord(Int_t w, Int_t n) {fTriggerWord[n] = w;}

   Int_t* GetHeader() {return &fTotalLength;}

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
   Int_t     fTotalLength;    ///< total length of block structure (w/o padding word)
   Int_t     fLength;         ///< length of raw data
   Int_t     fDspId;          ///< Dsp id
   Int_t     fTriggerWord[4]; ///< counter trigger word
   Int_t     fPadding;        ///< padding dummy word for 64 bits transfer

   static const Int_t fgkHeaderLength; ///< header length in word
 
   TClonesArray*  fDspHeaderArray;  ///< array of block header

   ClassDef(AliMUONBlockHeader,1)  // MUON block header for Tracker event
};
#endif
