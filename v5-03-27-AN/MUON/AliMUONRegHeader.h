#ifndef ALIMUONREGHEADER_H
#define ALIMUONREGHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONRegHeader
/// \brief MUON regional header for trigger
///
//  Author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONLocalStruct;

class AliMUONRegHeader : public TObject{
 
public:
   AliMUONRegHeader();
  AliMUONRegHeader(TRootIOCtor* dummy);
   AliMUONRegHeader(const AliMUONRegHeader& event);
   AliMUONRegHeader& operator=(const AliMUONRegHeader& event);

   virtual ~AliMUONRegHeader();

   /// Return darc word
   UInt_t   GetDarcWord()      const {return fDarcWord;}
   /// Return first reg word
   UInt_t   GetWord()          const {return fWord;}
   /// Return regional input
   UInt_t   GetInput(Int_t n)  const {return fInput[n];}
   /// Return L0
   UShort_t GetL0()            const {return fL0   & 0xFFFF;}
   /// Return mask
   UShort_t GetMask()          const {return fMask & 0xFFFF;}

   //word: phys type:1, reset: 6, serialNb:5, Id:4, version: 8, regional output:8
   //true for phys, false for soft
   /// Return RegPhysFlag
   Bool_t    GetRegPhysFlag() const {return (fWord & 0x80000000);} 
   /// Return ResetNb
   UChar_t   GetResetNb()    const {return (UChar_t)(fWord >> 25) &  0x3F;}
   /// Return SerialNb
   UChar_t   GetSerialNb()   const {return (UChar_t)(fWord >> 20) &  0x1F;}
   /// Return Id
   UChar_t   GetId()         const {return (UChar_t)(fWord >> 16) &  0x0F;}
   /// Return Version
   UChar_t   GetVersion()    const {return (UChar_t)(fWord >> 8)  &  0xFF;}
   /// Return Output
   UChar_t   GetOutput()     const {return (UChar_t)(fWord       &  0xFF);}

   //Darc Status: error:10, #fpag:3, MBZ:3, phys type:1, present:1, not_full:1
   // not_empty:1, L2Rej:1, L2Acc:1, L1:1, L0:1, #evt:4, busy:4
   /// Return ErrorBits
   UShort_t GetErrorBits()       const {return (UShort_t)(fDarcWord >> 22) &  0x3FF;}
   /// Return FPGANumber
   UChar_t  GetFPGANumber()      const {return (UChar_t)  (fDarcWord >> 19) &  0x7;}
   /// Return DarcPhysFlag
   Bool_t   GetDarcPhysFlag()    const {return (fDarcWord  &  0x8000);}
   /// Return PresentFlag
   Bool_t   GetPresentFlag()     const {return (fDarcWord  &  0x4000);}
   /// Return RamNotFullFlag
   Bool_t   GetRamNotFullFlag()  const {return (fDarcWord  &  0x2000);}
   /// Return RamNotEmptyFlag
   Bool_t   GetRamNotEmptyFlag() const {return (fDarcWord  &  0x1000);}
   /// Return L2RejStatus
   Bool_t   GetL2RejStatus()     const {return (fDarcWord  &  0x800);}
   /// Return L2AccStatus
   Bool_t   GetL2AccStatus()     const {return (fDarcWord  &  0x400);}
   /// Return L1Status
   Bool_t   GetL1Status()        const {return (fDarcWord  &  0x200);}
   /// Return L0Status
   Bool_t   GetL0Status()        const {return (fDarcWord  &  0x100);}
   /// Return EventInRam
   UChar_t  GetEventInRam()      const {return (UChar_t)  (fDarcWord >> 4)  &  0x4;}
   /// Return Busy
   UChar_t  GetBusy()            const {return (UChar_t)  (fDarcWord)       &  0x4;}

   /// Set darc word
   void    SetDarcWord(UInt_t w) {fDarcWord = w;}
   /// Set first reg word
   void    SetWord(UInt_t w)     {fWord = w;}
   /// Set local mask ("poids faible" 16 bits)
   void    SetMask(UShort_t m)   {fMask = m;}
   /// Set L0
   void    SetL0(UShort_t l)     {fL0   = (l & 0xFFFF);}
   /// Set regional input
   void    SetInput(UInt_t in, Int_t n) {fInput[n] = in;}
 

   /// Return header length in word
   Int_t   GetHeaderLength() const {return fgkHeaderLength;}
   /// Return end of regional info word
   UInt_t  GetEndOfReg()     const {return fgkEndOfReg;}
  /// Return error word
   UInt_t  GetErrorWord()    const {return fgkErrorWord;}

   /// Return header
   UInt_t* GetHeader() {return &fDarcWord;}

  // scalar methods
   /// Return regional clock
   UInt_t  GetClock()   const {return fClk;}
    /// Return regional ouput
   const UInt_t* GetScaler()  const {return fScaler;}
   /// Return regional hold (dead time)
   UInt_t  GetHold()    const {return fHold;}

   /// Return scalers
   UInt_t* GetScalers()    {return &fClk;}   

   /// get scaler length
   Int_t GetScalerLength()  const {return fgkScalerLength;} 

   /// set random numbers to fill variable
   void SetScalersNumbers();

   /// get TClonesArray
   TClonesArray*  GetLocalArray()  const {return fLocalArray;}

   /// get entries
   Int_t GetLocalEntries()  const {return fLocalArray->GetEntriesFast();}

   /// get entry
   AliMUONLocalStruct* GetLocalEntry(Int_t i) const {
     return (AliMUONLocalStruct*)fLocalArray->At(i);}

   // clear
   void Clear(Option_t* opt);

 private:
   
   // regional header
   UInt_t    fDarcWord;      ///< darc word
   UInt_t    fWord;          ///< first reg word
   UInt_t    fInput[2];      ///< regional input
   UShort_t  fMask;          ///< local mask ("poids faible" 16 bits)
   UShort_t  fL0;            ///< L0 counter (16 bits)

   // regional card scalers   
   UInt_t     fClk;        ///< regional clock
   UInt_t     fScaler[8];  ///< regional ouput
   UInt_t     fHold;       ///< regional hold (dead time)

   static const Int_t  fgkScalerLength;  ///< length of regional scaler in word
   static const Int_t  fgkHeaderLength;  ///< header length in word
   static const UInt_t fgkEndOfReg;      ///< end of regional info word
   static const UInt_t fgkErrorWord;     ///< error word when a regional board is missing
   
   TClonesArray* fLocalArray;   ///< array of local structure

   ClassDef(AliMUONRegHeader,4) //MUON regional header for trigger
};
#endif
