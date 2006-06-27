#ifndef ALIMUONREGHEADER_H
#define ALIMUONREGHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONRegHeader
/// \brief MUON regional header for trigger
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONLocalStruct;

class AliMUONRegHeader : public TObject{
 
public:
   AliMUONRegHeader();
   AliMUONRegHeader(const AliMUONRegHeader& event);
   AliMUONRegHeader& operator=(const AliMUONRegHeader& event);

   virtual ~AliMUONRegHeader();

   UInt_t  GetDarcWord()      const {return fDarcWord;}
   UInt_t  GetWord()          const {return fWord;}
   UInt_t  GetInput(Int_t n)  const {return fInput[n];}
   UInt_t  GetL0()            const {return ((fL0Mask >> 16) & 0xFFFF);}
   UInt_t  GetMask()          const {return (fL0Mask & 0xFFFF);}

   //word: phys type:1, reset: 6, serialNb:5, Id:4, version: 8, regional output:8
   //true for phys, false for soft
   Bool_t   GetRegPhysFlag() const {return (fWord & 0x80000000);} 
   Char_t   GetResetNb()     const {return (Char_t)(fWord >> 25) &  0x20;}
   Char_t   GetSerialNb()    const {return (Char_t)(fWord >> 20) &  0x1F;}
   Char_t   GetId()          const {return (Char_t)(fWord >> 16) &  0x0F;}
   Char_t   GetVersion()     const {return (Char_t)(fWord >> 8)  &  0xFF;}
   Char_t   GetOutput()      const {return (Char_t)(fWord)       &  0xFF;}

   //Darc Status: error:10, #fpag:3, MBZ:3, phys type:1, present:1, not_full:1
   // not_empty:1, L2Rej:1, L2Acc:1, L1:1, L0:1, #evt:4, busy:4
   UShort_t GetErrorBits()       const {return (UShort_t)(fDarcWord >> 21) &  0x3FF;}
   Char_t   GetFPGANumber()      const {return (Char_t)  (fDarcWord >> 18) &  0x7;}
   Bool_t   GetDarcPhysFlag()    const {return (fDarcWord  &  0x1000);}
   Bool_t   GetPresentFlag()     const {return (fDarcWord  &  0x800);}
   Bool_t   GetRamNotFullFlag()  const {return (fDarcWord  &  0x400);}
   Bool_t   GetRamNotEmptyFlag() const {return (fDarcWord  &  0x200);}
   Bool_t   GetL2RejStatus()     const {return (fDarcWord  &  0x100);}
   Bool_t   GetL2AccStatus()     const {return (fDarcWord  &  0x80);}
   Bool_t   GetL1Status()        const {return (fDarcWord  &  0x40);}
   Bool_t   GetL0Status()        const {return (fDarcWord  &  0x20);}
   Char_t   GetEventInRam()      const {return (Char_t)  (fDarcWord >> 4)  &  0x4;}
   Char_t   GetBusy()            const {return (Char_t)  (fDarcWord)       &  0x4;}

   void    SetDarcWord(UInt_t w) {fDarcWord = w;}
   void    SetWord(UInt_t w) {fWord = w;}
   void    SetInput(UInt_t in, Int_t n) {fInput[n] = in;}
   
   Int_t   GetHeaderLength() const {return fgkHeaderLength;}
   UInt_t  GetEndOfReg()     const {return fgkEndOfReg;}


   UInt_t* GetHeader() {return &fDarcWord;}

  // scalar methods
   UInt_t  GetClock()   const {return fClk;}
   const UInt_t* GetScaler()  const {return fScaler;}
   UInt_t  GetHold()    const {return fHold;}

   UInt_t* GetScalers()    {return &fClk;}   

   // get scaler length
   Int_t GetScalerLength()  const {return fgkScalerLength;} 

  // set random numbers to fill variable
   void SetScalersNumbers();

  // get TClonesArray
   TClonesArray*  GetLocalArray()  const {return fLocalArray;}

   // get entries
   Int_t GetLocalEntries()  const {return fLocalArray->GetEntriesFast();}

   // get entry
   AliMUONLocalStruct* GetLocalEntry(Int_t i) const {
     return (AliMUONLocalStruct*)fLocalArray->At(i);}

   // clear
   void Clear(Option_t* opt);

 private:
   
   // regional header
   UInt_t    fDarcWord;      ///< darc word
   UInt_t    fWord;          ///< first reg word
   UInt_t    fInput[2];      ///< regional input
   UInt_t    fL0Mask;        ///< L0 counter (16 bits) and local mask (16 bits)

   // regional card scalers   
   UInt_t     fClk;        ///< regional clock
   UInt_t     fScaler[8];  ///< regional ouput
   UInt_t     fHold;       ///< regional hold (dead time)

   static const Int_t  fgkScalerLength;  ///< length of regional scaler in word
   static const Int_t  fgkHeaderLength;  ///< header length in word
   static const UInt_t fgkEndOfReg;      ///< end of regional info word

   TClonesArray* fLocalArray;   ///< array of local structure

   ClassDef(AliMUONRegHeader,3)
};
#endif
