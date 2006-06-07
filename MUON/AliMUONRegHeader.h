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

   //MBZ:3,serialNb:5,Version:8,Id:4,MBZ:4,Out:8
   Char_t   GetSerialNb()  const {return (Char_t)(fWord >> 24) &  0x1F;}
   Char_t   GetVersion()   const {return (Char_t)(fWord >> 16) &  0xFF;}
   Char_t   GetId()        const {return (Char_t)(fWord >> 12) &  0x0F;}
   Char_t   GetOutput()    const {return (Char_t)(fWord)       &  0xFF;}

   void    SetDarcWord(UInt_t w) {fDarcWord = w;}
   void    SetWord(UInt_t w) {fWord = w;}
   void    SetInput(UInt_t in, Int_t n) {fInput[n] = in;}
   
   Int_t   GetHeaderLength() const {return fgkHeaderLength;}
   UInt_t  GetEndOfReg()     const {return fgkEndOfReg;}


   UInt_t* GetHeader() {return &fDarcWord;}

  // scalar methods
   UInt_t  GetL0()      const {return fL0;}
   UInt_t  GetClock()   const {return fClk;}
   const UInt_t* GetScaler()  const {return fScaler;}
   UInt_t  GetHold()    const {return fHold;}

   UInt_t* GetScalers()    {return &fL0;}   

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

   // regional card scalers   
   UInt_t     fL0;         ///< regional L0
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
