#ifndef ALIMUONDARCHEADER_H
#define ALIMUONDARCHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONDarcHeader
/// \brief MUON Darc header for Trigger
///
/// \author Christian Finck

#include <TObject.h>
#include <TClonesArray.h>

class AliMUONRegHeader;

class AliMUONDarcHeader : public TObject {
 
public:
   AliMUONDarcHeader();
   AliMUONDarcHeader(const AliMUONDarcHeader& event);
   AliMUONDarcHeader& operator=(const AliMUONDarcHeader& event);

   virtual ~AliMUONDarcHeader();


   UInt_t  GetWord()               const {return fWord;}
   Int_t   GetGlobalInput(Int_t n) const {return fGlobalInput[n];}
   Int_t   GetGlobalOutput()       const {return fGlobalOutput;}

   //DarcId:4,SerialNb:4,Version:8,EventType:4,GlobalFlag:4,MBZ:8;
   Char_t  GetDarcId()     const {return (Char_t)(fWord >> 28) &  0xF;}
   Char_t  GetSerialNb()   const {return (Char_t)(fWord >> 24) &  0xF;}
   Char_t  GetVersion()    const {return (Char_t)(fWord >> 16) &  0xFF;}
   Char_t  GetEventType()  const {return (Char_t)(fWord >> 12) &  0xF;}
   Char_t  GetGlobalFlag() const {return (Char_t)(fWord >>  8) &  0xF;}

   void    SetWord(UInt_t w) {fWord = w;}
   void    SetGlobalInput(Int_t in, Int_t n) {fGlobalInput[n] = in;}
   void    SetGlobalOutput(Int_t out) {fGlobalOutput = out;}

   Int_t   GetHeaderLength() const {return fgkHeaderLength;}

   UInt_t* GetHeader() {return &fWord;}
   Int_t*  GetGlobalInput()    {return &fGlobalInput[0];}

  // DARC get methods
   UInt_t  GetDarcL0R()     const {return fDarcL0R;}
   UInt_t  GetDarcL0U()     const {return fDarcL0U;}
   UInt_t  GetDarcL0P()     const {return fDarcL0P;}
   UInt_t  GetDarcL0S()     const {return fDarcL0S;}
   UInt_t  GetDarcClock()   const {return fDarcClk;}
   UInt_t  GetDarcHold()    const {return fDarcHold;}
   
   // don't use setting methods but memcpy
   UInt_t* GetGlobalScalers()  {return &fGlobalL0;}
   UInt_t* GetDarcScalers()    {return &fDarcL0R;} 

   // global get methods
   UInt_t  GetGlobalL0()      const {return fGlobalL0;}
   UInt_t  GetGlobalClock()   const {return fGlobalClk;}
   const UInt_t* GetGlobalScaler()  const {return fGlobalScaler;}
   UInt_t  GetGlobalHold()    const {return fGlobalHold;}
   UInt_t  GetGlobalSpare()   const {return fGlobalSpare;}

   Int_t GetGlobalScalerLength() const {return fgkGlobalScalerLength;}
   Int_t GetDarcScalerLength()   const {return fgkDarcScalerLength;} 

   UInt_t GetEndOfDarc()     const {return fgkEndOfDarc;} 
   UInt_t GetEndOfGlobal()   const {return fgkEndOfGlobal;} 

   // set random numbers to fill variable
   void SetScalersNumbers();

  // get TClonesArray
   TClonesArray*  GetRegHeaderArray()  const {return fRegHeaderArray;}

   // get entries
   Int_t GetRegHeaderEntries()  const {return fRegHeaderArray->GetEntriesFast();}

   // get entry
   AliMUONRegHeader* GetRegHeaderEntry(Int_t i) const  {
     return (AliMUONRegHeader*)fRegHeaderArray->At(i);}

   // clear
   void Clear(Option_t* opt);

 private:

   UInt_t    fWord;              // first word
   Int_t     fGlobalInput[4];    // global input
   Int_t     fGlobalOutput;      // global ouput

   static const Int_t fgkHeaderLength; // header length


 // global card scalers   
   UInt_t     fGlobalL0;         // global L0
   UInt_t     fGlobalClk;        // global clock
   UInt_t     fGlobalScaler[6];  // global ouput
   UInt_t     fGlobalHold;       // global hold (dead time)
   UInt_t     fGlobalSpare;      // global spare
   static const Int_t      fgkGlobalScalerLength;  // length of global scaler in word

   // DARC Scalers
   UInt_t     fDarcL0R;       // DARC L0 received
   UInt_t     fDarcL0U;       // DARC L0 used
   UInt_t     fDarcL0P;       // DARC Physical L0
   UInt_t     fDarcL0S;       // DARC Software (checking) L0
   UInt_t     fDarcClk;       // DARC clock
   UInt_t     fDarcHold;      // DARC hold (dead time)
   static const Int_t      fgkDarcScalerLength;  // length of DARC scaler in word

   static const UInt_t     fgkEndOfDarc;         // end of darc info word
   static const UInt_t     fgkEndOfGlobal;       // end of global info word

   TClonesArray* fRegHeaderArray; //container for regional header

   ClassDef(AliMUONDarcHeader,1)  // MUON DDL Trigger
};
#endif
