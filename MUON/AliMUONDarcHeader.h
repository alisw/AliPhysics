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


   UInt_t   GetWord()               const {return fWord;}
   Int_t    GetGlobalInput(Int_t n) const {return fGlobalInput[n];}
   UChar_t  GetGlobalOutput()       const {return (fGlobalOutput & 0xFF);}
   UShort_t GetGlobalConfig()       const {return ((fGlobalOutput >> 16) & 0xFFFF);}

   //MBZ:1, phys trig:1, type:3, ,SerialNb:4,Version:8,VME trig:1, 
   //GlobalFlag:1, CTP trig:1, DAQ:1, Reg pattern:8;

   Bool_t  GetEventType()  const {return (fWord &  0x40000000);}
   UChar_t GetDarcType()   const {return (UChar_t)(fWord >> 27) &  0x7;}
   UChar_t GetSerialNb()   const {return (UChar_t)(fWord >> 23) &  0xF;}
   UChar_t GetVersion()    const {return (UChar_t)(fWord >> 15) &  0xFF;}
   Bool_t  GetVMETrig()    const {return (fWord &  0x800);}
   Bool_t  GetGlobalFlag() const {return (fWord &  0x400);}
   Bool_t  GetCTPTrig()    const {return (fWord &  0x200);}
   Bool_t  GetDAQFlag()    const {return (fWord &  0x100);}
   UChar_t GetRegPattern() const {return (UChar_t)(fWord &  0xFF);}

   void    SetWord(UInt_t w) {fWord = w;}
   void    SetGlobalInput(Int_t in, Int_t n) {fGlobalInput[n] = in;}
   void    SetGlobalOutput(Int_t out) {fGlobalOutput = out;}

   Int_t   GetDarcHeaderLength()   const {return fgkDarcHeaderLength;}
   Int_t   GetGlobalHeaderLength() const {return fgkGlobalHeaderLength;}

   UInt_t* GetHeader() {return &fWord;}
   Int_t*  GetGlobalInput()    {return &fGlobalInput[0];}

  // DARC get methods
   UInt_t  GetDarcL0R()     const {return fDarcL0R;}
   UInt_t  GetDarcL1P()     const {return fDarcL1P;}
   UInt_t  GetDarcL1S()     const {return fDarcL1S;}
   UInt_t  GetDarcL2A()     const {return fDarcL2A;}
   UInt_t  GetDarcL2R()     const {return fDarcL2R;}
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

   UInt_t    fWord;              ///< first word
   Int_t     fGlobalInput[4];    ///< global input
   Int_t     fGlobalOutput;      ///< global ouput

   static const Int_t fgkDarcHeaderLength;   ///< darc header length
   static const Int_t fgkGlobalHeaderLength; ///< global header length


 // global card scalers   
   UInt_t     fGlobalL0;         ///< global L0
   UInt_t     fGlobalClk;        ///< global clock
   UInt_t     fGlobalScaler[6];  ///< global ouput
   UInt_t     fGlobalHold;       ///< global hold (dead time)
   UInt_t     fGlobalSpare;      ///< global spare
   static const Int_t      fgkGlobalScalerLength;  ///< length of global scaler in word

   // DARC Scalers
   UInt_t     fDarcL0R;       ///< DARC L0 received and used
   UInt_t     fDarcL1P;       ///< DARC L1 physics
   UInt_t     fDarcL1S;       ///< DARC L1 software
   UInt_t     fDarcL2A;       ///< DARC L2 accept
   UInt_t     fDarcL2R;       ///< DARC L2 reject
   UInt_t     fDarcClk;       ///< DARC clock
   UInt_t     fDarcHold;      ///< DARC hold (dead time)
   UInt_t     fDarcSpare;     ///< DARC Empty slot (for the moment)

   static const Int_t      fgkDarcScalerLength;  ///< length of DARC scaler in word

   static const UInt_t     fgkEndOfDarc;         ///< end of darc info word
   static const UInt_t     fgkEndOfGlobal;       ///< end of global info word

   TClonesArray* fRegHeaderArray; ///< container for regional header

   ClassDef(AliMUONDarcHeader,1)  // MUON DDL Trigger
};
#endif
