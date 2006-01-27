#ifndef ALIMUONSCALEREVENTTRIGGER_H
#define ALIMUONSCALEREVENTTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONScalerEventTrigger
/// \brief MUON scaler event trigger

#include <TObject.h>

class AliMUONScalerEventTrigger : public TObject{
 
public:
   AliMUONScalerEventTrigger();
   virtual ~AliMUONScalerEventTrigger(){;}

   // local get methods
   UInt_t  GetLocalL0()      const {return fLocalL0;}
   UInt_t  GetLocalHold()    const {return fLocalHold;}
   UInt_t  GetLocalClock()   const {return fLocalClk;}
   UChar_t GetLocalSwitch()  const {return (fLocalEOS >> 2) & 0x3FF;}
   UChar_t GetLocalComptXY() const {return  fLocalEOS & 3;}

   UShort_t GetXY1(Int_t n) {return  (n % 2 == 0) ? 
       (fLocalScaler[n] &  0xFFFF) : (fLocalScaler[n] >> 16) &  0xFFFF;}
   UShort_t GetXY2(Int_t n) {return  (n % 2 == 0) ? 
       (fLocalScaler[8 + n] &  0xFFFF) : (fLocalScaler[8 + n] >> 16) &  0xFFFF;}
   UShort_t GetXY3(Int_t n) {return  (n % 2 == 0) ? 
       (fLocalScaler[8*2 + n] &  0xFFFF) : (fLocalScaler[8*2 + n] >> 16) &  0xFFFF;}
   UShort_t GetXY4(Int_t n) {return  (n % 2 == 0) ? 
       (fLocalScaler[8*3 + n] &  0xFFFF) : (fLocalScaler[8*3 + n] >> 16) &  0xFFFF;}

   // regional get methods
   UInt_t  GetRegL0()      const {return fRegL0;}
   UInt_t  GetRegClock()   const {return fRegClk;}
   const UInt_t* GetRegScaler()  const {return fRegScaler;}
   UInt_t  GetRegHold()    const {return fRegHold;}

   // global get methods
   UInt_t  GetGlobalL0()      const {return fGlobalL0;}
   UInt_t  GetGlobalClock()   const {return fGlobalClk;}
   const UInt_t* GetGlobalScaler()  const {return fGlobalScaler;}
   UInt_t  GetGlobalHold()    const {return fGlobalHold;}
   UInt_t  GetGlobalSpare()   const {return fGlobalSpare;}

   // DARC get methods
   UInt_t  GetDarcL0R()     const {return fDarcL0R;}
   UInt_t  GetDarcL0U()     const {return fDarcL0U;}
   UInt_t  GetDarcL0P()     const {return fDarcL0P;}
   UInt_t  GetDarcL0S()     const {return fDarcL0S;}
   UInt_t  GetDarcClock()   const {return fDarcClk;}
   UInt_t  GetDarcHold()    const {return fDarcHold;}
   
   // don't use setting methods but memcpy
   UInt_t* GetLocalScalers()  {return &fLocalL0;} 
   UInt_t* GetRegScalers()    {return &fRegL0;}   
   UInt_t* GetGlobalScalers() {return &fGlobalL0;}
   UInt_t* GetDarcScalers()   {return &fDarcL0R;} 

   // get scaler length
   Int_t GetLocalScalerLength()  const {return fgkLocalScalerLength;} 
   Int_t GetRegScalerLength()    const {return fgkRegScalerLength;}   
   Int_t GetGlobalScalerLength() const {return fgkGlobalScalerLength;}
   Int_t GetDarcScalerLength()   const {return fgkDarcScalerLength;} 

   // set random numbers to fill variable
   void SetNumbers();

 private:

   // local card scalers   
   UInt_t     fLocalL0;        // local L0
   UInt_t     fLocalHold;      // local hold (dead time)
   UInt_t     fLocalClk;       // local clock

   UInt_t     fLocalLPtNTrig;  // local low Pt no trigger
   UInt_t     fLocalHPtNTrig;  // local high Pt no trigger

   UInt_t     fLocalLPtRTrig;  // local low Pt right trigger
   UInt_t     fLocalHPtRTrig;  // local high Pt right trigger

   UInt_t     fLocalLPtLTrig;  // local low Pt left trigger
   UInt_t     fLocalHPtLTrig;  // local high Pt left trigger

   UInt_t     fLocalLPtSTrig;  // local low Pt straight trigger
   UInt_t     fLocalHPtSTrig;  // local high Pt straight trigger

   UInt_t     fLocalScaler[8*4];   // local data
   UInt_t     fLocalEOS;           // contains switches conf. & flag for reading X (0) or Y (1) in fLocalScaler
   UInt_t     fLocalReset;         // reset signal
   static const Int_t      fgkLocalScalerLength;  // length of local scaler in word

   // regional card scalers   
   UInt_t     fRegL0;         // regional L0
   UInt_t     fRegClk;        // regional clock
   UInt_t     fRegScaler[8];  // regional ouput
   UInt_t     fRegHold;       // regional hold (dead time)
   static const Int_t      fgkRegScalerLength;  // length of regional scaler in word

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

   ClassDef(AliMUONScalerEventTrigger,1) 
};
#endif
