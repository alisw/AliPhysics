#ifndef ALIMUONLOCALSTRUCT_H
#define ALIMUONLOCALSTRUCT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONLocalStruct
/// \brief Rawdata local card structure for trigger
///
//  Author Christian Finck

#include <TObject.h>
#include <TMath.h>
#include <TArrayS.h>


class AliMUONLocalStruct : public TObject{
 
public:
   AliMUONLocalStruct();
   AliMUONLocalStruct(const AliMUONLocalStruct& event);
   AliMUONLocalStruct& operator=(const AliMUONLocalStruct& event);


   virtual ~AliMUONLocalStruct();

   // local board info
            /// Return local data
   UInt_t  GetData(Int_t n) const {return fData[n];}

            /// Return X2
   UShort_t GetX2() const {return (fData[0]     >> 16) &  0xFFFF;}
            /// Return X1
   UShort_t GetX1() const {return (fData[0])           &  0xFFFF;}
            /// Return X4
   UShort_t GetX4() const {return (fData[1] >> 16) &  0xFFFF;}
            /// Return X3
   UShort_t GetX3() const {return (fData[1])       &  0xFFFF;}
           /// return X pattern array
   void    GetXPattern(TArrayS& array) const;
           /// return Y pattern array
   void    GetYPattern(TArrayS& array) const;

            /// Return Y2
   UShort_t GetY2() const {return (fData[2] >> 16) &  0xFFFF;}
            /// Return Y1
   UShort_t GetY1() const {return (fData[2])       &  0xFFFF;}
            /// Return Y4
   UShort_t GetY4() const {return (fData[3] >> 16) &  0xFFFF;}
            /// Return Y3
   UShort_t GetY3() const {return (fData[3])       &  0xFFFF;}

            /// Return Id
   UChar_t  GetId()  const  {return fData[4] >> 19 &  0xF;}
            /// Return Dec
   UChar_t  GetDec() const  {return fData[4] >> 15 &  0xF;}
            /// Return TrigY
   Bool_t   GetTrigY() const {return (fData[4] >> 14 & 0x1);}
            /// Return TriggerY
   Bool_t   GetTriggerY() const {return !(GetTrigY() && GetYPos()==15);}
            /// Return Upos
   UChar_t  GetYPos() const {return fData[4] >> 10 &  0xF;}
            /// Get Sign of X deviation 
   UChar_t  GetSXDev() const {return fData[4] >> 9  &  0x1;}
            /// Get X deviation 
   UChar_t  GetXDev() const {return fData[4] >> 5  &  0xF;}
            /// Return TriggerX
   Bool_t   GetTriggerX() const {return !(GetSXDev() &&
					  !GetXDev() &&
					  GetXPos()==0);}   
            /// Return Xpos
   UChar_t  GetXPos() const {return fData[4]       &  0x1F;}

            /// Return LPT
   Int_t   GetLpt() const {return (GetDec() & 0x3);}
            /// Return HPT
   Int_t   GetHpt() const {return (GetDec() >> 2) & 0x3;}

            /// Set local data
   void    SetData(UInt_t d, Int_t n) {fData[n] = d;}

            /// Return data
   UInt_t* GetData() {return &fData[0];}

 // Scaler methods
            /// Return local L0
   UInt_t  GetL0()      const {return fL0;}
            /// Return local hold (dead time)
   UInt_t  GetHold()    const {return fHold;}
            /// Return local clock
   UInt_t  GetClock()   const {return fClk;}
            /// Return switch
   UShort_t GetSwitch()  const {return (fEOS >> 1) & 0x3FF;}
            /// Return ComptXY
   UChar_t GetComptXY() const {return  fEOS & 1;}

            /// Return XY1
   UShort_t GetXY1(Int_t n) const {return  (n % 2 == 1) ?
       (fScaler[TMath::Nint(Float_t(n/2))] &  0xFFFF) : 
       (fScaler[TMath::Nint(Float_t(n/2))] >> 16) &  0xFFFF;}

            /// Return XY2
   UShort_t GetXY2(Int_t n) const {return  (n % 2 == 1) ?
       (fScaler[8 + TMath::Nint(Float_t(n/2))] &  0xFFFF) : 
       (fScaler[8 + TMath::Nint(Float_t(n/2))] >> 16) &  0xFFFF;}

            /// Return XY3
   UShort_t GetXY3(Int_t n) const {return  (n % 2 == 1) ?
       (fScaler[8*2 + TMath::Nint(Float_t(n/2))] &  0xFFFF) : 
       (fScaler[8*2 + TMath::Nint(Float_t(n/2))] >> 16) &  0xFFFF;}

            /// Return XY4
   UShort_t GetXY4(Int_t n) const {return  (n % 2 == 1) ?
       (fScaler[8*3 + TMath::Nint(Float_t(n/2))] &  0xFFFF) : 
       (fScaler[8*3 + TMath::Nint(Float_t(n/2))] >> 16) &  0xFFFF;}

            /// Return scalers
   UInt_t* GetScalers()  {return &fL0;} 

   // get  length
            /// Return scaler length in word
   Int_t  GetScalerLength()  const {return fgkScalerLength;} 
            /// Return local info length in word
   Int_t  GetLength()        const {return fgkLength;} 
            /// Return end of local info word
   UInt_t GetEndOfLocal()    const {return fgkEndOfLocal;}
            /// Return Word for "empty" slots
   UInt_t GetDisableWord()   const {return fgkDisableWord;}

   // set random numbers to fill variable
   void SetScalersNumbers();

 private:
  
   // local info
   UInt_t    fData[5];  ///< local data
   
   // local card scalers   
   UInt_t     fL0;        ///< local L0
   UInt_t     fHold;      ///< local hold (dead time)
   UInt_t     fClk;       ///< local clock

   UInt_t     fLPtNTrig;  ///< local low Pt no trigger
   UInt_t     fHPtNTrig;  ///< local high Pt no trigger

   UInt_t     fLPtRTrig;  ///< local low Pt right trigger
   UInt_t     fHPtRTrig;  ///< local high Pt right trigger

   UInt_t     fLPtLTrig;  ///< local low Pt left trigger
   UInt_t     fHPtLTrig;  ///< local high Pt left trigger

   UInt_t     fLPtSTrig;  ///< local low Pt straight trigger
   UInt_t     fHPtSTrig;  ///< local high Pt straight trigger

   UInt_t     fScaler[8*4];   ///< local data
   UInt_t     fEOS;           ///< contains switches conf. & flag for reading X (0) or Y (1) in fScaler
   UInt_t     fReset;         ///< reset signal

   static const Int_t  fgkLength;       ///< local info length in word
   static const Int_t  fgkScalerLength; ///< scaler length in word
   static const UInt_t fgkEndOfLocal;   ///< end of local info word
   static const UInt_t fgkDisableWord;  ///< Word for "empty" slots

   ClassDef(AliMUONLocalStruct,3) //Rawdata local card structure for trigger
};
#endif
