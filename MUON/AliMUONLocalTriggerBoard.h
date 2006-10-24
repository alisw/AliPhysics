#ifndef ALIMUONLOCALTRIGGERBOARD_H
#define ALIMUONLOCALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim 
/// \class AliMUONLocalTriggerBoard
/// \brief Implementation of local trigger board objects
///
/// A local trigger board has as input a bit pattern and returns 
/// the local trigger response after comparison w/ a LUT
///
/// \author Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONTriggerLut;

class AliMUONLocalTriggerBoard : public AliMUONTriggerBoard 
{
   public: 

     enum ESwitch { 
       kX2d = 0,
       kX2m,
       kX2u,
       kOR0,
       kOR1,
       kENY,
       kZeroAllYLSB,
       kZeroDown,
       kZeroMiddle,
       kZeroUp };
       
      AliMUONLocalTriggerBoard();
      AliMUONLocalTriggerBoard(const char *name, Int_t a, AliMUONTriggerLut* lut);
      virtual ~AliMUONLocalTriggerBoard() {;}
      
      Bool_t HasLUT() const { return (fLUT != 0); }
      
      void SetLUT(AliMUONTriggerLut* lut) { fLUT = lut; }
      void SetCoinc44(Int_t coinc44=0) { fCoinc44 = coinc44; }
      
      virtual void     Setbit(Int_t strip, Int_t cathode, Int_t chamber);
      virtual void     SetbitM(Int_t strip, Int_t cathode, Int_t chamber);

      virtual void     Pattern(Option_t *option = "X Y") const; // default option displays X then Y bp

      virtual void     Reset();

      virtual void     SetSwitch(Int_t i, Int_t value) {fSwitch[i] = value;}

      virtual UShort_t GetSwitch(Int_t i) const {return fSwitch[i];}

      virtual void     SetTC(Bool_t con) {fTC = con;}

      virtual Bool_t   GetTC() const {return fTC;}

      virtual void     SetNumber(Int_t nb) {fNumber = nb;}

      virtual Int_t    GetNumber() const {return fNumber;}

      virtual void     Module(char *mod);

      virtual void     GetX34(UShort_t *X) const {for (Int_t i=0;i<2;i++) X[i] = fXY[0][i+2];}

      virtual void     SetX34(UShort_t *X) {for (Int_t i=0;i<2;i++) fXY[0][i+2] = X[i];}

      virtual void     GetY(UShort_t *Y) const {for (Int_t i=0;i<4;i++) Y[i] = fXY[1][i];}

      virtual void     SetY(UShort_t *Y) {for (Int_t i=0;i<4;i++) fXY[1][i] = Y[i];}

      virtual void     GetXY(UShort_t XY[2][4]) const {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) XY[i][j] = fXY[i][j];}

      virtual UShort_t GetXY(Int_t i, Int_t j) const {return fXY[i][j];} 

      virtual void     SetXY(UShort_t XY[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXY[i][j] = XY[i][j];}

      virtual void     Conf() const;

      virtual void     Response();

      virtual void     Mask(Int_t index, UShort_t value);

      virtual void     TrigX(Int_t ch1q[16], Int_t ch2q[16], Int_t ch3q[32], Int_t ch4q[32]);
      
      virtual void     Sort2x5(Int_t dev1[6], Int_t dev2[6],
                               Int_t minDev[6], Int_t &dev1GTdev2);
      
      virtual void     TrigY(Int_t y1[16], Int_t y2[16], Int_t y3[16], Int_t y4[16],
                             Int_t y3u[16], Int_t y3d[16], Int_t y4u[16], Int_t y4d[16]);

      virtual void     SetXYU(UShort_t V[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXYU[i][j] = V[i][j];}

      virtual void     SetXYD(UShort_t V[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXYD[i][j] = V[i][j];}

      virtual void     Scan(Option_t *option = "") const;

      virtual Int_t    GetI() const;

      virtual void     LocalTrigger();

      virtual Int_t    Triggered() const {return fOutput;}

      virtual Int_t    GetStripX11() const {return fStripX11;}

      virtual Int_t    GetStripY11() const {return fStripY11;}

      virtual Int_t    GetDev() const {return fDev;}
      
      virtual void     SetCrate(TString crate) {fCrate = crate;}

      virtual TString  GetCrate() const {return fCrate;}

   protected:
      AliMUONLocalTriggerBoard(const AliMUONLocalTriggerBoard& right);
      AliMUONLocalTriggerBoard&  operator = (const AliMUONLocalTriggerBoard& right);
      
      static const Int_t fgkCircuitId[234]; // old numbering (to be removed)

      virtual void     Resp(Option_t *option) const; // local trigger info before ("I") and after ("F") LUT

      virtual void     BP(Option_t *option) const;   // display X/Y bp

   private:

      Int_t    fNumber;           ///< Board number

      TString  fCrate;            ///< Crate name

      UShort_t fSwitch[10];       ///< Switch
      UShort_t fXY[2][4];         ///< Bit pattern
      UShort_t fXYU[2][4];        ///< Bit pattern UP
      UShort_t fXYD[2][4];        ///< Bit pattern DOWN
      UShort_t fMask[2][4];       ///< Mask

      Bool_t   fTC;               ///< Transverse connector

      Int_t    fStripX11;         ///< MT1 X position of the valid road 

      Int_t    fStripY11;         ///< MT1 Y position of the valid road

      Int_t    fDev;              ///< Deviation in [0;+30]

      Int_t    fLutLpt[2];        ///< Low Pt cuts estimated from LUT
      Int_t    fLutHpt[2];        ///< High Pt cuts estimated from LUT

//    Outputs of the local logic
      Int_t    fOutput;           ///< Board has triggered
      Int_t    fMinDevStrip[5];   ///< X (from algo)
      Int_t    fMinDev[5];        ///< Dev (from algo)
      Int_t    fCoordY[5];        ///< Y (from algo)
      
      AliMUONTriggerLut *fLUT;    //!< Pointer to trigger LUT, that we do not own.
      Int_t    fCoinc44;          ///< Coinc 44 (0/1 = coinc34/coinc44)
      
      ClassDef(AliMUONLocalTriggerBoard,1) 
};
#endif
