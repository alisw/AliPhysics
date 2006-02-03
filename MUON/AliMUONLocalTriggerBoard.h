#ifndef ALIMUONLOCALTRIGGERBOARD_H
#define ALIMUONLOCALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*-- Author: Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONTriggerLut;

class AliMUONLocalTriggerBoard : public AliMUONTriggerBoard 
{
   public: 
      AliMUONLocalTriggerBoard();
      AliMUONLocalTriggerBoard(const char *name, Int_t a);
      virtual ~AliMUONLocalTriggerBoard() {;}

      virtual void     Setbit(Int_t strip, Int_t cathode, Int_t chamber);

      virtual void     Pattern(Option_t *option = ""); // default option displays X then Y bp

      virtual void     Reset();

      virtual void     SetSwitch(Int_t i, Int_t value) {fSwitch[i] = value;}

      virtual UShort_t GetSwitch(Int_t i) {return fSwitch[i];}

      virtual void     SetTC(Bool_t con) {fTC = con;}

      virtual Bool_t   GetTC() {return fTC;}

      virtual void     Module(char *mod);

      virtual void     GetX34(UShort_t *X) {for (Int_t i=0;i<2;i++) X[i] = fXY[0][i+2];}

      virtual void     SetX34(UShort_t *X) {for (Int_t i=0;i<2;i++) fXY[0][i+2] = X[i];}

      virtual void     GetY(UShort_t *Y) {for (Int_t i=0;i<4;i++) Y[i] = fXY[1][i];}

      virtual void     SetY(UShort_t *Y) {for (Int_t i=0;i<4;i++) fXY[1][i] = Y[i];}

      virtual void     GetXY(UShort_t XY[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) XY[i][j] = fXY[i][j];}

      virtual UShort_t GetXY(Int_t i, Int_t j) {return fXY[i][j];} 

      virtual void     SetXY(UShort_t XY[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXY[i][j] = XY[i][j];}

      virtual void     Conf();

      virtual void     Response();

//    inclusive masking
      virtual void     Mask(UShort_t M[2][4]);

//    exclusive masking
      virtual void     Mask(char *in, UShort_t M);

      virtual void     TrigX(Int_t ch1q[16], Int_t ch2q[16], Int_t ch3q[32], Int_t ch4q[32], 
                             Int_t coinc44);
      
      virtual void     Sort2x5(Int_t dev1[6], Int_t dev2[6],
                               Int_t minDev[6], Int_t &dev1GTdev2);
      
      virtual void     TrigY(Int_t y1[16], Int_t y2[16], Int_t y3[16], Int_t y4[16],
                             Int_t y3u[16], Int_t y3d[16], Int_t y4u[16], Int_t y4d[16], 
                             Int_t coinc44);

      virtual void     SetXYU(UShort_t V[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXYU[i][j] = V[i][j];}

      virtual void     SetXYD(UShort_t V[2][4]) {for (Int_t i=0;i<2;i++) for (Int_t j=0;j<4;j++) fXYD[i][j] = V[i][j];}

      virtual void     Scan(Option_t *option = "");

      virtual Int_t    GetI();

      virtual void     LocalTrigger();

      virtual Int_t    Triggered() {return fOutput;}

      virtual Int_t    GetStripX11() {return fStripX11;}

      virtual Int_t    GetStripY11() {return fStripY11;}

      virtual Int_t    GetDev() {return fDev;}
      

   protected:

      virtual void     Resp(Option_t *option); // local trigger info before ("I") and after ("F") LUT

      virtual void     BP(Option_t *option);   // display X/Y bp

   private:

      Int_t    fNumber;

      UShort_t fSwitch[10], fXY[2][4], fXYU[2][4], fXYD[2][4], fMask[2][4];

//    Transverse connector
      Bool_t   fTC;

//    MT1 X position of the valid road 
      Int_t    fStripX11;

//    MT1 Y position of the valid road
      Int_t    fStripY11;

//    Deviation in [0;+30]
      Int_t    fDev;

//    Pt cuts estimated from LUT
      Int_t    fLutLpt[2];        
      Int_t    fLutHpt[2];
      Int_t    fLutApt[2];

//    Outputs of the local logic
      Int_t    fOutput;
      Int_t    fMinDevStrip[5];
      Int_t    fMinDev[5];
      Int_t    fCoordY[5];

      AliMUONTriggerLut *fLUT;
      
      ClassDef(AliMUONLocalTriggerBoard,1) 
};
#endif
