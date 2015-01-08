#ifndef AliPHOSCPVParam_h
#define AliPHOSCPVParam_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class provides a set of static methods to convert absolute number of pad to pair (X,Y) 
// and vice versa
// and some other
// Author - Mikhail Stolpovskiy, IHEP Protvino (2013)

#include <TNamed.h>        //base class
#include "AliPHOSCpv3GConnection.h"

class AliPHOSCpvParam :public TNamed  
{
 public:
  enum EChamberData{kMinCh=0,kMaxCh=0};      //Segmenation. CPV has only one chamber
  enum EPadxData{kPadPcX=128,kMinPx=0,kMaxPx=127};   //Segmentation structure along x
  enum EPadyData{kPadPcY=60 ,kMinPy=0,kMaxPy=59 };   //Segmentation structure along y 
  enum {
    kNRows       = 16,    // Number of rows (column controlers)
    kN3GAdd      = 10,    // Number of 3GASSIPLEXs in a row
    kNPadAdd     = 48,    // Number of pad row
    kNRowsPerSegment = 8, // Number of rows per segment
    kNDDL = 5,            // Number of already installed modules (ddls)
    kNModules = 5         // Number of modules (equals to the number of PHOS modules)
   };

  // x <=> phi
  // y <=> Z
  // But x-y is a local module axes (Int_t)

  static Bool_t IsValidAbs(Int_t abs);
  static Bool_t DecodeRawWord(Int_t ddl,Int_t rWord, Int_t & abs, Int_t & q, Int_t & eType);
  static Int_t Abs  (Int_t ddl,Int_t columnCtrl,Int_t gassiplex3,Int_t pad); // abs pad
  static Int_t A2DDL(Int_t abs) ;           // abs pad -> ddl
  static Int_t A2Mod(Int_t abs) ;           // abs -> number of module
  static Int_t DDL2Mod(Int_t ddl);
  static Int_t Mod2DDL(Int_t mod);
  static Int_t A2CC (Int_t abs) ;           // abs pad -> column controler
  static Int_t A23G (Int_t abs) ;           // abs pad -> number of 3gassiplex card
  static Int_t A2Pad(Int_t abs) ;           // abs pad -> number of pad in 3gassiplex
  static Int_t A2X  (Int_t abs) ;           // abs pad -> pad X
  static Int_t A2Y  (Int_t abs) ;           // abs pad -> pad Y
  static Int_t XY2A (Int_t ddl, Int_t x, Int_t y) ;    // pad X,Y -> abs pad  
  static Int_t X2CC (Int_t x) ;             // pad X -> number of column controller
  static Int_t Y23G (Int_t y) ;             // pad Y -> number of 3gassiplex card
  static Int_t XY2Pad (Int_t x, Int_t y) ;  // pad X,Y -> number of pad in 3gassiplex
  static Bool_t GetLimOfCConX( Int_t cc, Int_t &xmin, Int_t &xmax); // returns limits on X for column controler cc
  static Bool_t GetLimOf3GonY( Int_t g3, Int_t &ymin, Int_t &ymax); // returns limits on Y for 3gassiplex g3

  static Int_t A2fId(Int_t abs) ;           // returns number of channel with common PHOS+CPV numeration 

 private:
  // connection of channels of 3gassiplex to pads 
  static AliPHOSCpv3GConnection fConnection;
  ClassDef(AliPHOSCpvParam,1);           //CPV main parameters class
};
#endif
